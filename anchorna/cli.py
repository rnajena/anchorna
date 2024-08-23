# (C) 2024, Tom Eulenfeld, MIT license

import argparse
from copy import deepcopy
from importlib.resources import files
import logging
import logging.config
from pathlib import Path
import os
import subprocess
import sys
import tempfile
import time
import tomllib
from warnings import warn

from sugar import read
from anchorna.core import combine, cutout, find_my_anchors
from anchorna.io import jalview_features, read_anchors, load_selected_anchors
from anchorna.util import _apply_mode


EXAMPLE_TOML_CONFIG = """### Configuration for AnchoRNA in TOML format

fname = "pesti_example.gff" # sequence file with CDS definition
w = 5                      # word length
gseqid = "KC533775"        # ID of the guiding sequence, set to "none" to take first sequence in file
search_range = 100         # number of letters to search left and right from guiding fluke
scoring = "blosum62"       # substition matrix used for scoring,
                           # for possible values see https://rnajena-sugar.readthedocs.io/en/latest/search.html?q=submat
                           # or use your own file
score_add_word = 22        # score to add a word to word set, flukes with lower score are marked as "poor"
thr_quota_add_anchor = 1   # discard anchor if quota is not met (1=100%)
thr_score_add_anchor = 22  # score threshold for quota
score_use_fluke = 22       # score to use fluke in export, view, cutout
aggressive_remove = true   # turn aggressive remove mode on/off

removed_anchors_path = "removed_anchors_{}.gff"  # turn off with "none", alternatively delete this line
logfile = "anchorna.log"   # turn off with "none", alternatively delete this line
"""


LOGGING_DEFAULT_CONFIG = {
    'version': 1,
    'disable_existing_loggers': False,
    'capture_warnings': True,
    'formatters': {
        'file': {
            'format': ('%(asctime)s  %(levelname)s:%(name)s:%(message)s'),
            'datefmt': '%Y-%m-%d %H:%M:%S'
        },
        'console': {
            'format': '%(levelname)s:%(message)s'
        },
    },
    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
            'formatter': 'console',
            'level': 'DEBUG',
        },
        'file': {
            'class': 'logging.FileHandler',
            'formatter': 'file',
            'level': 'DEBUG',
            'filename': None,
        },
    },
    'loggers': {
        'anchorna': {
            'handlers': ['console', 'file'],
            'level': 'DEBUG',
            'propagate': False,
        },
        'py.warnings': {
            'handlers': ['console', 'file'],
            'level': 'DEBUG',
            'propagate': False,
        }
    }
}


def _configure_logging(loggingc, logfile=None):
    if loggingc is None:
        loggingc = deepcopy(LOGGING_DEFAULT_CONFIG)
        if logfile is None or logfile.lower() in ('none', 'null'):
            del loggingc['handlers']['file']
            loggingc['loggers']['anchorna']['handlers'] = ['console']
            loggingc['loggers']['py.warnings']['handlers'] = ['console']
        else:
            loggingc['handlers']['file']['filename'] = logfile
    logging.config.dictConfig(loggingc)
    logging.captureWarnings(loggingc.get('capture_warnings', False))


def _start_ipy(obj):
    from IPython import start_ipython
    print('Anchors loaded into anchors variable.')
    start_ipython(argv=[], user_ns={'anchors': obj}, display_banner=False)
    print('Bye')


def _tutorial_seqs(subset=False):
    seqs = read(files('anchorna.tests.data').joinpath('pesti56.gff.zip'))
    if subset:
        seqs = seqs[18:28]
        start, stop = zip(*[seq.fts.get('cds').loc_range for seq in seqs])
        start = max(start)
        stop = min(stop)
        stop = start + (stop-start) // 3 * 3 + 6
        for i in range(len(seqs)):
            seqs[i] = seqs[i][:start + 999] + seqs[i][stop:]
            seqs[i].fts[0].loc.stop -= stop-start - 999
    return seqs


def _cmd_create(conf, tutorial=False, tutorial_subset=False, no_cds=False):
    if tutorial and tutorial_subset:
        raise ValueError('Only one of tutorial or tutorial_subset are allowed')
    example = EXAMPLE_TOML_CONFIG
    if no_cds:
        example = example.replace('blosum62"', 'blosum62"       # use e.g. nuc for nucleotite sequences')
        example = example + 'no_cds = true              # directly use aa or nucleotite sequences, no translation, mode option is ignored\n'
    with open(conf, 'w') as f:
        f.write(example)
    if tutorial or tutorial_subset:
        seqs = _tutorial_seqs(subset=tutorial_subset)
        if no_cds:
            seqs = seqs['cds'].translate()
        path = os.path.dirname(conf)
        seqs.write(os.path.join(path, 'pesti_example.gff'))

def _cmd_go(fname, fname_anchor, pbar=True, continue_with=None,
            removed_anchors_path=None,
            logconf=None, logfile=None, no_logging=False,
            **kw):
    if no_logging:
        logging.disable(logging.CRITICAL)
    else:
        _configure_logging(logconf, logfile)
    log = logging.getLogger('anchorna')
    kw.pop('score_use_fluke', None)
    log.info(f'anchorna go is called with arguments {fname=}, {fname_anchor=}, '
             f'{continue_with=}, {removed_anchors_path=}')
    log.info(f'options={kw}')
    seqs = read(fname)
    if continue_with is not None:
        continue_with = read_anchors(continue_with)
    anchors, removed_anchors = find_my_anchors(
        seqs, continue_with=continue_with, **kw)
    # options = Options(**kw)
    rap = removed_anchors_path
    if removed_anchors is not None:
        if rap is None or rap.lower() in ('none', 'null'):
            log.debug('Do not save removed anchors')
        else:
            rap = rap.format(time.strftime('%Y_%m_%d-%H_%M_%S'))
            log.info(f'Save removed anchors in file {rap}')
            removed_anchors.write(rap)
    else:
        assert removed_anchors is None
    log.info(f'Write anchors to file {fname_anchor}')
    anchors.write(fname_anchor)

def _cmd_print(fname_anchor, verbose=False, mode='aa'):
    anchors = load_selected_anchors(fname_anchor)
    try:
        print(anchors.tostr(verbose=verbose, mode=mode))
    except BrokenPipeError:
        pass

def _cmd_load(fname_anchor):
    _start_ipy(read_anchors(fname_anchor))

def _cmd_export(fname_anchor, out, mode='aa', score_use_fluke=None, jalview=False):
    assert mode in ('seq', 'cds', 'aa')
    anchors = load_selected_anchors(fname_anchor)
    if jalview:
        outstr = jalview_features(anchors, mode=mode, score_use_fluke=score_use_fluke)
        if out is None:
            print(outstr)
        else:
            with open(out, 'w') as f:
                f.write(outstr)
    else:
        if out is None:
            out = sys.stdout
        anchors.write(out, mode=mode)

def _cmd_view(fname_anchor, fname, mode='aa', no_cds=False, align=None, score_use_fluke=None):
    assert mode in ('seq', 'cds', 'aa')
    with tempfile.TemporaryDirectory(prefix='anchorna') as tmpdir:
        fname_export = Path(tmpdir) / 'jalview_features.txt'
        fnameseq = Path(tmpdir) / 'aa_or_seq_or_cds.fasta'
        _cmd_export(fname_anchor, fname_export, mode=mode, score_use_fluke=score_use_fluke, jalview=True)
        seqs = read(fname)
        if mode != 'seq' and not no_cds:
            anchors = read_anchors(fname_anchor)
            offsets = {f.seqid: f.offset for anchor in anchors for f in anchor}
            for seq in seqs:
                seq.data = seq.data[offsets[seq.id]:]
            if mode == 'aa':
                seqs = seqs.translate(check_start=False)
        if align:
            anchors = read_anchors(fname_anchor)
            anchor = anchors[int(align.lower().removeprefix('a'))]
            start = max(_apply_mode(fluke.start, fluke.offset, mode=mode) for fluke in anchor)
            for seq in seqs:
                fluke = anchor.d[seq.id]
                seq.data = '-' * (start - _apply_mode(fluke.start, fluke.offset, mode=mode)) + seq.data
        seqs.write(fnameseq)
        subprocess.run(f'jalview {fnameseq} --features {fname_export}'.split())

def _cmd_combine(fname_anchor, out):
    lot_of_anchors = [load_selected_anchors(fn) for fn in fname_anchor]
    anchors = combine(lot_of_anchors)
    anchors.write(out)

def _cmd_cutout(fname, fname_anchor, pos1, pos2, out, fmt, mode='seq', score_use_fluke=None, no_cds=False):
    assert mode in ('seq', 'cds', 'aa')
    if no_cds:
        mode = 'aa'
    seqs = read(fname)
    anchors = load_selected_anchors(fname_anchor)
    seqs2 = cutout(seqs, anchors, pos1, pos2, mode=mode, score_use_fluke=score_use_fluke)
    if out is None:
        print(seqs2.tofmtstr(fmt or 'fasta'))
    else:
        seqs2.write(out, fmt)


def run(command, conf=None, pdb=False, **args):
    if command == 'create':
        return _cmd_create(conf, **args)
    if conf:
        try:
            with open(conf, 'rb') as f:
                conf = tomllib.load(f)
        except ValueError as ex:
            sys.exit('Error while parsing the configuration: %s' % ex)
        except IOError as ex:
            if command in ('go', 'export', 'view', 'cutout'):
                warn(ex)
            conf = {}
    else:
        conf = {}
    if command in ('export', 'view', 'cutout'):
        # ignore all config settings except the following
        conf = {k: conf[k] for k in ('fname', 'score_use_fluke', 'no_cds') if k in conf}
    # Populate args with conf, but prefer args
    conf.update(args)
    args = conf
    if pdb:
        import traceback
        import pdb
        def info(type, value, tb):
            traceback.print_exception(type, value, tb)
            print()
            pdb.pm()
        sys.excepthook = info
    # for no_cds, commands create, go, view and cutout expect this argument,
    # for other commands the default mode is already the correct one
    if getattr(args, 'pop' if command == 'export' else 'get')('no_cds', False):
        if 'mode' in args:
            raise ValueError('Mode may not be specified with --no-cds option')
    try:
        globals()[f'_cmd_{command}'](**args)
    except KeyError:
        raise ValueError(f'Unknown command: {command}')


def run_cmdline(cmd_args=None):
    """Main entry point from the command line"""
    from anchorna import __version__
    msg = 'AnchoRNA: Find anchors in short sequences of RNA/DNA'
    epilog = ('To get help on a subcommand run: anchorna command -h. '
              'Each time you provide a anchor file, you may load only '
              'selected anchors from the file using a dedicated syntax, '
              'fname_anchor|select or fname_anchor|select|remove, e.g. '
              'fname_anchor|4:10,12 only loads anchors 4:10, 12 (python indexing) - '
              'fname_anchor||a3 will remove the anchor a3 and fname_anchor|0:10|3 will '
              'select the first 10 anchors and remove a3. '
              'Numbers may be prepended with a single "a". '
              'The character - can be used to load a anchor file from stdin.')
    parser = argparse.ArgumentParser(description=msg, epilog=epilog)
    version = '%(prog)s ' + __version__
    parser.add_argument('--version', action='version', version=version)
    parser.add_argument('--pdb', action='store_true', help='Start the debugger upon exception')

    sub = parser.add_subparsers(title='commands', dest='command')
    sub.required = True

    p_create = sub.add_parser('create', help=(msg:='create example configuration or tutorial'), description=msg)
    p_go = sub.add_parser('go', help=(msg:='find anchors and write them into anchor file'), description=msg)
    p_print = sub.add_parser('print', help=(msg:='print contents of anchor file'), description=msg)
    p_load = sub.add_parser('load', help=(msg:='load anchors into IPython session'), description=msg)
    p_export = sub.add_parser('export', help=(msg:='export anchors to gff file (possibly using different mode than aa) or a JalView feature file'), description=msg)
    p_view = sub.add_parser('view', help=(msg:='view anchors in JalView (i.e. call export and start JalView)'), description=msg)
    msg = 'combine (and select or remove) anchors from one file or different files'
    msg2 = ('An expression consists of a file name fname_anchor or fname_anchor|selection to select anchors or fname_anchor||expression to remove anchors, or fname|selection|expression. '
            'Example call: anchorna combine f1.gff "f2.gff|a3:a5,a7" "f3.gff|a10" --out f5.gff, see also anchorna -h')
    p_combine = sub.add_parser('combine', help=msg, description=msg, epilog=msg2)
    msg = 'cut out sequences between two anchors and export them to new sequence file'
    msg2 = ('This command serves two purposes: '
            '1. Cutout sequences to put them again into the go command and '
            'find more anchors with updated options. '
            'Please use a file format which preserves the meta.offset attributes, e.g. sjson; use mode seq. '
            'Combine the found anchors with the anchorna combine command. '
            '2. Cutout sequences to use them in external tool. '
            'The cutout command expects two positions in the sequence. '
            'Each possition has 3 parts ABC where B and C are optional. '
            'Part A: Is a number or number prepended with letter a to specify the anchor number, '
            'use special words "start" and "end" for start or end of sequence, '
            'use special words "ATG" and "*" for start or stop codon of sequence (only allowed in mode "seq") '
            'Part B: One of the characters <, >, ^, for start, end or middle of word (anchor) specified in A, '
            'default is ^ (middle), must be ommitted for A=start or A=end. '
            'Part C: Additional character offset in the form +X or -X. '
            'Examples: anchorna cutout anchors.gff "a11<" "a12>+10", anchorna cutout anchors.gff a10 "*>"'
            )
    p_cutout = sub.add_parser('cutout', help=msg, description=msg, epilog=msg2)


    for p in (p_create, p_go, p_view, p_cutout):
        p.add_argument('-c', '--conf', default='anchorna.conf', help='configuration file to use (default: anchorna.conf)')
    p_create.add_argument('--tutorial', action='store_true', help='copy GFF sequence file for tutorial')
    p_create.add_argument('--no-cds', action='store_true', help='add the no_cds line to configuration, effect: directly use aa or nucleotide '
                                                                'sequences, no translation, --mode parameter not allowed')
    p_create.add_argument('--tutorial-subset', action='store_true', help=argparse.SUPPRESS, default=argparse.SUPPRESS)  # for testing purposes
    p_go.add_argument('fname_anchor', help='anchor file name (GFF output)')
    p_go.add_argument('--njobs', help='number of jobs to use', type=int, default=argparse.SUPPRESS)
    p_go.add_argument('--no-pbar', help='do not show progress bar', action='store_false', dest='pbar', default=argparse.SUPPRESS)
    p_go.add_argument('--no-remove', help='do not remove contradicting anchors', action='store_false', dest='remove', default=argparse.SUPPRESS)
    p_go.add_argument('--continue-with', help='only remove contradicting anchors from passed anchor list', default=argparse.SUPPRESS)
    p_go.add_argument('--no-logging', action='store_true', help=argparse.SUPPRESS, default=argparse.SUPPRESS)  # for testing purposes

    g = p_go.add_argument_group('optional arguments', description='Use these flags to overwrite values in the config file.')
    features = [(int, ('w', 'search-range')),
                (str, ('fname', 'gseqid', 'scoring', 'logfile', 'removed-anchors-path')),
                (float, ('score-add-word', 'thr-quota-add-anchor', 'thr-score-add-anchor'))]
    for type_, fs in features:
        for f in fs:
            g.add_argument('--' + f, default=argparse.SUPPRESS, type=type_)
    g.add_argument('--aggressive-remove', action=argparse.BooleanOptionalAction, default=argparse.SUPPRESS)
    g.add_argument('--no-cds', action='store_true', default=argparse.SUPPRESS)
    for p in (p_export, p_view, p_cutout):
        g = p.add_argument_group('optional arguments', description='Use these flags to overwrite values in the config file.')
        if p != p_export:
            g.add_argument('--fname', default=argparse.SUPPRESS)
        g.add_argument('--score-use-fluke', default=argparse.SUPPRESS, type=int)
        g.add_argument('--no-cds', action='store_true', default=argparse.SUPPRESS)


    p_export.add_argument('fname_anchor', help='anchor file name')
    p_export.add_argument('-o', '--out', help='output file name (by default prints to stdout)')
    p_export.add_argument('--jalview', action='store_true', help='export to JalView feature file instead of GFF file')
    p_view.add_argument('fname_anchor', help='anchor file name')
    p_view.add_argument('--align', help='align sequences at given anchor')
    p_combine.add_argument('fname_anchor', nargs='+', help='anchor file name')
    p_combine.add_argument('-o', '--out', help='output file name (by default prints to stdout)')

    p_cutout.add_argument('fname_anchor', help='anchor file name')
    p_cutout.add_argument('pos1', help='left position')
    p_cutout.add_argument('pos2', help='right position')

    p_cutout.add_argument('-o', '--out', help='output file name (by default prints fasta to stdout)')
    p_cutout.add_argument('--fmt', help='output format (default: autodetect from file extension)')

    p_print.add_argument('fname_anchor', help='anchor file name')
    p_print.add_argument('-v', '--verbose', help=msg, action='store_true')
    p_load.add_argument('fname_anchor', help='anchor file name')

    for p in (p_print, p_export, p_view, p_cutout):
        choices = ['seq', 'cds', 'aa']
        default = 'seq' if p == p_cutout else 'aa'
        msg = ('choose mode, seq: relative to original sequence, '
               'cds: accounting for offset (usually coding sequence), aa: translated cds '
               f'(default: {default})')
        p.add_argument('-m', '--mode', default=argparse.SUPPRESS, choices=choices, help=msg)

    # Get command line arguments and start run function
    args = parser.parse_args(cmd_args)
    run(**vars(args))

if __name__ == '__main__':
    run_cmdline()
