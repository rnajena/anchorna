# (C) 2024, Tom Eulenfeld, MIT license

import argparse
import sys
import contextlib
from pathlib import Path
import os
import logging
import json


@contextlib.contextmanager
def _changedir(path):
    origin = Path().absolute()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(origin)


def _start_ipy(obj):
    from IPython import start_ipython
    print('Anchors loaded into anchors variable.')
    start_ipython(argv=[], user_ns={'anchors': obj}, display_banner=False)
    print('Bye')


def parse_selection(anchors, selection):
    anchors2 = anchors[:0]
    for part in selection.split(','):
        if ':' in part:
            i, j = part.split(':')
            i = int(i.strip().removeprefix('a'))
            j = int(j.strip().removeprefix('a'))
        else:
            i = int(part.strip().removeprefix('a'))
            j = i + 1
        anchors2.extend(anchors[i:j])
    return anchors2

def load_anchors(fname):
    from anchorna import load_json
    if '|' not in fname:
        return load_json(fname)
    fname, selection = fname.split('|', 1)
    anchors = load_json(fname.strip())
    selection = selection.lower()
    if '|' not in selection:
        return parse_selection(anchors, selection)
    selection, rem = selection.split('|')
    anchors_remove = set(parse_selection(anchors, rem))
    anchors = parse_selection(anchors, selection) if selection.strip() else anchors
    anchors.data = [a for a in anchors if a not in anchors_remove]
    return anchors

def _cutout(fname_seq, fname_anchor,
            left_anchor, right_anchor, left, right,
            mode='default',
            # offset='nt',
            out=None,
            fmt=None):
    from sugar import read, BioBasket
    from anchorna.util import _apply_mode

    seqs = read(fname_seq).d
    anchors = load_anchors(fname_anchor)
    l = left_anchor.strip().lower()
    r = right_anchor.strip().lower()
    if l == 'start':
        i = 0
    else:
        l = anchors[int(l.removeprefix('a'))].d
        ids = l.keys()
    if r == 'end':
        j = None
    else:
        r = anchors[int(r.removeprefix('a'))].d
        ids = r.keys()
    if l == 'start' and r == 'end':
        raise ValueError("You don't want to cut out the whole sequence, do you?")
    seqs2 = BioBasket()
    for id_ in ids:
        if id_ not in seqs:
            import warnings
            warnings.warn(f'Id {id_} not present in sequences')
            continue
        if l != 'start':
            i1 = _apply_mode(l[id_].start, l[id_].offset, mode=mode)
            i2 = _apply_mode(l[id_].stop, l[id_].offset, mode=mode)
            i = i1 if left == 'full' else i2 if left == 'no' else (i1+i2)//2
        if r != 'end':
            j1 = _apply_mode(r[id_].start, r[id_].offset, mode=mode)
            j2 = _apply_mode(r[id_].stop, r[id_].offset, mode=mode)
            j = j2 if right == 'full' else j1 if right == 'no' else (j1+j2)//2
        seq = seqs[id_][i:j]
        offset = i if seq.type == 'nt' else 3 * i
        seq.meta.offset = seq.meta.get('offset', 0) + offset
        seqs2.append(seq)
    if out is None:
        print(seqs2.tofmtstr(fmt or 'fasta'))
    else:
        seqs2.write(out, fmt)


def run(command, pytest_args, mode='default', out=None, conf=None, pdb=False, **args):
    from sugar import read
    from anchorna import get_anchors, jalview_features, load_json, write_json
    from anchorna.util import AnchorList, ConfigJSONDecoder, Options

    if pdb:
        import traceback
        import pdb
        def info(type, value, tb):
            traceback.print_exception(type, value, tb)
            print()
            pdb.pm()
        sys.excepthook = info
    if command == 'print':
        anchors = load_anchors(args['fname'])
        try:
            print(anchors.tostr(verbose=args['verbose'], mode=mode))
        except BrokenPipeError:
            pass
    elif command == 'load':
        _start_ipy(load_json(args['fname']))
    elif command == 'go':
        logging.basicConfig(level=logging.DEBUG)
        if conf:
            try:
                with open(conf) as f:
                    conf = json.load(f, cls=ConfigJSONDecoder)
            except ValueError as ex:
                print('Error while parsing the configuration: %s' % ex)
                return
            except IOError as ex:
                print(ex)
                return
            # Populate args with conf, but prefer args
            conf.update(args)
        args = conf
        seqs = read(args.pop('fname_seqs'))
        out = args.pop('fname_anchors')
        options = Options(**args)
        anchors = get_anchors(seqs, options)
        write_json(anchors, out)
    elif command == 'combine':
        anchors = set()
        for fn in args['fname']:
            nans = set(load_anchors(fn))
            if len(nans & anchors) > 0:
                ids = ', '.join(a.id for a in nans & anchors)
                raise ValueError(f'Anchors {ids} exits in multiple files')
            anchors |= nans
        anchors = AnchorList(sorted(anchors, key=lambda a: a.ref.start))
        write_json(anchors, out)
    elif command == 'export':
        anchors = load_anchors(args['fname'])
        outstr = jalview_features(anchors, mode=mode)
        if out is None:
            print(outstr)
        else:
            with open(out, 'w') as f:
                f.write(outstr)
    elif command == 'cutout':
        _cutout(out=out, **args)
    elif command == 'test':
        try:
            import pytest
        except ImportError:
            msg = ("\nanchorna's test suite uses pytest. "
                   "Please install pytest before running the tests.")
            sys.exit(msg)
        path = Path(__file__).parent
        print(f'Run pytest in directory {path}')
        with _changedir(path):
            status = pytest.main(pytest_args)
        sys.exit(status)


def run_cmdline(cmd_args=None):
    """Main entry point from the command line"""
    from anchorna import __version__
    msg = 'anchoRNA: Find anchors in short sequences of RNA/DNA'
    epilog = ('To get help on a subcommand run: anchorna command -h. '
              'Each time you provide a anchor file, you may load only '
              'selected anchors from the file using a dedicated syntax, '
              'fname|select or fname|select|remove, e.g. '
              'fname|4:10,12 only loads anchors 4:10, 12 (python indexing) - '
              'fname||a3 will remove the anchor a3 and fname|0:10|3 will '
              'select the first 10 anchors and remove a3. '
              'Numbers may be prepended with a single "a". '
              'The character - can be used to load a anchor file from stdin.')
    parser = argparse.ArgumentParser(description=msg, epilog=epilog)
    version = '%(prog)s ' + __version__
    parser.add_argument('--version', action='version', version=version)

    sub = parser.add_subparsers(title='commands', dest='command')
    sub.required = True
    msg = 'find anchors and write them into anchor file'
    p_go = sub.add_parser('go', help=msg, description=msg)
    p_go.add_argument('fname_seqs', help='sequence file name')
    p_go.add_argument('fname_anchors', help='anchor file name (JSON output)')
    msg = 'Configuration file to use (default: conf.json)'
    p_go.add_argument('-c', '--conf', default='conf.json', help=msg)
    msg = 'Use these flags to overwrite values in the config file.'
    g = p_go.add_argument_group('optional arguments', description=msg)
    features = [(int, ('w', 'maxshift')),
                (str, ('refid', 'scoring')),
                (float, ('thr-score', 'thr-quota-score', 'thr-quota'))]
    for type_, fs in features:
        for f in fs:
            g.add_argument('--' + f, default=argparse.SUPPRESS, type=type_)

    msg = 'export anchors into feature file to view them in Jalview'
    p_export = sub.add_parser('export', help=msg, description=msg)
    p_export.add_argument('fname', help='anchor file name')
    p_export.add_argument('-o', '--out', help='output file name (by default prints to stdout)')


    msg = 'combine (and select or remove) anchors from one file or different files'
    msg2 = ('An expression consists of a file name fname or fname/s/selection to select acnhors or fname/r/expression to remove anchors. '
            'Example call: anchorna combine f1.json "f2.sjon|a3:a5,a7" "f3.json|a10|" --out f5.json')
    p_combine = sub.add_parser('combine', help=msg, description=msg, epilog=msg2)
    p_combine.add_argument('fname', nargs='+', help='anchor file name')
    p_combine.add_argument('-o', '--out', help='output file name (by default prints to stdout)')

    msg = 'cut out sequences between anchors and export them to new sequences file'
    msg2 = ('This command serves two purposes: '
            '1. Cutout sequences to put them again into the go command and find more anchors. '
            'In this case it might be important to use a file format that writes the offset metadata, '
            'e.g. sjson. '
            '2. Cutout sequences to use them in external tool.')
    p_cutout = sub.add_parser('cutout', help=msg, description=msg, epilog=msg2)
    p_cutout.add_argument('fname_seq', help='sequence file name')
    p_cutout.add_argument('fname_anchor', help='anchor file name')
    p_cutout.add_argument('left_anchor', help='left anchor, use start for start of sequence')
    p_cutout.add_argument('right_anchor', help='right anchor, use end for end of sequence')
    choices = ['full', 'no', 'mid']
    p_cutout.add_argument('-l', '--left', help='include left anchor, default: full', default='full', choices=choices)
    p_cutout.add_argument('-r', '--right', help='include right anchor, default: full', default='full', choices=choices)
    p_cutout.add_argument('-o', '--out', help='output file name (by default prints fasta to stdout)')
    p_cutout.add_argument('-f', '--fmt', help='format of file (default: autodetect)')
    # p_cutout.add_argument('--offset', help='', choices=['nt', 'aa'], default='nt')


    for p in (p_go, p_combine, p_cutout):
        msg = 'Start the debugger upon exception'
        p.add_argument('--pdb', action='store_true', help=msg)

    msg = 'print contents of anchor file'
    p_print = sub.add_parser('print', help=msg)
    p_print.add_argument('fname', help='anchor file name')
    p_print.add_argument('-v', '--verbose', help=msg, action='store_true')

    msg = 'load anchors into IPython session'
    p_load = sub.add_parser('load', help=msg)
    p_load.add_argument('fname', help='anchor file name')

    for p in (p_print, p_export, p_cutout):
        choices = ['default', 'offset', 'seq', 'seqoffset']
        msg = ('choose mode, offset: add offset, seq: multiply by 3, '
               'seqoffset: seq+offset ;)')
        p.add_argument('-m', '--mode', default='default', choices=choices, help=msg)

    # msg = 'run anchorna test suite'
    # msg2 = ('The test suite uses pytest. You can call pytest directly or use '
    #         'most of pytest cli arguments in the anchorna test call. '
    #         'See pytest -h.')
    # p_test = sub.add_parser('test', help=msg, description=msg2)


    # Get command line arguments and start run function
    args, pytest_args = parser.parse_known_args(cmd_args)
    if args.command != 'test':
        parser.parse_args(cmd_args)  # just call again to properly raise the error
    run(pytest_args=pytest_args, **vars(args))

if __name__ == '__main__':
    run_cmdline()
