# (C) 2024, Tom Eulenfeld, MIT license
import contextlib
from importlib.resources import files
import io
import os
from pathlib import Path
from subprocess import check_output
import sys
import tempfile
from unittest.mock import MagicMock, patch
from warnings import warn

import pytest
from sugar import read

from anchorna import cutout, read_anchors
from anchorna.cli import run_cmdline
from anchorna.io import load_json, write_json


_IDS = (  # Representative sequences of pesti virus
    'NC_076029 NC_001461 KC853440 JX419398 M96687 KT951841 ON165517 KX577637 '
    'MW054939 AF526381 NC_076032 NC_039237 MH231152 KT875135 '
    'MH806436 HG426490 MK599227 JN380086 NC_003678 NC_012812 NC_038912 NC_002657 OM817567 '
    'GU233732 KP343640 KC533775 AY805221 AF099102 AY775178 MG387218 '
    'NC_003679 MF102261 NC_018713 MZ664274 NC_077023 NC_024018 NC_023176 '
    'NC_035432 NC_077026 OM030320 NC_077024 ON024093 ON024108 '
    'NC_038964 NC_030653 MK216749 MH221025 MN099165 MN584738 '
    'NC_077015 OU592965 NC_077001 NC_025677 NC_077000 OM030319').split()


def create_example_seqs_file():
    """
    Function to create the example sequences file in the tests/data folder
    """
    from sugar.web import Entrez
    seqs = Entrez().get_basket(_IDS)
    for seq in seqs:
        seq.fts = seq.fts.select('cds')
    fname = files('anchorna.tests.data').joinpath('pesti55.gff')
    seqs.write(fname, archive='zip')


# @contextlib.contextmanager
# def __pseudo_tempdir(path):
#     yield path


# @contextlib.contextmanager
# def _changetmpdir(path=None):
#     origin = Path().resolve()
#     # ignore_cleanup_errors is necessary for windows
#     tmpdirmanager = tempfile.TemporaryDirectory(ignore_cleanup_errors=os.name == 'nt') if path is None else __pseudo_tempdir(path)
#     with tmpdirmanager as tmpdir:
#         try:
#             os.chdir(tmpdir)
#             yield Path(tmpdir)
#         finally:
#             os.chdir(origin)


@pytest.fixture()
def tmp_path_cd(tmp_path):
    """
    A fixture which changes the dir to a temporary directory
    """
    with contextlib.chdir(tmp_path):
        yield tmp_path


def check(cmd):
    args = cmd.split()
    assert args[0] == 'anchorna'
    with contextlib.redirect_stdout(io.StringIO()) as f:
        run_cmdline(args[1:])
    return f.getvalue()


def _fix_open_log_file_on_windows():
    # fix the following error observed in CI
    # FAILED ..\anchorna\tests\test_anchorna.py::test_anchorna_workflow_subset_poor
    # - PermissionError: [WinError 32] The process cannot access the file because
    # it is being used by another process: '...\\anchorna.log'
    import logging
    log = logging.getLogger('anchorna')
    for handler in log.handlers:
        handler.close()


def test_anchorna_script_help():
    """
    Test that anchorna can be called on the command line
    """
    assert b'Find anchors' in check_output('anchorna -h'.split())
    assert b'go' in check_output('anchorna go -h'.split())


def test_reproduce_anchor_file_subset(tmp_path_cd):
    tmpdir = tmp_path_cd
    check('anchorna create --tutorial-subset')
    check('anchorna go --no-pbar anchors.gff --no-logging --no-aggressive-remove')
    anchors = read_anchors('anchors.gff')
    fname = files('anchorna.tests.data').joinpath('anchors_subset.gff')
    try:
        anchors2 = read_anchors(fname)
    except FileNotFoundError:
        warn(f'Did not find test file {fname}, create it')
        import shutil
        shutil.copy(tmpdir / 'anchors.gff', fname)
        anchors2 = read_anchors(fname)
    assert anchors2 == anchors


def test_anchorna_workflow_subset(tmp_path_cd):
    """
    Tests the full anchorna workflow with a subset of the example sequences
    """
    tmpdir = tmp_path_cd
    assert '' == check('anchorna create')
    fname_seqs = tmpdir / 'pesti_example.gff'
    assert '' == check('anchorna create --tutorial-subset')
    assert '' == check('anchorna go --no-pbar anchors.gff')
    assert 'A11' in check('anchorna print anchors.gff')
    assert 'F1' in check('anchorna print anchors.gff -v')
    out1 = check('anchorna combine anchors.gff|a5:a10|a8')
    out2 = check('anchorna combine anchors.gff|a5,a6,a7,a9')
    assert out1 == out2
    assert 'anchor0' in check('anchorna export --fmt jalview anchors.gff')
    assert 'anchor0' in check('anchorna export --fmt jalview -m nt anchors.gff')
    assert 'anchor0' in check('anchorna export --fmt jalview -m aa anchors.gff')
    assert '' in check('anchorna export --fmt jalview -m cds anchors.gff -o test_jalview.txt')
    assert 'anchorna' in check('anchorna export anchors.gff')
    assert 'anchorna' in check('anchorna export -m nt anchors.gff')
    assert 'anchorna' in check('anchorna export -m aa anchors.gff')
    assert 'anchorna' in check('anchorna combine anchors.gff --convert-nt')
    assert '' in check('anchorna export -m cds anchors.gff -o test_anchor_export.gff')
    with pytest.raises(IOError):
        read_anchors('test_anchor_export.gff')
    read_anchors('test_anchor_export.gff', check_header=False)
    assert '1' in check('anchorna export --fmt dialign anchors.gff')
    assert '1' in check('anchorna export --fmt dialign -m nt anchors.gff')
    assert '1' in check('anchorna export --fmt dialign -m cds anchors.gff')
    with patch('subprocess.run'):  # we do not want to actually start jalview here
        assert '' == check('anchorna view anchors.gff')
        assert '' == check('anchorna view anchors.gff --align a1')

    anchors = read_anchors('anchors.gff')

    # test cutout
    assert '>' in check(f'anchorna cutout anchors.gff atg>+5 end-10 --fname {fname_seqs}')

    seqs = read(fname_seqs)
    seqs2 = cutout(seqs, anchors, 'start+10', 'a5^-5')
    seqs3 = cutout(seqs, anchors, 'a5^-5', '*>')
    seqs4 = cutout(seqs, anchors, '*>', 'end')
    assert str(seqs[0, 10:]) == seqs2[0].data + seqs3[0].data + seqs4[0].data

    fname = tmpdir / 'pesti_test_cutout.sjson'
    assert '' == check(f'anchorna cutout anchors.gff a0> a2< -o {fname}')
    assert '' == check(f'anchorna go --fname {fname} --no-pbar anchors_cutout.gff')
    assert '' == check('anchorna combine anchors.gff||a1 anchors_cutout.gff -o anchors_combined.gff')
    assert read_anchors('anchors_combined.gff') == anchors

    fname = tmpdir / 'pesti_test_cutout2.sjson'
    assert '' == check(f'anchorna cutout anchors.gff a6> a10< -o {fname}')
    assert '' == check(f'anchorna go --fname {fname} --no-pbar anchors_cutout2.gff --search-range=1000')
    assert '' == check('anchorna combine anchors.gff||a7:a10 anchors_cutout2.gff -o anchors_combined2.gff')
    assert read_anchors('anchors_combined2.gff') == read_anchors('anchors.gff')

    # check --no-remove option and --continue-with option
    assert '' == check('anchorna go --no-remove --no-pbar anchors2.gff')
    assert '' == check('anchorna go --continue-with anchors2.gff --no-pbar anchors3.gff')
    assert len(anchors) < len(read_anchors('anchors2.gff'))
    assert anchors == read_anchors('anchors3.gff')

    # test json
    json = tmpdir / 'anchors.json'
    write_json(anchors, json)
    assert load_json(json) == anchors

    _fix_open_log_file_on_windows()


def test_anchorna_workflow_subset_poor(tmp_path_cd):
    tmpdir = tmp_path_cd
    assert '' == check('anchorna create')
    fname_seqs = tmpdir / 'pesti_example.gff'
    assert '' == check('anchorna create --tutorial-subset')
    assert '' == check('anchorna go --thr-quota-add-anchor 0.5 --score-add-word 18 --no-pbar anchors.gff')
    assert 'A10' in check('anchorna print anchors.gff')
    assert '(poor)' in check('anchorna print anchors.gff -v')
    out1 = check('anchorna combine anchors.gff|a5:a10|a8')
    out2 = check('anchorna combine anchors.gff|a5,a6,a7,a9')
    assert out1 == out2
    assert 'anchor0' in check('anchorna export --fmt jalview anchors.gff')
    assert 'anchorna' in check('anchorna export anchors.gff')
    with patch('subprocess.run'):  # we do not want to actually start jalview here
        assert '' == check('anchorna view anchors.gff')
    anchors = read_anchors('anchors.gff')
    check(f'anchorna cutout anchors.gff atg>+5 end-10 --fname {fname_seqs}')
    seqs = read(fname_seqs)
    seqs2 = cutout(seqs, anchors, 'start+10', 'a5^-5')
    seqs3 = cutout(seqs, anchors, 'a5^-5', '*>')
    seqs4 = cutout(seqs, anchors, '*>', 'end')
    assert str(seqs[0, 10:]) == seqs2[0].data + seqs3[0].data + seqs4[0].data

    _fix_open_log_file_on_windows()


def test_anchorna_workflow_subset_no_cds(tmp_path_cd):
    tmpdir = tmp_path_cd
    assert '' == check('anchorna create')
    fname_seqs = tmpdir / 'pesti_example.gff'
    assert '' == check('anchorna create --tutorial-subset --no-cds')
    assert '' == check('anchorna go --no-pbar anchors.gff')
    anchors1 = read_anchors('anchors.gff')
    assert anchors1.no_cds

    assert '' == check('anchorna create --tutorial-subset')
    assert '' == check('anchorna go --no-pbar anchors_cds.gff')
    anchors2 = read_anchors('anchors_cds.gff')
    anchors2.no_cds = True
    for anchor in anchors2:
        for fluke in anchor:
            fluke.offset = 0
    assert anchors1 == anchors2

    assert 'A11' in check('anchorna print anchors.gff')
    assert 'F1' in check('anchorna print anchors.gff -v')
    out1 = check('anchorna combine anchors.gff|a5:a10|a8')
    out2 = check('anchorna combine anchors.gff|a5,a6,a7,a9')
    assert out1 == out2
    assert 'anchor0' in check('anchorna export --fmt jalview anchors.gff')
    assert 'anchor0' in check('anchorna export --fmt jalview -m nt anchors.gff')
    assert 'anchorna' in check('anchorna export anchors.gff')
    with patch('subprocess.run'):  # we do not want to actually start jalview here
        assert '' == check('anchorna view anchors.gff')
        assert '' == check('anchorna view anchors.gff --align a1')
    anchors = read_anchors('anchors.gff')
    with pytest.raises(ValueError):
        check(f'anchorna cutout anchors.gff atg>+5 end-10 --fname {fname_seqs}')
    seqs = read(fname_seqs)
    seqs2 = cutout(seqs, anchors, 'start+10', 'a5^-5')
    seqs3 = cutout(seqs, anchors, 'a5^-5', '*>')
    seqs4 = cutout(seqs, anchors, '*>', 'end')
    assert str(seqs[0, 10:]) == seqs2[0].data + seqs3[0].data + seqs4[0].data

    _fix_open_log_file_on_windows()


@pytest.mark.slowtest
def test_reproduce_anchor_file_complete(tmp_path_cd):
    tmpdir = tmp_path_cd
    check('anchorna create --tutorial')
    check('anchorna go --no-pbar anchors.gff --no-logging --njobs=-1')
    anchors = read_anchors('anchors.gff')
    fname = files('anchorna.tests.data').joinpath('anchors_complete.gff')
    try:
        anchors2 = read_anchors(fname)
    except FileNotFoundError:
        warn(f'Did not find test file {fname}, create it')
        import shutil
        shutil.copy(tmpdir / 'anchors.gff', fname)
        anchors2 = read_anchors(fname)
    assert anchors2 == anchors


@pytest.mark.slowtest
def test_tutorial(tmp_path_cd):
    fname = files('anchorna').joinpath('../README.md')
    if not os.path.exists(fname):
        pytest.skip('README.md only available in dev install')
    with open(fname) as f:
        readme = f.read()
    i1 = readme.find('anchorna go')
    i2 = readme.find('```', i1) + i1
    tutorial = 'anchorna create --tutorial\n' + readme[i1:i2]
    tutorial = tutorial.replace('"A??>" "A??<"', '"A40>" "A41<"')
    sys.modules['IPython'] = MagicMock()
    with patch('subprocess.run'):
        for line in tutorial.splitlines():
            if line.startswith('anchorna'):
                line = line.replace('| anchorna view -', '').replace('"', '')
                line = line.split('#')[0].strip()
                check(line)


def test_export_stockholm(tmp_path_cd):
    stockholmf = ('# STOCKHOLM 1.0\n'
                  'S1 ---gGTATACG--\n'
                  'S2 -ggggtatacc--\n'
                 #'   ....|-A0-|...\n'
                  )
    gcrow = 'GC AnchoRNA ....|-A0-|...'
    anchorf = (
        '##gff-version 3\n'
        '#AnchoRNA anchor file\n'
        '#offset S1 1\n'
        '#offset S2 0\n'
        'S1	anchorna	anchor	1	2	28	+	.	word=BLA;median_score=22;Name=A0\n'
        'S2	anchorna	fluke	2	3	28	+	.	word=UPS;median_score=22;Name=A0_S2\n'
    )
    with open('seqs.stk', 'w') as f:
        f.write(stockholmf)
    with open('anchors.gff', 'w') as f:
        f.write(anchorf)
    check('anchorna create')
    check('anchorna export anchors.gff --fname seqs.stk --fmt stockholm --out seqs_gc.stk')
    with open('seqs_gc.stk') as f:
        content = f.read()
    assert gcrow in content


def test_anchorna_multiple_cds(tmp_path_cd):
    assert '' == check('anchorna create')
    assert '' == check('anchorna create --tutorial-subset')
    assert '' == check('anchorna go --no-pbar anchors.gff')
    out1 = check('anchorna print anchors.gff')
    out1v = check('anchorna print anchors.gff -v --mode nt')
    # we create two CDS which have the same anchors
    anchors = read_anchors('anchors.gff')
    fts1 = anchors[0:1].convert2fts(mode='nt')
    fts2 = anchors[3:4].convert2fts(mode='nt')
    fts3 = anchors[4:5].convert2fts(mode='nt')
    fts4 = anchors[-1:].convert2fts(mode='nt')
    for ft1, ft2 in zip(fts1, fts2):
        assert ft1.seqid == ft2.seqid
        ft1.loc.stop = ft2.loc.stop
        ft1.type = 'CDS'
    for ft1, ft2 in zip(fts3, fts4):
        assert ft1.seqid == ft2.seqid
        ft1.loc.stop = ft2.loc.stop
        ft1.type = 'CDS'
    seqs = read('pesti_example.gff')
    seqs.fts = fts1
    seqs.write('cds1.gff')
    seqs.fts = fts3
    seqs.write('cds2.gff')
    # run anchorna individually on two CDS and combine results
    with pytest.warns():  # first codon is not a start codon
        assert '' == check('anchorna go --fname cds1.gff --no-pbar anchors_cds1.gff')
        assert '' == check('anchorna go --fname cds2.gff --no-pbar anchors_cds2.gff')
    assert '' == check('anchorna combine anchors_cds1.gff anchors_cds2.gff -o anchors_both_cds.gff')
    assert '' == check('anchorna combine anchors_cds1.gff anchors_cds2.gff -o anchors_both_cds_nt.gff --convert-nt')
    out2 = check('anchorna print anchors_both_cds.gff')
    out2v = check('anchorna print anchors_both_cds.gff -v --mode nt')
    out3v = check('anchorna print anchors_both_cds_nt.gff -v')
    # for print with mode aa, we have the wrong offset
    for line1, line2 in zip(out1.strip().splitlines(), out2.strip().splitlines()):
        assert line1.split('+')[1] == line2.split('+')[1]
    # perfect for mode nt
    assert out2v == out1v
    assert out3v == out1v


def test_anchorna_antisense(tmp_path_cd):
    assert '' == check('anchorna create')
    assert '' == check('anchorna create --tutorial-subset')
    assert '' == check('anchorna go --no-pbar anchors.gff')
    out1v = check('anchorna print anchors.gff -v')
    seqs1 = read('pesti_example.gff')
    seqs2 = seqs1.copy().rc(update_fts=True)
    seqs2.write('pestirc.gff')
    assert '' == check('anchorna go --no-pbar anchorsrc.gff --fname pestirc.gff')
    out2v = check('anchorna print anchorsrc.gff -v')
    assert out2v == out1v
    anchors1 = read_anchors('anchors.gff')
    anchors2 = read_anchors('anchorsrc.gff')
    c1 = cutout(seqs1, anchors1, '3', '10')
    c2 = cutout(seqs2, anchors2, '10', '3')
    for s in c1 + c2:
        del s.meta.offset
    assert c2.rc() == c1
    c1 = cutout(seqs1.copy()['cds'].translate(), anchors1, '3', '10^', mode='aa')
    c2 = cutout(seqs2.copy()['cds'].translate(), anchors2, '3', '10^', mode='aa')
    for s in c1 + c2:
        del s.meta.fts
    assert c2 == c1
    c1 = cutout(seqs1.copy()['cds'], anchors1, '3', '10^', mode='cds')
    c2 = cutout(seqs2.copy()['cds'], anchors2, '3', '10^', mode='cds')
    for s in c1 + c2:
        del s.meta.fts
    assert c2 == c1


def test_anchorna_multiple_cds_antisense(tmp_path_cd):
    assert '' == check('anchorna create')
    assert '' == check('anchorna create --tutorial-subset')
    assert '' == check('anchorna go --no-pbar anchors.gff')
    # we create two CDS, one on + one on - strand, both have the same anchors
    anchors = read_anchors('anchors.gff')
    fts1 = anchors[0:1].convert2fts(mode='nt')
    fts2 = anchors[3:4].convert2fts(mode='nt')
    for ft1, ft2 in zip(fts1, fts2):
        assert ft1.seqid == ft2.seqid
        ft1.loc.stop = ft2.loc.stop
        ft1.type = 'CDS'
    seqs = read('pesti_example.gff')
    seqs.fts = fts1
    seqs = (seqs + seqs.copy().rc(update_fts=True)).merge(update_fts=True)
    seqs.write('combined.gff')
    seqs1 = seqs.copy()
    seqs1.fts = seqs1.fts.select(strand='+')
    seqs1.write('cds1.gff')
    seqs2 = seqs.copy()
    seqs2.fts = seqs2.fts.select(strand='-')
    seqs2.write('cds2.gff')
    # run anchorna individually on two CDS and combine results
    with pytest.warns():  # first codon is not a start codon
        assert '' == check('anchorna go --fname combined.gff --no-pbar anchors_combined.gff')
        assert '' == check('anchorna go --fname cds1.gff --no-pbar anchors_cds1.gff')
        assert '' == check('anchorna go --fname cds2.gff --no-pbar anchors_cds2.gff')
    with pytest.raises(ValueError, match='Cannot combine'):
        check('anchorna combine anchors_cds1.gff anchors_cds2.gff')
    assert '' == check('anchorna combine anchors_cds1.gff anchors_cds2.gff -o anchors_both_cds_nt.gff --convert-nt')
    out2 = check('anchorna print anchors_combined.gff --mode nt')
    out2v = check('anchorna print anchors_combined.gff -v --mode nt')
    out3 = check('anchorna print anchors_both_cds_nt.gff')
    out3v = check('anchorna print anchors_both_cds_nt.gff -v')
    assert out2 in out3
    assert out2v in out3v

    anchors = read_anchors('anchors_both_cds_nt.gff')
    c1 = cutout(seqs, anchors, '0', '2', mode='aa')
    c2 = cutout(seqs, anchors, '5', '7', mode='aa')
    assert c2.rc() == c1
    assert '' == check('anchorna cutout anchors_both_cds_nt.gff 5 7 --fname combined.gff -o cutout.fasta')
    c3 = read('cutout.fasta')
    for s in c1 + c3:
        s.meta = {}
    assert c3.rc() == c1
