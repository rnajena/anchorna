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
    'MW054939 AF526381 JQ799141 NC_076032 NC_039237 MH231152 KT875135 '
    'MH806436 HG426490 MK599227 JN380086 NC_038912 NC_002657 OM817567 '
    'GU233732 KP343640 KC533775 AY805221 AF099102 AY775178 MG387218 '
    'NC_003679 MF102261 NC_024018 NC_023176 NC_003678 NC_012812 NC_018713 '
    'NC_025677 NC_038964 NC_030653 MK216749 MH221025 MN099165 MN584738 '
    'NC_035432 NC_077026 MZ664274 NC_077023 NC_077024 ON024093 ON024108 '
    'NC_077001 NC_077000 NC_077015 OU592965 OM030319 OM030320').split()


def create_example_seqs_file():
    """
    Function to create the example sequences file in the tests/data folder
    """
    from sugar.web import Entrez
    seqs = Entrez().get_basket(_IDS)
    for seq in seqs:
        seq.fts = seq.fts.select('cds')
    seqs.write('data/pesti56.gff', archive='zip')


@contextlib.contextmanager
def _changetmpdir():
    origin = Path().resolve()
    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            os.chdir(tmpdir)
            yield Path(tmpdir)
        finally:
            os.chdir(origin)


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


def test_reproduce_anchor_file_subset():
    with _changetmpdir() as tmpdir:
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


def test_anchorna_workflow_subset():
    """
    Tests the full anchorna workflow with a subset of the example sequences
    """
    with _changetmpdir() as tmpdir:
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


def test_anchorna_workflow_subset_poor():
    with _changetmpdir() as tmpdir:
        assert '' == check('anchorna create')
        fname_seqs = tmpdir / 'pesti_example.gff'
        assert '' == check('anchorna create --tutorial-subset')
        assert '' == check('anchorna go --thr-quota-add-anchor 0.5 --score-add-word 18 --no-pbar anchors.gff')
        assert 'A11' in check('anchorna print anchors.gff')
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


def test_anchorna_workflow_subset_no_cds():
    with _changetmpdir() as tmpdir:
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
def test_reproduce_anchor_file_complete():
    with _changetmpdir() as tmpdir:
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
def test_tutorial():
    fname = files('anchorna').joinpath('../README.md')
    if not os.path.exists(fname):
        pytest.skip('README.md only available in dev install')
    with open(fname) as f:
        readme = f.read()
    i1 = readme.find('anchorna go')
    i2 = readme.find('```', i1) + i1
    tutorial = 'anchorna create --tutorial\n' + readme[i1:i2]
    tutorial = tutorial.replace('"A??>" "A??<"', '"A33>" "A34<"')
    sys.modules['IPython'] = MagicMock()
    with _changetmpdir():
        with patch('subprocess.run'):
            for line in tutorial.splitlines():
                if line.startswith('anchorna'):
                    line = line.replace('| anchorna view -', '').replace('"', '')
                    line = line.split('#')[0].strip()
                    check(line)
