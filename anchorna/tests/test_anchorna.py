# (C) 2024, Tom Eulenfeld, MIT license
import contextlib
import io
import os
from pathlib import Path
from subprocess import check_output
import tempfile
from unittest.mock import patch

from sugar import read

from anchorna import cutout, read_anchors
from anchorna.cli import run_cmdline
from anchorna.io import load_json, load_selected_anchors, write_json


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
def _changedir(path):
    origin = Path().resolve()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(origin)

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
        assert 'anchor0' in check('anchorna export anchors.gff')
        assert 'anchor0' in check('anchorna export -m seq anchors.gff')
        assert 'anchor0' in check('anchorna export -m aa anchors.gff')
        assert 'anchor0' in check('anchorna export -m cds anchors.gff')
        with patch('subprocess.run'):  # we do not want to actually start jalview here
            assert '' == check('anchorna view anchors.gff')
            assert '' == check('anchorna view anchors.gff --align a1')

        anchors = read_anchors('anchors.gff')
        assert len(anchors) == 15

        # test cutout
        assert '>' in check(f'anchorna cutout anchors.gff atg>+5 end-10 --fname {fname_seqs}')

        seqs = read(fname_seqs)
        seqs2 = cutout(seqs, anchors, 'start+10', 'a5-5')
        seqs3 = cutout(seqs, anchors, 'a5-5', '*>')
        seqs4 = cutout(seqs, anchors, '*>', 'end')
        assert str(seqs[0, 10:]) == str(seqs2[0] + seqs3[0] + seqs4[0])

        fname = tmpdir / 'pesti_test_cutout.sjson'
        assert '' == check(f'anchorna cutout anchors.gff a0> a2< -o {fname}')
        assert '' == check(f'anchorna go --fname {fname} --no-pbar anchors_cutout.gff')
        assert '' == check('anchorna combine anchors.gff||a1 anchors_cutout.gff -o anchors_combined.gff')
        assert read_anchors('anchors_combined.gff') == anchors

        fname = tmpdir / 'pesti_test_cutout2.sjson'
        assert '' == check(f'anchorna cutout anchors.gff a6> a10< -o {fname}')
        assert '' == check(f'anchorna go --fname {fname} --no-pbar anchors_cutout2.gff --maxshift=1000')
        assert '' == check('anchorna combine anchors.gff||a7:a10 anchors_cutout2.gff -o anchors_combined2.gff')
        assert read_anchors('anchors_combined2.gff') == load_selected_anchors('anchors.gff||a7')  # A7 is overlapping with A& and therefore not found

        # TODO: add test with missing flukes

        # check --no-remove option and --continue-with option
        assert '' == check('anchorna go --no-remove --no-pbar anchors2.gff')
        assert '' == check('anchorna go --continue-with anchors2.gff --no-pbar anchors3.gff')
        assert len(anchors) < len(read_anchors('anchors2.gff'))
        assert anchors == read_anchors('anchors3.gff')

        # test json
        json = tmpdir / 'anchors.json'
        write_json(anchors, json)
        assert load_json(json) == anchors


def test_anchorna_script_help():
    """
    Test that anchorna can be called on the command line
    """
    assert b'Find anchors' in check_output('anchorna -h'.split())
    assert b'go' in check_output('anchorna go -h'.split())
