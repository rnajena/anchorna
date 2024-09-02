# (C) 2024, Tom Eulenfeld, MIT license
import json
import logging
import sys

import matplotlib.colors as mcolors
from matplotlib.colors import to_hex, to_rgb
from sugar import read_fts

from anchorna.util import _apply_mode, fts2anchors, Anchor, AnchorList, Fluke


log = logging.getLogger('anchorna')


def write_anchors(anchors, fname, mode=None):
    """
    Write anchors to GFF file

    Offsets are stored as comments.

    :param str mode: if specified, transform indices and export to GFF file
       The result is a file which should not be read in again with anchorna.
    """
    from anchorna import __version__
    fts = anchors.convert2fts()
    offsets = {ft.seqid: ft.meta._gff.pop('offset') for ft in fts}
    offsets_header = ''.join(f'#offset {seqid} {offset}\n' for seqid, offset in offsets.items())
    if mode is None:
        header_cds = '#no_cds\n' if anchors.no_cds else (
            '# Indices are given for amino acids, the offset specifies the offset of index 0 from this file\n'
            '# to the beginning of the original sequence (in nucleotides).\n'
            )

        header = (
            f'#AnchoRNA anchor file\n'
            f'# written with AnchoRNA v{__version__}\n' + header_cds + offsets_header)
    else:
        header = f'# Anchors exported by AnchoRNA v{__version__} with mode {mode}\n'
    if fname is None:
        try:
            print(fts.tofmtstr('gff', header=header))
        except BrokenPipeError:
            pass
    else:
        fts.write(fname, 'gff', header=header)


def read_anchors(fname, check_header=True):
    """
    Read anchors from GFF file

    Offsets are restored from comments.
    """
    comments = []
    fts = read_fts(fname, 'gff', comments=comments)
    offsets = {}
    no_cds = False
    for i, line in enumerate(comments):
        if i == 0:
            if check_header and not line.startswith('##gff-version 3'):
                raise IOError(f'{fname} not a valid GFF file')
        elif i == 1:
            if check_header and not line.startswith('#AnchoRNA'):
                raise IOError(f'{fname} not a valid anchor file')
        elif line.startswith('#no_cds'):
            no_cds = True
        elif line.startswith('#offset'):
            seqid, offset = line.split()[1:]
            offsets[seqid] = int(offset)
        elif line.startswith('#'):
            continue
    for ft in fts:
        ft.meta._gff.offset = offsets[ft.seqid]
    return fts2anchors(fts, no_cds=no_cds)


def _parse_selection(anchors, selection):
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


def load_selected_anchors(fname):
    """
    Read anchors and select or remove some of them
    """
    if '|' not in fname:
        return read_anchors(fname)
    fname, selection = fname.split('|', 1)
    anchors = read_anchors(fname.strip())
    selection = selection.lower()
    if '|' not in selection:
        return _parse_selection(anchors, selection)
    selection, rem = selection.split('|')
    anchors_remove = set(_parse_selection(anchors, rem))
    anchors = _parse_selection(anchors, selection) if selection.strip() else anchors
    anchors.data = [a for a in anchors if a not in anchors_remove]
    return anchors


class _AnchorJSONEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, (Anchor, AnchorList, Fluke)):
            obj = {k: v for k, v in o.__dict__.items()}
            obj['_cls'] = type(o).__name__
            return obj
        else:
            # Let the base class default method raise the TypeError
            return json.JSONEncoder.default(self, o)


def _json_hook(d):
    if cls := d.pop('_cls', None):
        return globals()[cls](**d)
    else:
        return d


def load_json(fname):
    """
    Load anchors from JSON file
    """

    if fname in ('-', None):
        return json.loads(sys.stdin.read(), object_hook=_json_hook)
    else:
        with open(fname) as f:
            return json.load(f, object_hook=_json_hook)


def write_json(obj, fname):
    """
    Write anchors to JSON file
    """
    if fname in ('-', None):
        print(json.dumps(obj, cls=_AnchorJSONEncoder))
    else:
        with open(fname, 'w') as f:
            json.dump(obj, f, cls=_AnchorJSONEncoder)


def _make_rgb_transparent(rgb, bg_rgb, alpha):
    return [alpha * c1 + (1 - alpha) * c2
            for (c1, c2) in zip(to_rgb(rgb), to_rgb(bg_rgb))]


cols = list(mcolors.TABLEAU_COLORS.values())
to_hex(cols[0], 0.3)

def jalview_features(anchors, mode='aa', score_use_fluke=None):
    """
    Write anchors to Jalview feature file
    """
    assert mode in ('nt', 'cds', 'aa')
    anchors = sorted(anchors, key=lambda a: a.guide.start)
    content = []
    header = []
    for k, a in enumerate(anchors):
        poor = sum(f.poor for f in a)
        for j, f in enumerate(a):
            if score_use_fluke is not None and f.score < score_use_fluke:
                continue
            # anchor2	C_CSFV_KC533775	-1	130	150	anchorsim	1.0
            # anhcor100	D_BDV_NC_003679	-1	10	20	anchorxx
            wlen = a.guide.len
            score_ = max(1, f.median_score) / a.maxscore
            c = to_hex(_make_rgb_transparent(cols[k % len(cols)], 'white', score_)).strip('#')
            al = f'anchor{k}_s{f.median_score}'
            header.append(
                f'{al}\t{c}\n'
                )
            i = _apply_mode(f.start, f.offset, mode)
            j = _apply_mode(f.stop, f.offset, mode)
            content.append(f'{f.word[:5]} w{wlen} poor:{poor}\t{f.seqid}\t-1\t{i+1}\t{j}\t{al}\n')
    header.append('\nSTARTFILTERS\nENDFILTERS\n\n')
    return ''.join(header) + ''.join(content)
