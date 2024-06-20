# (C) 2024, Tom Eulenfeld, MIT license
import json
import logging
import sys
from warnings import warn

import matplotlib.colors as mcolors
from matplotlib.colors import to_hex, to_rgb
from sugar import read_fts

from anchorna.util import _apply_mode, fts2anchors, Anchor, AnchorList, Fluke, Options


log = logging.getLogger('anchorna')


def write_anchors(anchors, fname):
    """
    Write anchors to GFF file

    Offsets are stored as comments.
    """
    from anchorna import __version__
    fts = anchors.convert2fts()
    offsets = {ft.seqid: ft.meta._gff.pop('offset') for ft in fts}
    offsets_header = ''.join(f'#offset {seqid} {offset}\n' for seqid, offset in offsets.items())
    header = f'#AnchoRNA anchor file\n#written with AnchoRNA v{__version__}\n' + offsets_header
    if fname is None:
        try:
            print(fts.tofmtstr('gff', header=header))
        except BrokenPipeError:
            pass
    else:
        fts.write(fname, 'gff', header=header)


def read_anchors(fname):
    """
    Read anchors from GFF file

    Offsets are restored from comments.
    """
    comments = []
    fts = read_fts(fname, 'gff', comments=comments)
    offsets = {}
    for i, line in enumerate(comments):
        if i == 0:
            if not line.startswith('##gff-version 3'):
                warn(f'{fname} not a valid GFF file')
        elif i == 1:
            if not line.startswith('#AnchoRNA'):
                warn(f'{fname} not a valid anchor file')
        elif line.startswith('#offset'):
            seqid, offset = line.split()[1:]
            offsets[seqid] = int(offset)
        elif line.startswith('#'):
            continue
        else:
            break
    for ft in fts:
        ft.meta._gff.offset = offsets[ft.seqid]
    return fts2anchors(fts)


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
        if isinstance(o, (Anchor, AnchorList, Fluke, Options)):
            obj = {k: v for k, v in o.__dict__.items()}
            obj['_cls'] = type(o).__name__
            return obj
        else:
            # Let the base class default method raise the TypeError
            return json.JSONEncoder.default(self, o)


class _ConfigJSONDecoder(json.JSONDecoder):
    """Decode JSON config with comments stripped"""
    def decode(self, s):
        s = '\n'.join(l.split('#', 1)[0] for l in s.split('\n'))
        return super(_ConfigJSONDecoder, self).decode(s)


def _json_hook(d):
    if cls := d.pop('_cls', None):
        return globals()[cls](**d)
    else:
        return d


def load_json(fname):
    """
    Load anchors from JSON file
    """

    if fname == '-':
        return json.loads(sys.stdin.read(), object_hook=_json_hook)
    else:
        with open(fname) as f:
            return json.load(f, object_hook=_json_hook)


def write_json(obj, fname):
    """
    Write anchors to JSON file
    """
    if fname is None:
        print(json.dumps(obj, cls=_AnchorJSONEncoder))
    else:
        with open(fname, 'w') as f:
            json.dump(obj, f, cls=_AnchorJSONEncoder)


def _make_rgb_transparent(rgb, bg_rgb, alpha):
    return [alpha * c1 + (1 - alpha) * c2
            for (c1, c2) in zip(to_rgb(rgb), to_rgb(bg_rgb))]


cols = list(mcolors.TABLEAU_COLORS.values())
to_hex(cols[0], 0.3)

def jalview_features(anchors, mode='aa'):
    """
    Write anchors to Jalview feature file
    """
    assert mode in ('seq', 'cds', 'aa')
    anchors = sorted(anchors, key=lambda a: a.ref.start)
    content = []
    header = []
    #for winlen, *_, aaref, r in anchors:
    for k, a in enumerate(anchors):
        #for seqid, (score, i, j, _) in r.items():
        for j, f in enumerate(a):
            # anchor2	C_CSFV_KC533775	-1	130	150	anchorsim	1.0
            # anhcor100	D_BDV_NC_003679	-1	10	20	anchorxx
            wlen = a.ref.len
            score_ = f.score / a.maxscore
            c = to_hex(_make_rgb_transparent(cols[k % len(cols)], 'white', score_)).strip('#')

            al = f'anchor{k}_s{f.score}'
            miss = f'misses:{len(a.missing)}'
            header.append(
                f'{al}\t{c}\n'
                )
            i = _apply_mode(f.start, f.offset, mode)
            j = _apply_mode(f.stop, f.offset, mode)
            content.append(
                f'{a.ref.word[:5]} w{wlen} {miss}\t{f.seqid}\t-1\t{i+1}\t{j}\t{al}\n')
                 #f'{a.ref.str[:5]} w{wlen} {miss}\t{f.id}\t-1\t{i+1}\t{i+wlen}\t{al}\n' if mode == 'default' else
                 #f'{a.ref.str[:5]} w{wlen} {miss}\t{f.id}\t-1\t{i+o+1}\t{i+o+wlen}\t{al}\n' if mode == 'offset' else
                 #f'{a.ref.str[:5]} w{wlen} {miss}\t{f.id}\t-1\t{3*i+1}\t{3*i+3*wlen}\t{al}\n' if mode == 'seq' else
                 #f'{a.ref.str[:5]} w{wlen} {miss}\t{f.id}\t-1\t{3*i+o+1}\t{3*i+o+3*wlen}\t{al}\n' if mode == 'seqoffset' else
                 #None)
    header.append('\nSTARTFILTERS\nENDFILTERS\n\n')
    return ''.join(header) + ''.join(content)
