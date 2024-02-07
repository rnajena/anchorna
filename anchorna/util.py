# (C) 2024, Tom Eulenfeld, MIT license

import collections
import json
import itertools
import logging
from statistics import median
import sys

from sugar import Attr

log = logging.getLogger('anchorna')


class AnchorJSONEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, (Options, Fluke, Anchor, AnchorList)):
            obj = {k: v for k, v in o.__dict__.items()}
            obj['_cls'] = type(o).__name__
            return obj
        else:
            # Let the base class default method raise the TypeError
            return json.JSONEncoder.default(self, o)


class ConfigJSONDecoder(json.JSONDecoder):
    """Decode JSON config with comments stripped"""
    def decode(self, s):
        s = '\n'.join(l.split('#', 1)[0] for l in s.split('\n'))
        return super(ConfigJSONDecoder, self).decode(s)


def json_hook(d):
    if cls := d.pop('_cls', None):
        return globals()[cls](**d)
    else:
        return d

def load_json(fname):
    if fname == '-':
        return json.loads(sys.stdin.read(), object_hook=json_hook)
    else:
        with open(fname) as f:
            return json.load(f, object_hook=json_hook)


def write_json(obj, fname):
    if fname is None:
        print(json.dumps(obj, cls=AnchorJSONEncoder))
    else:
        with open(fname, 'w') as f:
            json.dump(obj, f, cls=AnchorJSONEncoder)

class Options(Attr):
    pass


class Fluke(Attr):
    @property
    def len(self):
        return self.stop - self.start



def _apply_mode(i, o, mode, islen=False):
    if mode == 'default' or mode == 'offset' and islen:
        return i
    elif mode == 'seq' or mode == 'seqoffset' and islen:
        return 3 * i
    elif mode == 'offset':
        return i + o
    elif mode == 'seqoffset':
        return 3 * i + o
    raise ValueError('mode not allowed')



class Anchor(collections.UserList):
    def __init__(self, data=None, **kw):
        # if hasattr(data, 'meta'):
        #     meta = data.meta
        # elif 'meta' in data:
        #     meta = data['meta']
        # elif meta is None:
        #     meta = {}
        super().__init__(data)
        for k, v in kw.items():
            setattr(self, k, v)
    def __str__(self):
        return self.tostr()

    @property
    def id(self):
        return f'A {self.ref.start}+{self.ref.len}'

    def __hash__(self):
        return hash((self.ref.start, self.ref.len))

    def todict(self):
        return {f.id: f for f in self}
    @property
    def d(self):
        return self.todict()
    def tostr(self, i='', verbose=False, mode='default'):
        ind = _apply_mode(self.ref.start, self.ref.offset, mode=mode)
        len_ = _apply_mode(self.ref.len, self.ref.offset, mode=mode, islen=True)
        out = f'A{i} {ind}+{len_}  minscore {self.minscore}  misses {len(self.missing)}  {self.ref.str}'
        if not verbose:
            return out
        flukes = [f'  F{j} {_apply_mode(f.start, f.offset, mode=mode)}  score {f.score}  {f.str}  {f.id}' for j, f in enumerate(self)]
        return '\n'.join([out]+flukes)
    @property
    def ref(self):
        return self.d[self.refid]
    #worst, words
    @property
    def minscore(self):
        return min(f.score for f in self)
    @property
    def maxscore(self):
        return max(f.score for f in self)
    @property
    def medscore(self):
        return median(f.score for f in self)

    def overlaps_with(self, a2):
        key = lambda f: f.id
        a1 = self
        return (max(a1.ref.start, a2.ref.start) <= min(a1.ref.stop, a2.ref.stop) and
                {f.id for f in a1.missing} == {f.id for f in a2.missing} and
                {f.id for f in a1} == {f.id for f in a2} and
                all(f1.start - f2.start == a1.ref.start - a2.ref.start for f1, f2 in zip(sorted(a1, key=key), sorted(a2, key=key))))

    def join_with(self, a2):
        a1 = self
        if not a1.overlaps_with(a2):
            raise ValueError('Cannot join achors which do not overlap')
        if a1.ref.start <= a2.ref.start and a1.ref.stop >= a2.ref.stop:
            # a2 is contained in a1
            return a1
        elif a1.ref.start >= a2.ref.start and a1.ref.stop <= a2.ref.stop:
            # a1 is contained in a2
            return a2
        else:
            # overlap
            key = lambda f: f.id
            flukes = []
            for f1, f2 in zip(sorted(a1, key=key), sorted(a2, key=key)):
                assert f1.id == f2.id
                start = min(f1.start, f2.start)
                stop = max(f1.stop, f2.stop)
                assert f1.offset == f2.offset
                l = max(f1.stop, f2.stop) - start
                overlaplen = f1.len + f2.len - l
                score = f1.score * (f1.len - overlaplen / 2) / l + f2.score * (f2.len - overlaplen / 2) / l
                nstr = f1.str + f2.str[overlaplen:] if f1.start <= f2.start else f2.str + f1.str[overlaplen:]
                assert len(nstr) == l
                fluke = Fluke(id=f1.id, score=round(score, 3), start=start, stop=stop, offset=f1.offset, str=nstr)
                flukes.append(fluke)
            # todo recalc scores
            return Anchor(flukes, refid=a1.refid, missing=a1.missing)

    def contradicts(self, a2):
        a1 = self
        key = lambda f: f.id
        if a1.ref.start > a2.ref.start:
            a1, a2 = a2, a1
        return not all(f1.start <= f2.start for f1, f2 in zip(sorted(a1, key=key), sorted(a2, key=key)))



class AnchorList(collections.UserList):
    def __init__(self, data=None, **kw):
        # if hasattr(data, 'meta'):
        #     meta = data.meta
        # elif 'meta' in data:
        #     meta = data['meta']
        # elif meta is None:
        #     meta = {}
        super().__init__(data)
        for k, v in kw.items():
            setattr(self, k, v)
        # self.meta = Meta(meta)
    def __str__(self):
        return self.tostr()

    def _repr_pretty_(self, p, cycle):
        if cycle:
            p.text('...')
        else:
            p.text(str(self))

    def tostr(self, verbose=False, mode='default'):
        return '\n'.join(a.tostr(i=i, verbose=verbose, mode=mode) for i, a in enumerate(self))

    def merge_neighbor_anchors(self):
        while True:
            for a1, a2 in itertools.combinations(sorted(self, key=lambda a: a.minscore, reverse=True), 2):
                if a1.overlaps_with(a2):
                    break
            else:
                break
            # merge a1 and a2
            self.data.remove(a1)
            self.data.remove(a2)
            a = a1.join_with(a2)
            self.data.append(a)
        return self

    def remove_contradicting_anchors(self):
        while True:
            for a1, a2 in itertools.combinations(self, 2):
                if a1.contradicts(a2):
                    break
            else:
                break
            # remove one of a1, a2
            if a1.minscore < a2.minscore:
                a1, a2 = a2, a1
            log.debug(f'Remove anchor {a2.ref.start}+{a2.ref.len} with min score {a2.minscore}, '
                      f'keep anchor {a1.ref.start}+{a1.ref.len} with min score {a1.minscore}')
            self.data.remove(a2)
        return self


import matplotlib.colors as mcolors
from matplotlib.colors import to_hex, to_rgb

def make_rgb_transparent(rgb, bg_rgb, alpha):
    return [alpha * c1 + (1 - alpha) * c2
            for (c1, c2) in zip(to_rgb(rgb), to_rgb(bg_rgb))]

cols = list(mcolors.TABLEAU_COLORS.values())
to_hex(cols[0], 0.3)

def jalview_features(anchors, mode='default'):
    assert mode in ('default', 'offset', 'seq', 'seqoffset')
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
            c = to_hex(make_rgb_transparent(cols[k % len(cols)], 'white', score_)).strip('#')

            al = f'anchor{k}_s{f.score}'
            miss = f'misses:{len(a.missing)}'
            header.append(
                f'{al}\t{c}\n'
                )
            i = _apply_mode(f.start, f.offset, mode)
            j = _apply_mode(f.stop, f.offset, mode)
            content.append(
                f'{a.ref.str[:5]} w{wlen} {miss}\t{f.id}\t-1\t{i+1}\t{j}\t{al}\n')
                 #f'{a.ref.str[:5]} w{wlen} {miss}\t{f.id}\t-1\t{i+1}\t{i+wlen}\t{al}\n' if mode == 'default' else
                 #f'{a.ref.str[:5]} w{wlen} {miss}\t{f.id}\t-1\t{i+o+1}\t{i+o+wlen}\t{al}\n' if mode == 'offset' else
                 #f'{a.ref.str[:5]} w{wlen} {miss}\t{f.id}\t-1\t{3*i+1}\t{3*i+3*wlen}\t{al}\n' if mode == 'seq' else
                 #f'{a.ref.str[:5]} w{wlen} {miss}\t{f.id}\t-1\t{3*i+o+1}\t{3*i+o+3*wlen}\t{al}\n' if mode == 'seqoffset' else
                 #None)
    header.append('\nSTARTFILTERS\nENDFILTERS\n\n')
    return ''.join(header) + ''.join(content)





# def write_anchortxt(anchors, fname):
#     with open(fname, 'w') as f:
#         print('# winlen score index(aa) index(nt) aa nt', file=f)
#         for a in anchors:
#             print(*a[:-1], file=f)


# def __write_jalview_features_simple(anchors, fnameaa, fnament, fnamecds):
#     header = ('anchor\tlabel\n'
#               'anchorsim\tscore|bad9f7|105191|absolute|0.7|1\n\n'
#               'STARTFILTERS\n'
#               'ENDFILTERS\n\n')
#     anchors = sorted(anchors, key=lambda a: a.ref.start)
#     contentnt = []
#     contentcds = []
#     contentaa = []
#     #for winlen, *_, aaref, r in anchors:
#     for a in anchors:
#         #for seqid, (sim, i, j, _) in r.items():
#         for f in a:
#             # anchor2	C_CSFV_KC533775	-1	130	150	anchorsim	1.0
#             # anhcor100	D_BDV_NC_003679	-1	10	20	anchorxx
#             wlen = a.ref.len
#             i = f.start
#             j = f.seqstart
#             contentaa.append(
#                   #f'anchorsim\t{seqid}\t-1\t{i+1}\t{i+winlen}\tanchorsim\t{sim}\n'
#                  f'{a.ref.str[:5]} w{wlen}\t{f.id}\t-1\t{i+1}\t{i+wlen}\tanchor\n')
#             contentnt.append(
#                   #f'anchorsim\t{seqid}\t-1\t{j+1}\t{j+3*winlen}\tanchorsim\t{sim}\n'
#                  f'{a.ref.str[:5]} w{wlen}\t{f.id}\t-1\t{j+1}\t{j+3*wlen}\tanchor\n')
#             contentcds.append(
#                   #f'anchorsim\t{seqid}\t-1\t{j+1}\t{j+3*winlen}\tanchorsim\t{sim}\n'
#                  f'{a.ref.str[:5]} w{wlen}\t{f.id}\t-1\t{3*i+1}\t{3*i+3*wlen}\tanchor\n')

#     if fnameaa:
#         with open(fnameaa, 'w') as f:
#             f.write(header + ''.join(contentaa))
#     if fnament:
#         with open(fnament, 'w') as f:
#             f.write(header + ''.join(contentnt))
#     if fnamecds:
#         with open(fnamecds, 'w') as f:
#             f.write(header + ''.join(contentcds))


# def __write_jalview_features(anchors, fnameaa, fnament, fnamecds):
#     anchors = sorted(anchors, key=lambda a: a.ref.start)
#     contentnt = []
#     contentcds = []
#     contentaa = []
#     header = []
#     #for winlen, *_, aaref, r in anchors:
#     for k, a in enumerate(anchors):
#         #for seqid, (score, i, j, _) in r.items():
#         for j, f in enumerate(a):
#             # anchor2	C_CSFV_KC533775	-1	130	150	anchorsim	1.0
#             # anhcor100	D_BDV_NC_003679	-1	10	20	anchorxx
#             wlen = a.ref.len
#             i = f.start
#             j = f.seqstart
#             score_ = f.score / a.maxscore
#             c = to_hex(make_rgb_transparent(cols[k % len(cols)], 'white', score_)).strip('#')

#             al = f'anchor{k}_s{f.score}'
#             miss = f'misses:{len(a.missing)}'
#             header.append(
#                 f'{al}\t{c}\n'
#                 )
#             contentaa.append(
#                  f'{a.ref.str[:5]} w{wlen} {miss}\t{f.id}\t-1\t{i+1}\t{i+wlen}\t{al}\n')
#             contentnt.append(
#                  f'{a.ref.str[:5]} w{wlen} {miss}\t{f.id}\t-1\t{j+1}\t{j+3*wlen}\t{al}\n')
#             contentcds.append(
#                  f'{a.ref.str[:5]} w{wlen} {miss}\t{f.id}\t-1\t{3*i+1}\t{3*i+3*wlen}\t{al}\n')

#     header.extend(['STARTFILTERS\n', 'ENDFILTERS\n\n'])
#     header = ''.join(header)

#     if fnameaa:
#         with open(fnameaa, 'w') as f:
#             f.write(header + ''.join(contentaa))
#     if fnament:
#         with open(fnament, 'w') as f:
#             f.write(header + ''.join(contentnt))
#     if fnamecds:
#         with open(fnamecds, 'w') as f:
#             f.write(header + ''.join(contentcds))


# def __write_jalview_features(anchors, fnameaa, fnament, fnamecds):
#     anchors = sorted(anchors, key=lambda a: a.ref.start)
#     contentnt = []
#     contentcds = []
#     contentaa = []
#     header = []
#     #for winlen, *_, aaref, r in anchors:
#     for k, a in enumerate(anchors):
#         #for seqid, (score, i, j, _) in r.items():
#         for j, f in enumerate(a):
#             # anchor2	C_CSFV_KC533775	-1	130	150	anchorsim	1.0
#             # anhcor100	D_BDV_NC_003679	-1	10	20	anchorxx
#             wlen = a.ref.len
#             i = f.start
#             j = f.seqstart
#             score_ = f.score / a.maxscore
#             c = to_hex(make_rgb_transparent(cols[k % len(cols)], 'white', score_)).strip('#')

#             al = f'anchor{k}_s{f.score}'
#             miss = f'misses:{len(a.missing)}'
#             header.append(
#                 f'{al}\t{c}\n'
#                 )
#             contentaa.append(
#                  f'{a.ref.str[:5]} w{wlen} {miss}\t{f.id}\t-1\t{i+1}\t{i+wlen}\t{al}\n')
#             contentnt.append(
#                  f'{a.ref.str[:5]} w{wlen} {miss}\t{f.id}\t-1\t{j+1}\t{j+3*wlen}\t{al}\n')
#             contentcds.append(
#                  f'{a.ref.str[:5]} w{wlen} {miss}\t{f.id}\t-1\t{3*i+1}\t{3*i+3*wlen}\t{al}\n')

#     header.extend(['STARTFILTERS\n', 'ENDFILTERS\n\n'])
#     header = ''.join(header)

#     if fnameaa:
#         with open(fnameaa, 'w') as f:
#             f.write(header + ''.join(contentaa))
#     if fnament:
#         with open(fnament, 'w') as f:
#             f.write(header + ''.join(contentnt))
#     if fnamecds:
#         with open(fnamecds, 'w') as f:
#             f.write(header + ''.join(contentcds))
