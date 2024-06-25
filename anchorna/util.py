# (C) 2024, Tom Eulenfeld, MIT license

import collections
import itertools
import logging
from statistics import median
from warnings import warn

from sugar import Attr
from sugar.data import submat


log = logging.getLogger('anchorna')


def corrscore(seq1, seq2, gap='-', sm=submat('blosum62')):
    return sum(sm[nt1][nt2] for nt1, nt2 in zip(seq1, seq2) if nt1 != gap and nt2 != gap)


class Options(Attr):
    pass


class Fluke(Attr):
    @property
    def len(self):
        return self.stop - self.start


def _apply_mode(i, o, mode, islen=False):
    if mode == 'aa':
        return i
    elif mode == 'cds' or mode == 'seq' and islen:
        return 3 * i
    elif mode == 'seq':
        return 3 * i + o
    raise ValueError('mode not allowed')


class Anchor(collections.UserList):
    def __init__(self, data=None, **kw):
        super().__init__(data)
        for k, v in kw.items():
            setattr(self, k, v)

    def __str__(self):
        return self.tostr()

    @property
    def id(self):
        return f'A {self.ref.start}+{self.ref.len}'

    def __hash__(self):
        return hash((self.ref.start, self.ref.len, self.ref.offset))

    def todict(self):
        return {f.seqid: f for f in self}

    @property
    def d(self):
        return self.todict()

    def tostr(self, i='', verbose=False, mode='aa'):
        ind = _apply_mode(self.ref.start, self.ref.offset, mode=mode)
        len_ = _apply_mode(self.ref.len, self.ref.offset,
                           mode=mode, islen=True)
        poor = sum(f.poor for f in self)
        out = f'A{i} {ind}+{len_}  minscore {self.minscore}  poor {poor}  {self.ref.word}'
        if not verbose:
            return out
        flukes = [f'  F{j} {_apply_mode(f.start, f.offset, mode=mode)}  medscore {f.median_score}  {f.word}  {f.seqid}' + '  (poor)' * f.poor
                  for j, f in enumerate(self)]
        return '\n'.join([out]+flukes)

    def sort(self, key=None, **kw):
        if key is None:
            def key(f): return (0, '') if f.seqid == self.refid else (f.poor, f.seqid)
        self.data = sorted(self.data, key=key, **kw)
        return self

    @property
    def ref(self):
        return self.d[self.refid]

    @property
    def minscore(self):
        return min(f.median_score for f in self)

    @property
    def maxscore(self):
        return max(f.median_score for f in self)

    @property
    def medscore(self):
        return median(f.median_score for f in self)

    @property
    def poor_flukes(self):
        data = [f for f in self if f.poor]
        return Anchor(data, refid=self.refid)

    @property
    def good_flukes(self):
        data = [f for f in self if not f.poor]
        return Anchor(data, refid=self.refid)

    def overlaps_in_a_simple_way_with(self, a2):
        a1 = self.good_flukes
        a2 = a2.good_flukes
        return (max(a1.ref.start, a2.ref.start) <= min(a1.ref.stop, a2.ref.stop) and
                all((f1.start - f2.start == a1.ref.start - a2.ref.start and f1.poor == f2.poor) for f1, f2 in zip(a1.sort(), a2.sort())))

    def join_with(self, a2):
        a1 = self
        if not a1.overlaps_in_a_simple_way_with(a2):
            raise ValueError('Cannot join anchors which do not overlap')
        if a1.ref.start <= a2.ref.start and a1.ref.stop >= a2.ref.stop:
            # a2 is contained in a1
            return a1
        elif a1.ref.start >= a2.ref.start and a1.ref.stop <= a2.ref.stop:
            # a1 is contained in a2
            return a2
        else:
            # overlap
            flukes = []
            correctl = None
            for f1, f2 in zip(a1.sort(), a2.sort()):
                assert f1.seqid == f2.seqid
                assert f1.poor == f2.poor
                assert f1.offset == f2.offset
                if f1.start > f2.start:
                    f1, f2 = f2, f1
                start = f1.start
                stop = f2.stop
                l = stop - start
                overlaplen = f1.len + f2.len - l
                if not f1.poor:
                    if correctl is None:
                        correctl = l
                    assert correctl == l
                if correctl != l:
                    # there is a gap between the flukes!
                    # that might only be the case for "poor" flukes below the thresshold
                    # adapt the word
                    assert f1.poor and f2.poor
                    l = correctl
                    overlaplen = f1.len + f2.len - l
                    if f1.score >= f2.score:
                        stop = start + l
                    else:
                        start = stop - l
                # old way to calculate the score: estimation in the following way, we do not recalculate it
                #score = f1.score * (1 - overlaplen / 2 / f1.len) + f2.score * (1 - overlaplen / 2 / f2.len)
                nword = f1.word + f2.word[overlaplen:] if f1.start <= f2.start else f2.word + f1.word[overlaplen:]
                assert len(nword) == l
                fluke = Fluke(seqid=f1.seqid, score=None,
                              start=start, stop=stop, offset=f1.offset, word=nword,
                              poor=f1.poor)
                flukes.append(fluke)
            anchor = Anchor(flukes, refid=a1.refid)
            anchor._calculate_fluke_scores()
            return anchor


    def _calculate_fluke_scores(self):
        words = {fluke.word for fluke in self if not fluke.poor}
        for fluke in self:
            scores = [corrscore(fluke.word, w) for w in words]
            fluke.score = max(scores)
            fluke.median_score = median(scores)


    def contradicts(self, a2):
        a1 = self
        if a1.ref.start > a2.ref.start:
            a1, a2 = a2, a1
        a1 = a1.good_flukes
        a2 = a2.good_flukes
        return not all(f1.start <= f2.start for f1, f2 in zip(a1.sort(), a2.sort()))


class AnchorList(collections.UserList):
    def __init__(self, data=None):
        super().__init__(data)

    def __str__(self):
        return self.tostr()

    def _repr_pretty_(self, p, cycle):
        if cycle:
            p.text('...')
        else:
            p.text(str(self))

    def tostr(self, verbose=False, mode='aa'):
        return '\n'.join(a.tostr(i=i, verbose=verbose, mode=mode) for i, a in enumerate(self))

    def sort(self, key=lambda a: a.ref.start, **kw):
        self.data = sorted(self, key=key, **kw)
        return self

    def merge_neighbor_anchors(self):
        while True:
            for a1, a2 in itertools.combinations(sorted(self, key=lambda a: a.minscore, reverse=True), 2):
                if a1.overlaps_in_a_simple_way_with(a2):
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
        anchors = sorted(self, key=lambda a: a.minscore, reverse=True)
        removed_anchors = []
        while len(anchors) > 0:
            a1 = anchors.pop(0)
            for a2 in anchors:
                assert a1 != a2 and a1.minscore >= a2.minscore
                if a1.contradicts(a2):
                    log.debug(f'Remove anchor {a2.ref.start}+{a2.ref.len} with min score {a2.minscore}, '
                              f'keep anchor {a1.ref.start}+{a1.ref.len} with min score {a1.minscore}')
                    self.data.remove(a2)
                    anchors.remove(a2)
                    removed_anchors.append(a2)
        return AnchorList(removed_anchors).sort()

    def convert2fts(self):
        return anchors2fts(self)

    def write(self, fname):
        from anchorna.io import write_anchors
        return write_anchors(self, fname)


def anchors2fts(anchors):
    """
    Convert anchors to sugar.FeatureList object

    Write some important attributes to Feature.meta._gff metadata.
    """
    from sugar import Feature, FeatureList
    fts = []
    for i, a in enumerate(anchors):
        for j, f in enumerate(a.sort()):  # fluke
            if a.refid == f.seqid:
                assert j == 0
                ftype = 'anchor'
                name = f'A{i}'
            else:
                assert j > 0
                ftype = 'fluke'
                name=f'A{i}_{f.seqid}'
            ft = Feature(ftype, start=f.start, stop=f.stop,
                         meta=dict(name=name, seqid=f.seqid, score=f.score,
                                   _gff=Attr(source='anchorna', word=f.word,
                                             median_score=f.median_score)))
            if hasattr(f, 'offset'):
                ft.meta._gff.offset = f.offset
            if f.poor:
                ft.meta._gff.poor = 1
            fts.append(ft)
    return FeatureList(fts)


def fts2anchors(fts):
    """
    Convert sugar.FeatureList object to anchors

    Read some important attributes from Feature.meta._gff metadata.
    """
    anchors = []
    flukes = []
    refid = None
    for ft in fts:
        if ft.type == 'anchor':
            if len(flukes) > 0:
                anchors.append(Anchor(flukes, refid=refid))
                flukes = []
            refid = ft.seqid
            aname = ft.meta.name
        if ft.type in ('anchor', 'fluke'):
            if not ft.meta.name.startswith(aname):
                raise ValueError(f'Fluke with name {ft.meta.name} not part of anchor {aname}')
            start, stop = ft.loc_range
            fluke = Fluke(seqid=ft.seqid, score=ft.meta.score, start=start, stop=stop,
                          word=ft.meta._gff.word, poor=hasattr(ft.meta._gff, 'poor'),
                          median_score = float(ft.meta._gff.get('median_score', 1)))
            if hasattr(ft.meta._gff, 'offset'):
                fluke.offset = ft.meta._gff.offset
            flukes.append(fluke)
        else:
            warn(f'Cannot convert feature of type {ft.type}')
    if len(flukes) > 0:
        anchors.append(Anchor(flukes, refid=refid))
    return AnchorList(anchors)


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
