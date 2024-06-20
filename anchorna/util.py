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
        super().__init__([fluke for fluke in data if fluke.above_thres])
        for k, v in kw.items():
            setattr(self, k, v)
        self.missing = [fluke for fluke in data if not fluke.above_thres]

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
        out = f'A{i} {ind}+{len_}  minscore {self.minscore}  misses {len(self.missing)}  {self.ref.word}'
        if not verbose:
            return out
        flukes = [f'  F{j} {_apply_mode(f.start, f.offset, mode=mode)}  score {f.score}  {f.word}  {f.seqid}' + '  (below thresshold)' * (not f.above_thres)
                  for j, f in enumerate(self)]
        flukesm = [f'  -- {_apply_mode(f.start, f.offset, mode=mode)}  score {f.score}  {f.word}  {f.seqid}' + '  (below thresshold)' * (not f.above_thres)
                  for j, f in enumerate(self.missing)]
        return '\n'.join([out]+flukes + flukesm)

    @property
    def ref(self):
        return self.d[self.refid]

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
        def key(f): return '' if f.seqid == a1.refid else f.seqid
        a1 = self
        return (max(a1.ref.start, a2.ref.start) <= min(a1.ref.stop, a2.ref.stop) and
                {f.seqid for f in a1.missing} == {f.seqid for f in a2.missing} and
                {f.seqid for f in a1} == {f.seqid for f in a2} and
                all(f1.start - f2.start == a1.ref.start - a2.ref.start for f1, f2 in zip(sorted(a1, key=key), sorted(a2, key=key))))

    def join_with(self, a2):
        a1 = self
        if not a1.overlaps_with(a2):
            raise ValueError('Cannot join anchors which do not overlap')
        if a1.ref.start <= a2.ref.start and a1.ref.stop >= a2.ref.stop:
            # a2 is contained in a1
            return a1
        elif a1.ref.start >= a2.ref.start and a1.ref.stop <= a2.ref.stop:
            # a1 is contained in a2
            return a2
        else:
            # overlap
            def key(f): return not f.above_thres, '' if f.seqid == a1.refid else f.seqid
            flukes = []
            correctl = None
            for f1, f2 in zip(sorted(a1.data + a1.missing, key=key), sorted(a2.data + a2.missing, key=key)):
                assert f1.seqid == f2.seqid
                assert f1.above_thres == f2.above_thres
                assert f1.offset == f2.offset
                start = min(f1.start, f2.start)
                stop = max(f1.stop, f2.stop)
                l = max(f1.stop, f2.stop) - start
                if f1.above_thres:
                    if correctl is None:
                        correctl = l
                    assert correctl == l
                overlaplen = f1.len + f2.len - l
                if correctl != l:
                    # there is a gap between the flukes!
                    # that might only be the case for "missing" flukes below the thresshold
                    # adapt the word
                    assert not f1.above_thres and not f2.above_thres
                    l = correctl
                    stop = start + l
                    overlaplen = f1.len + f2.len - l
                # old way to calculate the score: estimation in the following way, we do not recalculate it
                #score = f1.score * (1 - overlaplen / 2 / f1.len) + f2.score * (1 - overlaplen / 2 / f2.len)
                nword = f1.word + f2.word[overlaplen:] if f1.start <= f2.start else f2.word + f1.word[overlaplen:]
                assert len(nword) == l
                fluke = Fluke(seqid=f1.seqid, score=None,
                              start=start, stop=stop, offset=f1.offset, word=nword,
                              above_thres=f1.above_thres or f2.above_thres)
                flukes.append(fluke)
            anchor = Anchor(flukes, refid=a1.refid)
            anchor._calculate_scores()
            return anchor


    def _calculate_scores(self):
        words = {fluke.word for fluke in self if fluke.above_thres}
        for fluke in self.data + self.missing:
            fluke.score = median(corrscore(fluke.word, w) for w in words)


    def contradicts(self, a2):
        a1 = self
        def key(f): return '' if f.seqid == a1.refid else f.seqid
        if a1.ref.start > a2.ref.start:
            a1, a2 = a2, a1
        return not all(f1.start <= f2.start for f1, f2 in zip(sorted(a1, key=key), sorted(a2, key=key)))


class AnchorList(collections.UserList):
    def __init__(self, data=None, **kw):
        super().__init__(data)
        for k, v in kw.items():
            setattr(self, k, v)

    def __str__(self):
        return self.tostr()

    def _repr_pretty_(self, p, cycle):
        if cycle:
            p.text('...')
        else:
            p.text(str(self))

    def tostr(self, verbose=False, mode='aa'):
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
        removed_anchors = []
        a1s = sorted(self, key=lambda a: a.minscore, reverse=True)
        a2s = sorted(self, key=lambda a: a.minscore)
        while True:
            for a1, a2 in itertools.product(a1s, a2s):
                if a1 == a2 or a1.minscore < a2.minscore:
                    continue
                if a1.contradicts(a2):
                    break
            else:
                break
            # remove a2
            # remove one of a1, a2
            # if a1.minscore < a2.minscore:
            #     a1, a2 = a2, a1
            log.debug(f'Remove anchor {a2.ref.start}+{a2.ref.len} with min score {a2.minscore}, '
                      f'keep anchor {a1.ref.start}+{a1.ref.len} with min score {a1.minscore}')
            self.data.remove(a2)
            a2s.remove(a2)
            a1s.remove(a2)
            removed_anchors.append(a2)
        return AnchorList(removed_anchors)

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
        def key(f): return '' if f.seqid == a.refid else f.seqid
        data = sorted(a.data + a.missing, key=key)
        for j, f in enumerate(data):  # fluke
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
                                   _gff=Attr(source='anchorna', word=f.word)))
            if hasattr(f, 'offset'):
                ft.meta._gff.offset = f.offset
            if not f.above_thres:
                ft.meta._gff.below_thres = 1
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
                          word=ft.meta._gff.word, above_thres=not hasattr(ft.meta._gff, 'below_thres'))
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
