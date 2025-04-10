# (C) 2024, Tom Eulenfeld, MIT license
"""
`.Fluke`, `.Anchor` and `.AnchorList` classes providing useful attributes and methods
"""

import collections
from copy import deepcopy
import logging
from statistics import median
from warnings import warn

from sugar.core.meta import Attr
from sugar.data import submat

log = logging.getLogger('anchorna')


def corrscore(seq1, seq2, gap='-', sm=None):
    """Similarity score between two words"""
    if sm is None:
        from sugar.data import submat
        sm = submat('blosum62')
    return sum(sm[nt1][nt2] for nt1, nt2 in zip(seq1, seq2) if nt1 != gap and nt2 != gap)


class Fluke(Attr):
    """
    A fluke is a word position on a single sequence and part of an `Anchor`

    Properties: id, seqid, score, median_score, start, stop,
    offset, word, poor.
    """
    @property
    def len(self):
        return self.stop - self.start

    @property
    def strand(self):
        return self.get('strand', '+')

    def _apply_mode(self, mode):
        """
        Transform index, allowed modes are ``'aa', 'cds', 'nt'``
        """
        fluke = self
        assert mode in ('aa', 'cds', 'nt')
        o = fluke.offset
        i1 = fluke.start
        i2 = fluke.stop
        assert i1 < i2
        if mode in ('aa', None):
            return i1, i2
        elif mode == 'cds':
            return 3 * i1, 3 * i2
        elif mode == 'nt' and fluke.strand == '-':
            return o - 3 * i2, o - 3 * i1
        else:
            return o + 3 * i1, o + 3 * i2


class Anchor(collections.UserList):
    """
    A single Anchor is a list of `Flukes <Fluke>`, one for each sequence

    Some properties: data (list of flukes), id, gseqid, guide,
    minscore (aka score), maxscore, medscore.
    """
    def __init__(self, data=None, **kw):
        super().__init__(data)
        for k, v in kw.items():
            setattr(self, k, v)

    def __str__(self):
        return self.tostr()

    @property
    def id(self):
        return f'A {self.guide.start}+{self.guide.len}'

    @property
    def strand(self):
        strands = {f.strand for f in self}
        assert len(strands) <= 1
        if len(strands) == 1:
            return strands.pop()

    @strand.setter
    def strand(self, value):
        for f in self:
            f.strand = value

    def __hash__(self):
        return hash((self.guide.start, self.guide.len, self.guide.offset))

    def todict_seqid(self):
        return {f.seqid: f for f in self}

    @property
    def sid(self):
        return self.todict_seqid()

    @property
    def seqids(self):
        return [f.seqid for f in self]

    def tostr(self, i='', verbose=False, mode='aa'):
        start, stop = self.guide._apply_mode(mode)
        poor = sum(f.poor for f in self)
        out = f'A{i} {start}+{stop-start}  minscore {self.minscore}  poor {poor}  {self.guide.word}'
        if not verbose:
            return out
        flukes = [f'  F{j} {f._apply_mode(mode)[0]}  medscore {f.median_score}  {f.word}  {f.seqid}' + '  (poor)' * f.poor
                  for j, f in enumerate(self)]
        return '\n'.join([out]+flukes)

    def sort(self, key=None, **kw):
        if key is None:
            def key(f): return (False, '') if f.seqid == self.gseqid else (f.poor, f.seqid)
        self.data = sorted(self.data, key=key, **kw)
        return self

    @property
    def guide(self):
        return self.sid[self.gseqid]

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
        return Anchor(data, gseqid=self.gseqid)

    @property
    def good_flukes(self):
        data = [f for f in self if not f.poor]
        return Anchor(data, gseqid=self.gseqid)

    def nicely_overlaps_with(self, a2):
        a1 = self.good_flukes
        a2 = a2.good_flukes
        return (max(a1.guide.start, a2.guide.start) <= min(a1.guide.stop, a2.guide.stop) and
                set(a1.seqids) == set(a2.seqids) and
                all((f1.start - f2.start == a1.guide.start - a2.guide.start and f1.poor == f2.poor) for f1, f2 in zip(a1.sort(), a2.sort())))

    def join_with(self, a2, scoring=None):
        a1 = self
        if not a1.nicely_overlaps_with(a2):
            raise ValueError('Cannot join anchors which do not overlap')
        if a1.guide.start <= a2.guide.start and a1.guide.stop >= a2.guide.stop:
            # a2 is contained in a1
            return a1
        elif a1.guide.start >= a2.guide.start and a1.guide.stop <= a2.guide.stop:
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
                    # that might only be the case for "poor" flukes below the threshold
                    # adapt the word
                    assert f1.poor and f2.poor
                    l = correctl
                    overlaplen = f1.len + f2.len - l
                    if f1.score >= f2.score:
                        stop = start + l
                    else:
                        start = stop - l
                nword = f1.word + f2.word[overlaplen:] if f1.start <= f2.start else f2.word + f1.word[overlaplen:]
                assert len(nword) == l
                fluke = Fluke(seqid=f1.seqid, score=None, median_score=None,
                              start=start, stop=stop, offset=f1.offset, word=nword,
                              poor=f1.poor)
                flukes.append(fluke)
            anchor = Anchor(flukes, gseqid=a1.gseqid)
            if scoring:
                anchor._calculate_fluke_scores(scoring)
            return anchor

    def _calculate_fluke_scores(self, scoring):
        words = {fluke.word for fluke in self if not fluke.poor}
        for fluke in self:
            scores = [corrscore(fluke.word, w, sm=submat(scoring)) for w in words]
            fluke.score = max(scores)
            fluke.median_score = median(scores)

    def contradicts(self, a2, aggressive=True):
        a1 = self
        if a1.guide.start > a2.guide.start:
            a1, a2 = a2, a1
        a1 = a1.good_flukes
        a2 = a2.good_flukes
        attr = 'stop' if aggressive else 'start'
        return not all(getattr(f1, attr) <= f2.start for f1, f2 in zip(a1.sort(), a2.sort()))

    def _convert_nt(self):
        for fluke in self:
            start, stop = fluke._apply_mode('nt')
            fluke.start, fluke.stop, fluke.offset = start, stop, 0
            fluke.strand = '+'


class AnchorList(collections.UserList):
    """
    Collection of `Anchors <Anchor>` with useful methods
    """
    def __init__(self, data=None, no_cds=False):
        super().__init__(data)
        if no_cds:
            self._no_cds = True

    def __str__(self):
        return self.tostr()

    def _repr_pretty_(self, p, cycle):
        if cycle:
            p.text('...')
        else:
            p.text(str(self))

    @property
    def no_cds(self):
        return getattr(self, '_no_cds', False)

    @no_cds.setter
    def no_cds(self, value):
        if value:
            self._no_cds = True
        else:
            raise ValueError('Not allowed to unset no_cds')

    def tostr(self, verbose=False, mode='aa'):
        return '\n'.join(a.tostr(i=i, verbose=verbose, mode=mode) for i, a in enumerate(self))

    def copy(self):
        return deepcopy(self)

    def sort(self, key=lambda a: a.guide.start, **kw):
        self.data = sorted(self, key=key, **kw)
        return self

    def merge_overlapping_anchors(self, scoring=None):
        """
        Remove overlapping anchors, step B of ``anchorna go``
        """
        self.sort()
        already_merged = set()
        ndata = []
        for i, a1 in enumerate(self.data):
            if a1 in already_merged:
                continue
            for a2 in self.data[i+1:]:
                if a2 in already_merged:
                    continue
                if a1.guide.stop < a2.guide.start:
                    break
                if a1.nicely_overlaps_with(a2):
                    a1 = a1.join_with(a2, scoring=scoring)
                    already_merged.add(a2)
            ndata.append(a1)
        self.data = ndata
        return self

    def remove_contradicting_anchors(self, aggressive=True):
        """
        Remove contradicting anchors, step C of ``anchorna go``
        """
        # The runtime of this method can be enhanced by using an interval tree or a nested containment list.
        # But this method is not the bottleneck at all.
        anchors = sorted(self, key=lambda a: a.minscore, reverse=True)
        remove_anchors = set()
        for i, a1 in enumerate(anchors):
            if a1 in remove_anchors:
                continue
            for a2 in anchors[i+1:]:
                if a2 in remove_anchors:
                    continue
                assert a1 != a2 and a1.minscore >= a2.minscore
                if a1.contradicts(a2, aggressive=aggressive):
                    log.debug(f'Remove anchor {a2.guide.start}+{a2.guide.len} with min score {a2.minscore}, '
                              f'keep anchor {a1.guide.start}+{a1.guide.len} with min score {a1.minscore}')
                    remove_anchors.add(a2)
        self.data = [anchor for anchor in self if anchor not in remove_anchors]
        return AnchorList(remove_anchors, no_cds=self.no_cds).sort()

    def convert2fts(self, **kw):
        """
        Convert anchors to `~sugar.core.fts.FeatureList` object, see `.anchors2fts()`
        """
        return anchors2fts(self, **kw)

    def _convert_nt(self):
        for anchor in self:
            anchor._convert_nt()
        self.no_cds = True

    def write(self, fname, **kw):
        """
        Write anchors to GFF file, see `.write_anchors()`
        """
        from anchorna.io import write_anchors
        return write_anchors(self, fname, **kw)


def anchors2fts(anchors, mode=None):
    """
    Convert anchors to `~sugar.core.fts.FeatureList` object

    Write some important attributes to ``Feature.meta._gff`` metadata.

    :param mode: Optionally convert anchors to ``'nt'`` or ``'cds'`` indices
    """
    from sugar import Feature, FeatureList
    fts = []
    for i, a in enumerate(anchors):
        for j, f in enumerate(a.sort()):  # fluke
            if a.gseqid == f.seqid:
                assert j == 0
                ftype = 'anchor'
                name = f'A{i}'
            else:
                assert j > 0
                ftype = 'fluke'
                name=f'A{i}_{f.seqid}'
            start, stop = f._apply_mode(mode)
            ft = Feature(ftype, start=start, stop=stop, strand=f.strand,
                         meta=dict(name=name, seqid=f.seqid, score=f.score,
                                   _gff=Attr(source='anchorna', word=f.word,
                                             median_score=f.median_score)))
            if hasattr(f, 'offset'):
                ft.meta._gff.offset = f.offset
            if f.poor:
                ft.meta._gff.poor = 1
            fts.append(ft)
    return FeatureList(fts)


def fts2anchors(fts, no_cds=False):
    """
    Convert `~sugar.core.fts.FeatureList` object to anchors

    Read some important attributes from ``Feature.meta._gff`` metadata.
    """
    anchors = []
    flukes = []
    gseqid = None
    for ft in fts:
        if ft.type == 'anchor':
            if len(flukes) > 0:
                anchors.append(Anchor(flukes, gseqid=gseqid))
                flukes = []
            gseqid = ft.seqid
            aname = ft.meta.name
        if ft.type in ('anchor', 'fluke'):
            if not ft.meta.name.startswith(aname):
                raise ValueError(f'Fluke with name {ft.meta.name} not part of anchor {aname}')
            start, stop = ft.locs.range
            fluke = Fluke(seqid=ft.seqid, score=ft.meta.score,
                          start=start, stop=stop, strand=str(ft.loc.strand),
                          word=ft.meta._gff.word, poor=hasattr(ft.meta._gff, 'poor'),
                          median_score = float(ft.meta._gff.get('median_score', 1)))
            if hasattr(ft.meta._gff, 'offset'):
                fluke.offset = ft.meta._gff.offset
            flukes.append(fluke)
        else:
            warn(f'Cannot convert feature of type {ft.type}')
    if len(flukes) > 0:
        anchors.append(Anchor(flukes, gseqid=gseqid))
    return AnchorList(anchors, no_cds=no_cds)
