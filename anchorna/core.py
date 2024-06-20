# (C) 2023, Tom Eulenfeld, MIT license

import logging
from statistics import median
from warnings import warn

from sugar import BioBasket
from sugar.data import submat
from tqdm import tqdm

from anchorna.util import _apply_mode, corrscore, Anchor, AnchorList, Fluke


log = logging.getLogger('anchorna')


def maxes(a, key=None, default=None):
    """
    Like max function, but returns list of *all* maximal values
    """
    if len(a) == 0:
        return default
    if key is None:
        key = lambda x: x
    kmax = key(a[0])
    max_list = []
    for s in a:
        k = key(s)
        if k > kmax:
            kmax = k
            max_list = [s]
        elif k == kmax:
            max_list.append(s)
    return max_list


def shift_and_find_best_word(seq, words, starti, w, sm, maxshift=None, maxshift_right=None):
    """
    Find position of the most similar word in a sequence compared to a set of words

    :return: tuple with similarity and index
    """
    seq = str(seq)
    if maxshift is None:
        maxshift = len(seq)
    if maxshift_right is None:
        maxshift_right = maxshift
    siminds = maxes([(max(corrscore(seq[i:i+w], word, sm=sm) for word in words), i)
                      for i in range(max(0, starti-maxshift_right),
                                     min(len(seq)-w, starti+maxshift+1))],
                     default=(0, None))
    if len(siminds) > 1:
        log.warning(f'multiple max values exist: {siminds=}')
    sim, ind = siminds[0]
    return sim, ind


def anchor_for_words(i, words, aas, w, refid, maxshift, thr_score, thr_quota_score, thr_quota,
                     scoring, method='default'):
    """
    Find the most similar not yet found word in a list of sequences and add it to words set

    If no good word is found anymore, return a single anchor.
    """
    assert thr_quota_score >= thr_score
    winlen = w
    aaref = [aa for aa in aas if aa.id == refid][0]
    assert aaref == aas[0]
    res = []
    nfails = 0
    for aa in aas:
        # if aa.id == seqrefid:
        #     continue
        maxshiftl = maxshift + max(0, len(aa) - len(aaref))
        maxshiftr = maxshift + max(0, len(aaref) - len(aa))
        score, start = shift_and_find_best_word(aa, words, i, winlen, submat(scoring), maxshiftl, maxshiftr)
        if start is None:
            raise ValueError('zero length sequence - that should not happen')
        if score < thr_quota_score:
            nfails += 1
            # short way out 50 % fail
            if nfails / len(aas) >= 0.5:
                return
        res.append((score, start, str(aa)[start:start+winlen], aa, score >= thr_score))
    if method != 'simple':
        for score, start, word, aa_, above_thres in sorted(res, reverse=True):
            if above_thres and word not in words:
                words.add(word)
                return 'continue'
    if nfails / len(aas) > 1 - thr_quota:
        return
    flukes = []
    for _, start, word, aa, above_thres in res:
        # score = median(corrscore(word, w) for w in words)  # this is done in anchor._calculate_scores()
        fluke=Fluke(seqid=aa.id, score=None, start=start, offset=aa.meta.get('offset'), stop=start+winlen, word=word, above_thres=above_thres)
        flukes.append(fluke)
    assert flukes[0].seqid == refid
    anchor = Anchor(flukes, refid=refid)
    anchor._calculate_scores()
    return anchor


def find_anchors_winlen(aas, options, indexrange=None, anchors=None, pbar=False):
    """
    Find multiple anchors in aa sequences for a specific word length
    """
    for i, aa in enumerate(aas):
        if aa.id == options.refid:
            break
    else:
        raise ValueError(f'No sequence with ref id {options.refid}')
    aas.insert(0, aas.pop(i))
    aaref = aas[0]
    if indexrange is None:
        indexrange = list(range(len(aaref)-options.w))
    if anchors is None:
        anchors = []
    if pbar:
        desc = '{:3d} anchors found, check candidates'
        pbar = tqdm(desc=desc.format(0), total=len(indexrange))
    for i in indexrange:
        words = {str(aaref)[i:i+options.w]}
        anchor = 'continue'
        while anchor == 'continue':
            anchor = anchor_for_words(i, words, aas, **options)
        if anchor is not None:
            anchors.append(anchor)
        if pbar:
            if pbar.update():
                pbar.set_description(desc.format(len(anchors)))
    anchors = sorted(anchors, key=lambda a: a.ref.start)
    return AnchorList(anchors, options=options)


def find_my_anchors(seqs, options, remove=True, continue_with=None, **kw):
    """
    Find and return anchors in CDS region of nucleotide sequences
    """
    if continue_with is None:
        all_offset = all('offset' in seq.meta for seq in seqs)
        all_cds = all('features' in seq.meta and seq.fts.get('cds') for seq in seqs)
        if all_offset:
            log.info('Found offsets in sequence file, translate full sequence')
            aas = seqs.translate(complete=True)
            for aa in aas:
                if '*' in str(aa):
                    log.error(f'Stop codon in aa sequence {aa.id}')
        elif all_cds:
            log.info('Found CDS features in sequence file, translate CDS')
            aas = seqs['cds'].translate()
            for aa in aas:
                aa.meta.offset = aa.fts.get('cds').loc.start
        else:
            log.warning('Did not found CDS annotation or offset for at least one sequence, do not translate')
            aas = seqs
        log.debug('result of translation are {}'.format(aas.tostr(h=0)))
        log.info(f'find anchors for word length {options.w}')
        anchors = find_anchors_winlen(aas, options, **kw)
        log.info(f'found {len(anchors)} anchors for word length {options.w}')
        anchors = anchors.merge_neighbor_anchors()
        log.info(f'merged into {len(anchors)} anchors')
    else:
        anchors = continue_with
    if remove:
        removed_anchors = anchors.remove_contradicting_anchors()
        removed_anchors.data = sorted(removed_anchors, key=lambda a: a.ref.start)
        log.info(f'After removal of contradicting anchors {len(anchors)} anchors left')
    else:
        removed_anchors = None
    anchors.data = sorted(anchors, key=lambda a: a.ref.start)
    assert all([f[0].seqid == options.refid for f in anchors])
    return anchors, removed_anchors


def _split_cutout_pos(pos, mode, seqs, anchors):
    """
    Split a string, .i.e. "A10>+5", into its three parts A, B, C
    """

    if (sign:='+') in pos or (sign:='-') in pos:
        pos, C = pos.split(sign)
        C = int(sign + C)
    else:
        C = 0
    if (B := pos[-1]) in '<>^':
        pos = pos[:-1]
        if pos in ('start', 'end'):
            raise ValueError('<>^ alignment characters not allowed for start, end')
    else:
        B = '^'
    if pos == 'start' and C < 0:
        raise ValueError('C<0 not allowed for start')
    if pos == 'end' and C > 0:
        raise ValueError('C>0 not allowed for end')
    if pos in ('atg', '*') and mode != 'seq':
        raise ValueError(f'{pos} only allowed in mode seq')
    A = pos
    if is_real_anchor := (A not in ('start', 'end', 'atg', '*')):
        A = anchors[int(A.removeprefix('a'))].d
        if len(ids := set(A.keys()) - set(seqs.keys())) > 0:
            warn(f'Some anchor ids {ids} not present in sequences')
    return A, B, C, is_real_anchor


def _transform_cutout_index(A, B, C, id_, seq, mode):
    """
    Transform a fluke and position given by A,B,C to and an integer index
    """
    if A == 'start':
        i1 = i2 = 0
    elif A == 'end':
        i1 = i2 = len(seq)
    elif A == 'atg':
        assert mode == 'seq'
        i1 = seq.fts.get('cds').loc.start
        i2 = i1 + 3
    elif A == '*':
        assert mode == 'seq'
        i2 = seq.fts.get('cds').loc.stop
        i1 = i2 - 3
    else:
        i1 = _apply_mode(A[id_].start, A[id_].offset, mode=mode)
        i2 = _apply_mode(A[id_].stop, A[id_].offset, mode=mode)
    i = i1 if B == '<' else i2 if B == '>' else (i1+i2)//2
    i += C
    return i


def cutout(seqs, anchors, pos1, pos2, mode='seq'):
    """
    Cutout subsequences from pos1 to pos2 (i.e. between two anchors)
    """
    seqs = seqs.d
    la, lb, lc, l_is_real_anchor = _split_cutout_pos(pos1.strip().lower(), mode, seqs, anchors)
    ra, rb, rc, r_is_real_anchor = _split_cutout_pos(pos2.strip().lower(), mode, seqs, anchors)
    seqs2 = BioBasket()
    for id_ in seqs:
        if (l_is_real_anchor and id_ not in la) or (r_is_real_anchor and id_ not in ra):
            # no fluke for this sequence
            continue
        i = _transform_cutout_index(la, lb, lc, id_, seqs[id_], mode)
        j = _transform_cutout_index(ra, rb, rc, id_, seqs[id_], mode)
        seq2 = seqs[id_][i:j]
        if mode == 'seq':
            seq2.meta.offset = i
            seq2.meta.pop('features', None)
        seqs2.append(seq2)
    return seqs2


def combine(lot_of_anchors):
    """
    Combine lists of anchors into a single AnchorList

    Deal with possibly different offset values.
    """
    anchors = set()
    offsets = {f.seqid: f.offset for anchor in lot_of_anchors[0] for f in anchor}
    for nans in lot_of_anchors:
        if len(set(nans) & anchors) > 0:
            ids = ', '.join(a.id for a in nans & anchors)
            raise ValueError(f'Anchors {ids} exits in multiple files')
        for anchor in nans:
            for f in anchor:
                doff = f.offset - offsets[f.seqid]
                f.offset = offsets[f.seqid]
                f.start += doff // 3
                f.stop += doff // 3
        anchors |= set(nans)
    return AnchorList(sorted(anchors, key=lambda a: a.ref.start))
