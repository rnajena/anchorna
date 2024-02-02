# (C) 2023, Tom Eulenfeld, MIT license

import logging
from statistics import median

from sugar.data import submat
from tqdm import tqdm

from anchorna.util import Fluke, Anchor, AnchorList


log = logging.getLogger('anchorna')

def corrscore(seq1, seq2, gap='-', sm=submat('blosum62')):
    return sum(sm[nt1][nt2] for nt1, nt2 in zip(seq1, seq2) if nt1 != gap and nt2 != gap)


def maxes(a, key=None, default=None):
    """
    Similar to max function, but returns list of *all* maximal values
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


# def find_similar(seq, motif, starti, window, maxshift=None, maxshift_right=None):
#     seq = str(seq)
#     motif = str(motif)
#     if maxshift is None:
#         maxshift = len(seq)
#     if maxshift_right is None:
#         maxshift_right = maxshift
#     siminds = maxes([(corrscore(seq[i:i+window], motif), i)
#                       for i in range(max(0, starti-maxshift_right),
#                                      min(len(seq)-window, starti+maxshift+1))],
#                      default=(0, None))
#     if len(siminds) > 1:
#         log.warning(f'multiple max values exist: {siminds=}')
#     sim, ind = siminds[0]
#     return sim, ind


def find_words(seq, words, starti, w, sm, maxshift=None, maxshift_right=None):
    """
    Find occurance of a word similar to words in `words`
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


# def anchors_for_winlen_simple(winlen, aas, refid, maxshift, pthres, thres, exclude_thres,
#                               indexrange=None, anchors=None,
#                               show_progress=False):
#     assert thres >= exclude_thres
#     aaref = [aa for aa in aas if aa.id == refid][0]
#     if indexrange is None:
#         indexrange = list(range(len(aaref)-winlen))
#     if anchors is None:
#         anchors = []
#     for i in tqdm(list(indexrange)) if show_progress else indexrange:
#         scores = []
#         nfails = 0
#         anchor = Anchor(data=[], missing=[], refid=refid)
#         for aa in aas:
#             # if aa.id == seqrefid:
#             #     continue
#             maxshiftl = maxshift + max(0, len(aa) - len(aaref))
#             maxshiftr = maxshift + max(0, len(aaref) - len(aa))
#             score, start = find_similar(aa, aaref[i:i+winlen], i, winlen, maxshiftl, maxshiftr)
#             if start is None:
#                 log.error('zero length aa string - that should not happen')
#                 break
#             else:
#                 seqstart = aa.meta.offset+3*start
#             if score <= thres:
#                 nfails += 1
#                 if nfails / len(aas) > 1 - pthres:
#                     break
#             fluke=Fluke(id=aa.id, score=round(score, 3), start=start, seqstart=seqstart, stop=start+winlen, str=str(aa[start:start+winlen]))
#             if score <= exclude_thres:
#                 anchor.missing.append(fluke)
#             else:
#                 anchor.data.append(fluke)
#                 scores.append(score)
#         else:
#             anchors.append(anchor)
#     return AnchorList(anchors)



def _anchors_for_words(i, words, aas, w, refid, maxshift, thr_score, thr_quota_score, thr_quota,
                       scoring, method='default'):
    assert thr_quota_score >= thr_score
    winlen = w
    aaref = [aa for aa in aas if aa.id == refid][0]
    res = []
    mis = []
    nfails = 0
    for aa in aas:
        # if aa.id == seqrefid:
        #     continue
        maxshiftl = maxshift + max(0, len(aa) - len(aaref))
        maxshiftr = maxshift + max(0, len(aaref) - len(aa))
        score, start = find_words(aa, words, i, winlen, submat(scoring), maxshiftl, maxshiftr)
        if start is None:
            raise ValueError('zero length sequence - that should not happen')
        if score < thr_quota_score:
            nfails += 1
            # short way out 50 % fail
            if nfails / len(aas) >= 0.5:
                return
        if score >= thr_score:
            res.append((score, start, str(aa)[start:start+winlen], aa))
        else:
            mis.append((score, start, str(aa)[start:start+winlen], aa))
    if method != 'simple':
        for score, start, word, aa_ in sorted(res, reverse=True):
            if word not in words:
                words.add(word)
                return 'continue'
    if nfails / len(aas) > 1 - thr_quota:
        return
    flukes = []
    mflukes = []
    for score, start, word, aa in res:
        score = median(corrscore(word, w) for w in words)
        fluke=Fluke(id=aa.id, score=round(score, 3), start=start, offset=aa.meta.get('offset'), stop=start+winlen, str=word)
        flukes.append(fluke)
    for score, start, word, aa in mis:
        score = median(corrscore(word, w) for w in words)
        fluke=Fluke(id=aa.id, score=round(score, 3), start=start, offset=aa.meta.get('offset'), stop=start+winlen, str=word)
        mflukes.append(fluke)
    return Anchor(flukes, missing=mflukes, refid=refid)


def get_anchors_winlen(aas, options, indexrange=None, anchors=None, show_progress=False):
    if options.refid not in aas.ids:
        raise ValueError(f'No sequence with ref id {options.refid}')
    if indexrange is None:
        aaref = [aa for aa in aas if aa.id == options.refid][0]
        indexrange = list(range(len(aaref)-options.w))
    if anchors is None:
        anchors = []
    for i in tqdm(list(indexrange)) if show_progress else indexrange:
        words = {str(aaref)[i:i+options.w]}
        anchor = 'continue'
        while anchor == 'continue':
            anchor = _anchors_for_words(i, words, aas, **options)
        if anchor is not None:
            anchors.append(anchor)
    anchors = sorted(anchors, key=lambda a: a.ref.start)
    return AnchorList(anchors, options=options)


def get_anchors(seqs, options, remove=True):
    log.info(f'find anchors for word length {options.w}')
    anchors = get_anchors_winlen(seqs, options, show_progress=True)
    log.info(f'found {len(anchors)} anchors for word length {options.w}')
    anchors = anchors.merge_neighbor_anchors()
    log.info(f'merged into {len(anchors)} anchors')
    if remove:
        anchors = anchors.remove_contradicting_anchors()
        log.info(f'After removal of contradicting anchors {len(anchors)} anchors left')
    anchors.data = sorted(anchors, key=lambda a: a.ref.start)
    return anchors
