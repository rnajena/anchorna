# (C) 2023, Tom Eulenfeld, MIT license

from functools import partial
from heapq import heappush, heappop
import logging
import multiprocessing
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


def shift_and_find_best_word(seq, word, starti, w, sm, maxshift=None, maxshift_right=None):
    """
    Find position of the most similar word in a sequence

    :return: tuple with similarity and index
    """
    seq = str(seq)
    if maxshift is None:
        maxshift = len(seq)
    if maxshift_right is None:
        maxshift_right = maxshift
    siminds = maxes([(corrscore(seq[i:i+w], word, sm=sm), i)
                      for i in range(max(0, starti-maxshift_right),
                                     min(len(seq)-w, starti+maxshift+1))],
                     default=(0, None))
    sim, ind = min((abs(ind - starti), sim, ind) for sim, ind in siminds)[1:]
    if len(siminds) > 1:
        log.warning(f'multiple max values exist: {siminds=}, select {sim=}, {ind=} with smallest shift')
    return sim, ind


def anchor_at_pos(i, aas, w, refid, maxshift,
                  score_add_word, thr_quota_add_anchor, thr_score_add_anchor,
                  scoring):
    """
    Find an anchor for a specific position i in the reference sequence refid

    Return anchor or None, for descitption of keyowrds,
    see example configuration file.

    1) Add fluke at position i for refid to anchor, set word to refword, set words set to {word}
    2) Add all other ids to todo list
    3) Until todo list is empty
      a) find best (score, index j) with word for each sequence in todo and add (score, j, seqid) to heap
      b) pop (score, j, seqid) pair with highest score from heap, until empty
        - if seqid not in todo -> continue 3b)
        - add to anchor, remove seqid from todos
        - if score < thr_score_add_anchor, check if thr_quota_add_anchor can still be fulfilled,
          otherwise return None (no anchor found)
        - if new word not in words, add it to words and set as new word, break loop 3b)
    4) Recalculate score, create and return anchor
    """
    winlen = w
    aaref = [aa for aa in aas if aa.id == refid][0]
    assert aaref == aas[0]
    refword = str(aaref)[i:i+winlen]
    refscore = corrscore(refword, refword, sm=submat(scoring))
    if refscore < thr_score_add_anchor:
        log.warning(f'The score of word {refword} with itself is {refscore}, smaller than {thr_score_add_anchor=}')
        return
    words = {refword}
    todo = set(aas[1:].ids)
    aas = aas.d
    toadd = []
    res = [(refscore, i, refid)]
    nfails = 0
    while len(todo) > 0:
        for id_ in todo:
            aa = aas[id_]
            maxshiftl = maxshift + max(0, len(aa) - len(aaref))
            maxshiftr = maxshift + max(0, len(aaref) - len(aa))
            score, j = shift_and_find_best_word(aa, refword, i, winlen, submat(scoring), maxshiftl, maxshiftr)
            if j is None:
                raise ValueError('zero length sequence - that should not happen')
            heappush(toadd, (-score, j, id_))
        while len(toadd) > 0:
            score, j, id_ = heappop(toadd)
            score = -score
            if id_ not in todo:
                continue
            if score < thr_score_add_anchor:
                nfails += 1
                if nfails / len(aas) > 1 - thr_quota_add_anchor:
                    return
            refword = str(aas[id_])[j:j+winlen]
            res.append((score, j, id_))
            todo.discard(id_)
            if refword not in words and score >= score_add_word:
                words.add(refword)
                break
        else:
            assert len(todo) == 0
    flukes = []
    for score, j, id_ in res:
        aa = aas[id_]
        fluke=Fluke(seqid=aa.id, score=None, start=j, offset=aa.meta.get('offset'),
                    stop=j+winlen, word=str(aa)[j:j+winlen], poor=score<score_add_word)
        flukes.append(fluke)
    assert flukes[0].seqid == refid
    anchor = Anchor(flukes, refid=refid)
    anchor._calculate_fluke_scores()
    return anchor


def _start_parallel_jobs(tasks, do_work, results, njobs=0, pbar=None):
    if results is None:
        results = []
    if njobs == 0:
        log.info('sequential processing')
        mymap = map(do_work, tasks)
    else:
        cpus = multiprocessing.cpu_count()
        njobs = min(cpus, njobs) if njobs > 0 else cpus + njobs
        log.info(f'use {njobs} cores in parallel')
        pool = multiprocessing.Pool(njobs)
        mymap = pool.imap_unordered(do_work, tasks)
    if pbar:
        desc = '{:3d} anchors found, check candidates'
        pbar = tqdm(desc=desc.format(0), total=len(tasks))
    for res in mymap:
        if res is not None:
            results.append(res)
        if pbar:
            if pbar.update():
                pbar.set_description(desc.format(len(results)))
    if njobs != 0:
        pool.terminate()
    return results


def find_anchors_winlen(aas, w, refid, indexrange=None, anchors=None, njobs=0, pbar=None, **kw):
    """
    Find multiple anchors in aa sequences for a specific word length
    """
    for i, aa in enumerate(aas):
        if aa.id == refid:
            break
    else:
        raise ValueError(f'No sequence with ref id {refid}')
    aas.insert(0, aas.pop(i))
    aaref = aas[0]
    if indexrange is None:
        indexrange = list(range(len(aaref)-w))
    do_work = partial(anchor_at_pos, aas=aas, w=w, refid=refid, **kw)
    anchors = _start_parallel_jobs(indexrange, do_work, anchors, njobs=njobs, pbar=pbar)
    assert all([f[0].seqid == refid for f in anchors])
    return AnchorList(anchors).sort()


def find_my_anchors(seqs, remove=True, aggressive_remove=True,
                    continue_with=None, **kw):
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
        log.info('find anchors for specified word length')
        anchors = find_anchors_winlen(aas, **kw)
        log.info(f'found {len(anchors)} anchors')
        anchors = anchors.merge_neighbor_anchors()
        log.info(f'merged into {len(anchors)} anchors')
    else:
        anchors = continue_with
    if remove:
        removed_anchors = anchors.remove_contradicting_anchors(aggressive=aggressive_remove)
        log.info(f'After removal of contradicting anchors {len(anchors)} anchors left')
    else:
        removed_anchors = None
    return anchors.sort(), removed_anchors


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


def cutout(seqs, anchors, pos1, pos2, mode='seq', score_use_fluke=None):
    """
    Cutout subsequences from pos1 to pos2 (i.e. between two anchors)
    """
    seqs = seqs.d
    la, lb, lc, l_is_real_anchor = _split_cutout_pos(pos1.strip().lower(), mode, seqs, anchors)
    ra, rb, rc, r_is_real_anchor = _split_cutout_pos(pos2.strip().lower(), mode, seqs, anchors)
    seqs2 = BioBasket()
    for id_ in seqs:
        if (l_is_real_anchor and id_ not in la) or (r_is_real_anchor and id_ not in ra):
            raise ValueError(f'No fluke found for sequence {id_}, this should not happen, contact devs')
        if score_use_fluke is not None and (
                (l_is_real_anchor and la[id_].score < score_use_fluke) or
                (r_is_real_anchor and ra[id_].score < score_use_fluke)):
            warn(f'do not use sequence {id_}, score is below {score_use_fluke=}')
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
            raise ValueError(f'Anchors {ids} exist in multiple files')
        for anchor in nans:
            for f in anchor:
                doff = f.offset - offsets[f.seqid]
                f.offset = offsets[f.seqid]
                f.start += doff // 3
                f.stop += doff // 3
        anchors |= set(nans)
    return AnchorList(anchors).sort()
