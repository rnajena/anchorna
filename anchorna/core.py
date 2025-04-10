# (C) 2024, Tom Eulenfeld, MIT license
"""
Find anchors with `.find_my_anchors()`

The module also provides `.combine()` to combine/select/remove anchors and
`.cutout()` to cut out subsequences.
"""

from functools import partial
from heapq import heappush, heappop
import logging
import concurrent.futures
# import multiprocessing.pool
from warnings import warn

from sugar import BioBasket
from sugar.data import submat
from tqdm import tqdm

from anchorna.util import corrscore, Anchor, AnchorList, Fluke

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
                                     min(len(seq)-w+1, starti+maxshift+1))],
                     default=(0, None))
    sim, ind = min((abs(ind - starti), sim, ind) for sim, ind in siminds)[1:]
    if len(siminds) > 1:
        log.warning(f'multiple max values exist: {siminds=}, select {sim=}, {ind=} with smallest shift')
    return sim, ind


def anchor_at_pos(i, aas, w, gseqid, search_range,
                  score_add_word, thr_quota_add_anchor, thr_score_add_anchor,
                  scoring):
    """
    Find an anchor for a specific position i in the guiding sequence gseqid

    Return anchor or None, for description of options,
    see example configuration file and ``anchorna go -h``.

    1) Add fluke at position i for gseqid to anchor, set word to gword, set words set to {word}
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
    # if thr_score_add_anchor < score_add_word:  # there might be anchors with only poor flukes
    #     raise ValueError('Please set thr_score_add_anchor >= score_add_word.')
    winlen = w
    gaa = [aa for aa in aas if aa.id == gseqid][0]
    assert gaa == aas[0]
    gword = str(gaa)[i:i+winlen]
    gscore = corrscore(gword, gword, sm=submat(scoring))
    if gscore < thr_score_add_anchor:
        return
    words = {gword}
    todo = set(aas[1:].ids)
    aas = aas.todict()
    toadd = []
    res = [(gscore, i, gseqid)]
    nfails = 0
    while len(todo) > 0:
        for id_ in todo:
            aa = aas[id_]
            maxshiftl = search_range + max(0, len(aa) - len(gaa))
            maxshiftr = search_range + max(0, len(gaa) - len(aa))
            score, j = shift_and_find_best_word(aa, gword, i, winlen, submat(scoring), maxshiftl, maxshiftr)
            if j is None:
                raise ValueError('zero length sequence - that should not happen')
            # The lowest object is popped from heap.
            # Therefore, we add -score as first element.
            # abs(j-i) is the tiebreaker
            heappush(toadd, (-score, abs(j-i), score, j, id_))
        while len(toadd) > 0:
            _, _, score, j, id_ = heappop(toadd)
            if id_ not in todo:
                continue
            if score < thr_score_add_anchor:
                nfails += 1
                if nfails / len(aas) > 1 - thr_quota_add_anchor:
                    return
            gword = str(aas[id_])[j:j+winlen]
            res.append((score, j, id_))
            todo.discard(id_)
            if gword not in words and score >= score_add_word:
                words.add(gword)
                break
        else:
            assert len(todo) == 0
    flukes = []
    for score, j, id_ in res:
        aa = aas[id_]
        fluke=Fluke(seqid=aa.id, score=None, start=j, offset=aa.meta.get('offset'),
                    stop=j+winlen, word=str(aa)[j:j+winlen], poor=score<score_add_word)
        flukes.append(fluke)
    assert flukes[0].seqid == gseqid
    anchor = Anchor(flukes, gseqid=gseqid)
    anchor._calculate_fluke_scores(scoring)
    return anchor


def _start_parallel_jobs(tasks, do_work, results, njobs=0, pbar=True, threaded=False):
    if results is None:
        results = []
    if njobs == 0:
        log.info('sequential processing')
        mymap = map(do_work, tasks)
    else:
        try:
            from os import process_cpu_count
        except ImportError:
            from os import cpu_count as process_cpu_count
        if threaded:
            # Background: Python is moving towards a free threaded version
            # When running the tutorial with python 3.13.2, the speed-up using processes was much better.
            # Need to check with a later python version if this problem persists.
            try:
                from sys import _is_gil_enabled
            except ImportError:
                gil = True
            else:
                gil = _is_gil_enabled()
            if gil:
                msg = 'GIL detected, using threads is not feasible with this version of Python'
                raise RuntimeError(msg)
        Executor = partial(concurrent.futures.ThreadPoolExecutor, thread_name_prefix='anchorna-go') if threaded else concurrent.futures.ProcessPoolExecutor
        # Pool = multiprocessing.pool.ThreadPool if threaded else multiprocessing.pool.Pool
        njobs = njobs if njobs > 0 else process_cpu_count() + njobs
        log.info(f"use {njobs} {'threads' if threaded else 'processes'} in parallel")
        executor = Executor(max_workers=njobs)
        mymap = executor.map(do_work, tasks)
        # pool = Pool(njobs)
        # mymap = pool.imap_unordered(do_work, tasks)
    if pbar:
        desc = '{:3d} anchors found, check positions'
        pbar = tqdm(desc=desc.format(0), total=len(tasks))
    for res in mymap:
        if res is not None:
            results.append(res)
        if pbar:
            if pbar.update():
                pbar.set_description(desc.format(len(results)))
    if njobs != 0:
        executor.shutdown()
        # pool.terminate()
    return results


def find_anchors_winlen(aas, w, gseqid, indexrange=None, anchors=None, njobs=0, pbar=True, threaded=False, **kw):
    """
    Find multiple anchors in aa sequences for a specific word length

    Calls `anchor_at_pos()` for each position in the guiding sequence,
    possibly in parallel.
    """
    if str(gseqid).lower() in ('none', 'null'):
        gseqid = aas[0].id
    for i, aa in enumerate(aas):
        if aa.id == gseqid:
            break
    else:
        raise ValueError(f'No guiding sequence with id {gseqid}')
    aas.insert(0, aas.pop(i))
    gaa = aas[0]
    if indexrange is None:
        indexrange = list(range(len(gaa)-w+1))
    do_work = partial(anchor_at_pos, aas=aas, w=w, gseqid=gseqid, **kw)
    anchors = _start_parallel_jobs(indexrange, do_work, anchors, njobs=njobs, pbar=pbar, threaded=threaded)
    assert all([f[0].seqid == gseqid for f in anchors])
    return AnchorList(anchors).sort()


def find_my_anchors(seqs, remove=True, aggressive_remove=True,
                    continue_with=None, no_cds=False, scoring=None,
                    **kw):
    """
    Find and return anchors in CDS region of nucleotide sequences

    This function is called by the ``anchorna go`` command.
    For a description of arguments see the example configuration file and
    the CLI help.
    The function applies three steps:

        | A Find anchors of predefined word length in all sequences,
            this is done in the `find_anchors_winlen()` function,
            unresolved kwargs are passed on,
        | B Merge overlapping anchors with `.AnchorList.merge_overlapping_anchors()` and
        | C Remove conflicting anchors with `.AnchorList.remove_contradicting_anchors()`.
    """
    if continue_with is None:
        all_offset = all('offset' in seq.meta for seq in seqs)
        if no_cds:
            aas = seqs
            if not all_offset:
                for aa in aas:
                    aa.meta.offset = 0
        else:
            cds = [seq.fts.get('cds') for seq in seqs]
            if all_offset:
                log.info('Found offsets in sequence file, translate full sequence')
                aas = seqs.translate(complete=True)
            elif all(cds):
                log.info('Found CDS features in sequence file, translate CDS')
                aas = seqs['cds'].translate(complete=True, final_stop=False, warn=True)
                for aa in aas:
                    loc = aa.fts.get('cds').loc
                    aa.meta.offset = loc.stop if loc.strand == '-' else  loc.start
            else:
                raise ValueError('Did not found CDS annotation or offset for at least one sequence')
            log.debug('Result of translation are {}'.format(aas.tostr(h=0)))
            for aa in aas:
                if '*' in str(aa):
                    log.warning(f'Stop codon in the middle of sequence {aa.id}')
        log.info('Find anchors for specified word length')
        anchors = find_anchors_winlen(aas, scoring=scoring, **kw)
        log.info(f'Found {len(anchors)} anchors')
        anchors = anchors.merge_overlapping_anchors(scoring=scoring)
        log.info(f'Merged into {len(anchors)} anchors')
        if no_cds:
            anchors.no_cds = True
        else:
            if all(cds):
                # TODO cutout command not yet working with - strand
                strands = {ft.loc.strand for ft in cds}
                if len(strands) > 1:
                    raise ValueError('CDS features are on a different strand')
                strand = str(strands.pop())
                for anchor in anchors:
                    anchor.strand = strand
    else:
        anchors = continue_with
    if remove:
        removed_anchors = anchors.remove_contradicting_anchors(aggressive=aggressive_remove)
        log.info(f'After removal of contradicting anchors {len(anchors)} anchors left')
    else:
        removed_anchors = None
    return anchors.sort(), removed_anchors


def _split_cutout_pos(pos, mode, seqs, anchors, defaultB='^'):
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
        B = defaultB
    if pos == 'start' and C < 0:
        raise ValueError('C<0 not allowed for start')
    if pos == 'end' and C > 0:
        raise ValueError('C>0 not allowed for end')
    if pos in ('atg', '*') and mode != 'nt':
        raise ValueError(f'{pos} only allowed in mode seq')
    A = pos
    if is_real_anchor := (A not in ('start', 'end', 'atg', '*')):
        A = anchors[int(A.removeprefix('a'))].todict_seqid()
        if len(ids := set(A.keys()) - set(seqs.keys())) > 0:
            warn(f'Some anchor ids {ids} not present in sequences')
    return A, B, C, is_real_anchor


def _transform_cutout_index(A, B, C, id_, seq, mode):
    """
    Transform a fluke and position given by A,B,C to an integer index
    """
    if A == 'start':
        i1 = i2 = 0
    elif A == 'end':
        i1 = i2 = len(seq)
    elif A == 'atg':
        assert mode == 'nt'
        i1 = seq.fts.get('cds').loc.start
        i2 = i1 + 3
    elif A == '*':
        assert mode == 'nt'
        i2 = seq.fts.get('cds').loc.stop
        i1 = i2 - 3
    else:
        i1, i2 = A[id_]._apply_mode(mode)
    i = i1 if B == '<' else i2 if B == '>' else (i1+i2)//2
    i += C
    return i


def cutout(seqs, anchors, pos1, pos2, mode='nt', score_use_fluke=None, gap=None, update_fts=False):
    """
    Cutout subsequences from pos1 to pos2 (i.e. between two anchors)

    For help about this command and about how to define positions, see ``anchorna cutout -h``.
    """
    seqs = seqs.d
    la, lb, lc, l_is_real_anchor = _split_cutout_pos(pos1.strip().lower(), mode, seqs, anchors, defaultB='<')
    ra, rb, rc, r_is_real_anchor = _split_cutout_pos(pos2.strip().lower(), mode, seqs, anchors, defaultB='>')
    seqs2 = BioBasket()
    for id_ in seqs:
        if (l_is_real_anchor and id_ not in la) or (r_is_real_anchor and id_ not in ra):
            raise ValueError(f'No fluke found for sequence {id_}, this should not happen, contact developers')
        if score_use_fluke is not None and (
                (l_is_real_anchor and la[id_].score < score_use_fluke) or
                (r_is_real_anchor and ra[id_].score < score_use_fluke)):
            warn(f'do not use sequence {id_}, score is below {score_use_fluke=}')
            continue
        i = _transform_cutout_index(la, lb, lc, id_, seqs[id_], mode)
        j = _transform_cutout_index(ra, rb, rc, id_, seqs[id_], mode)
        seq2 = seqs[id_].sl(gap=gap, update_fts=update_fts)[i:j]
        if mode == 'nt':
            seq2.meta.offset = i
            if not update_fts:
                seq2.meta.pop('fts', None)
        seqs2.append(seq2)
    return seqs2


def combine(lot_of_anchors, convert_nt=False):
    """
    Combine lists of anchors into a single AnchorList

    Deal with possibly different offset values.
    Is called by ``anchorna combine``.
    """
    no_cds_set = set(anchors.no_cds for anchors in lot_of_anchors)
    if len(no_cds_set) > 1:
        raise ValueError('Some anchors calculated with no_cds option, some without')
    no_cds = no_cds_set.pop()
    if convert_nt:
        if no_cds:
            raise ValueError('--convert-nt option can only be used for CDS')
        for anchors in lot_of_anchors:
            anchors._convert_nt()
        no_cds = True
    offsets = {f.seqid: f.offset for anchor in lot_of_anchors[0] for f in anchor}
    div = 1 if no_cds else 3
    anchors = set()
    for nans in lot_of_anchors:
        if len(set(nans) & anchors) > 0:
            ids = ', '.join(a.id for a in set(nans) & anchors)
            raise ValueError(f'Anchors {ids} exist in multiple files')
        for anchor in nans:
            if anchor.strand == '-':
                raise ValueError('Cannot combine anchors on - strand, use --convert-nt option')
            for f in anchor:
                doff = f.offset - offsets[f.seqid]
                f.offset = offsets[f.seqid]
                if doff % div != 0:
                    raise ValueError('Cannot combine anchors with offset difference module not equal 3, use --convert-nt option')
                f.start += doff // div
                f.stop += doff // div
        anchors |= set(nans)
    return AnchorList(anchors, no_cds=no_cds).sort()
