#!/usr/bin/env python
"""Super-simple converter from blasr m4 alignments to pbdagcon 'pre'
alignments. For use in the pre-assembler dagcon workflow.
"""
#################################################################################$$
# Copyright (c) 2011-2015, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted (subject to the limitations in the
# disclaimer below) provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above
#  copyright notice, this list of conditions and the following
#  disclaimer in the documentation and/or other materials provided
#  with the distribution.
#
#  * Neither the name of Pacific Biosciences nor the names of its
#  contributors may be used to endorse or promote products derived
#  from this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
# GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
# BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
#################################################################################$$


import sys
import heapq
import string # pylint: disable=W0402
from itertools import ifilter
from collections import namedtuple, defaultdict
import numpy as np
from pbcore.io.FastaIO import FastaReader

# qname tname score pctsimilarity qstrand qstart qend qseqlength tstrand tstart
# ... tend tseqlength mapqv
#
# store only fields we need
__m4fields__ = [0, 1, 2, 5, 6, 8, 9, 10, 11]
M4RECORD = namedtuple(
    'M4RECORD', 'qname tname score qstart qend tstrand tstart tend tseqlength')

__tuplfy__ = M4RECORD._make # pylint: disable=W0212

# dna compliment
__rc__ = string.maketrans('actgACTG', 'tgacTGAC')


def parse_m4(rec):
    """Parse in the m4 file, returning a list of records"""
    return [y for (x, y) in enumerate(rec.split()) if x in __m4fields__]


def rating(rec):
    """Rates the alignment for by length and score (revisit at some point)"""
    score = -int(rec.score)
    alen = int(rec.tend) - int(rec.tstart)
    return score + alen


def schwartzian(rec):
    """Provides a schwartzian transform for the given record, used for 
    sorting
    """
    flds = rec.split()
    return (flds[1], float(flds[2]), rec)


def sort_targ_score(recs):
    """Sorts the list in place by target id (string), then score (float)"""
    recs[:] = [schwartzian(x) for x in recs]
    recs.sort()
    recs[:] = [rec for (target, score, rec) in recs] # pylint: disable=W0612


def rescore(recs):
    """Rescore alignments using coverage based statistics"""
    prev = ""
    cov = np.zeros(1)
    for idx, rec in enumerate(recs):
        fields = rec.split()
        rec = __tuplfy__(fields)
        if rec.tname != prev:
            prev = rec.tname
            cov = np.zeros(int(rec.tseqlength), dtype=np.float16)

        if rec.tstrand:
            start = int(rec.tseqlength) - int(rec.tend)
            end = int(rec.tseqlength) - int(rec.tstart)
        else:
            start = int(rec.tstart)
            end = int(rec.tend)

        cov[start:end] += 1
        score = np.sum(1/cov[start:end])
        fields[2] = str(-score)
        recs[idx] = " ".join(fields)


def bestn_true(recstr, myq):
    """Checks if the record falls inside bestn (used when blasr is chunked)"""
    rec = __tuplfy__(recstr.split())
    rate = rating(rec)
    return rate in myq[rec.qname[32:]].top


class AlnLimiter(object): # pylint: disable=R0903
    """Functor that returns alignments until some count is reached. Alignments
    should be sorted. 
    """
    def __init__(self, limit=76):
        self.count = 0
        self.target = ''
        self.limit = limit

    def __call__(self, rec):
        target = rec.split()[1]
        if target != self.target:
            self.count = 0
            self.target = target
        self.count += 1
        return self.count < self.limit


class TopAlignments(object): # pylint: disable=R0903
    """Tracks the top alignments for a given query, used for bestn calc"""
    bestn = 10

    def __init__(self):
        self.top = [0] * TopAlignments.bestn

    def __call__(self):
        return  # noop

    def add(self, aln):
        """Adds an alignment to a bounded list, kicking out another if 
        necessary
        """
        heapq.heappushpop(self.top, aln)


def main(): # pylint: disable=R0914
    """Drives the program"""
    mym4 = sys.argv[1]
    allm4 = sys.argv[2]
    reads = sys.argv[3]
    TopAlignments.bestn = int(sys.argv[4])

    # tracks bestn
    my_queries = defaultdict(TopAlignments)
    my_m4recs = []

    # load my m4 chunk
    m4h = open(mym4)
    rec_add = my_m4recs.append
    for line in m4h:
        flds = parse_m4(line)
        rec = __tuplfy__(flds)
        rate = rating(rec)
        my_queries[rec.qname[32:]].add(rate)
        rec_add(' '.join(flds))

    m4h.close()

    # if we're chunked locate relevant alignments
    if mym4 != allm4:
        # assuming fofn here
        m4files = [x.rstrip() for x in open(allm4) if x.rstrip() != mym4]
        for m4f in m4files:
            m4h = open(m4f)
            for recstr in m4h:
                rec = __tuplfy__(parse_m4(recstr))
                if rec.qname[32:] in my_queries:
                    rate = rating(rec)
                    my_queries[rec.qname[32:]].add(rate)
            m4h.close()

        # remove alignments that fall outside of bestn
        my_m4recs[:] = [x for x in my_m4recs if bestn_true(x, my_queries)]

    # sort by target name/score
    sort_targ_score(my_m4recs)

    # rescore based on coverage
    rescore(my_m4recs)

    # sort one more time be new score
    sort_targ_score(my_m4recs)

    # take a max number of alignments for each target
    limiter = AlnLimiter()
    my_m4recs[:] = [x for x in ifilter(limiter, my_m4recs)]

    # load only related sequences
    seqs = {}
    frea = FastaReader(reads)
    for fent in frea:
        if fent.name[32:] in my_queries:
            seqs[fent.name] = fent.sequence

    # may or may not help
    del my_queries

    # generate pre-alignments
    for recstr in my_m4recs:
        rec = __tuplfy__(recstr.split())

        # Bug 24538, rare case missing self hit
        if rec.tname not in seqs:
            msg = "Warning: skipping query %s target %s\n"
            sys.stderr.write(msg % (rec.qname, rec.tname))
            continue

        qst = int(rec.qstart)
        qnd = int(rec.qend)
        qseq = seqs[rec.qname][qst:qnd]
        strand = '-' if rec.tstrand == '1' else '+'
        tst = int(rec.tstart)
        tnd = int(rec.tend)
        if strand == '+':
            tseq = seqs[rec.tname][tst:tnd]
        else:
            tseq = seqs[rec.tname].translate(__rc__)[::-1][tst:tnd]

        print ' '.join([rec.qname, rec.tname, strand,
                       rec.tseqlength, str(tst), str(tnd), qseq, tseq])

if __name__ == '__main__':
    sys.exit(main())
