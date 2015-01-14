#!/usr/bin/env python

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

# Filters for unique, highest scoring subread query/target pairs from an m4
# file. Helps get rid of chimeras, at the cost of some yield.

import sys
from collections import namedtuple

M4Record = namedtuple('M4Record', ('qname tname score pctsimilarity qstrand '
                                   'qstart qend qseqlength tstrand tstart '
                                   'tend tseqlength mapqv'))


class Count(object):
    """Tracks record count for original and filtered"""
    def __init__(self):
        self.orig = 0
        self.filt = 0

    def __repr__(self):
        return "Record count: original=%i, filtered=%i\n" % \
            (self.orig, self.filt)


def printUniq(qgroup, count):
    top = dict()
    for q in qgroup:
        m = M4Record._make(q.split())
        k = "%s%s" % (m.qname, m.tname)
        if k in top:
            n = M4Record._make(top[k].split())
            if int(m.score) < int(n.score):
                top[k] = q
        else:
            top[k] = q

    for r in top.values():
        count.filt += 1
        print r,

    qgroup[:] = []


def main():
    m4file = sys.argv[1]
    m4Hndl = open(m4file)
    qgroup = []
    curr = ''
    count = Count()
    for rec in m4Hndl:
        count.orig += 1
        m = M4Record._make(rec.split())
        if curr != m.qname:
            printUniq(qgroup, count)
            qgroup.append(rec)
            curr = m.qname
        else:
            qgroup.append(rec)

    printUniq(qgroup, count)
    sys.stderr.write(str(count))

if __name__ == '__main__':
    sys.exit(main())
