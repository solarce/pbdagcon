#!/bin/bash

#################################################################################$$
# Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
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


# Simple pbdagcon workflow script.  Written for the benefit of running via
# smrtpipe so I can communicate pipe errors to the task.  We're overcoming
# the limitation of smrtpipe forcing tasks to run serially, enabling a new
# level of pipelining that's extremely efficient in an imperfect world ...
# However, direct file I/O is faster by default.

tmp=${tmp-"/tmp"}

trap "rm -f $tmp/aln.$$.pre" EXIT SIGINT

echo "Generating pre-alignments"
echo "m4topre.py $mym4 $allm4 $subreads ${bestn-24} > $tmp/aln.$$.pre"

# generate pre-alignments to a tmp directory
m4topre.py $mym4 $allm4 $subreads ${bestn-24} > $tmp/aln.$$.pre || exit $?

echo "Correcting reads"
# pipe it to consensus and generate fasta
pbdagcon -c ${cov-8} -a -j ${nproc-15} $tmp/aln.$$.pre | tee ${fasta-"corrected.fa"} | \
# generate a fastq
awk '{if($0~/>/){sub(/>/,"@",$0);print;}else{l=length($0);q="";while(l--){q=q "9"}printf("%s\n+\n%s\n",$0,q)}}' > ${fastq-"corrected.fq"}


# check the status of each piped command and exit non-zero if found
for exitval in ${PIPESTATUS[*]}
do
    if [ $exitval -gt 0 ]
    then
        exit $exitval
    fi
done


exit 0;
