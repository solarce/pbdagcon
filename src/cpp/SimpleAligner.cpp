// Copyright (c) 2011-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following
// disclaimer in the documentation and/or other materials provided
// with the distribution.
//
// * Neither the name of Pacific Biosciences nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.


#include <vector>
#include <stdint.h>
#include <cstring>
#include <string>
#include <cassert>
#include <algorithm>
#include "Alignment.hpp"
#include "SimpleAligner.hpp"

SimpleAligner::SimpleAligner() {
    config_.indelRate = 0.3;
    config_.indel = 5;
    config_.match = 0;
    config_.sdpIndel = 5;
    config_.sdpIns = 5;
    config_.sdpDel = 10;
    config_.kmer = 11;
    config_.bandSize = 10;
    tupleMetrics_.Initialize(config_.kmer);
    distScoreFn_.del = config_.indel;
    distScoreFn_.ins = 4;
    distScoreFn_.InitializeScoreMatrix(SMRTDistanceMatrix);
}

void SimpleAligner::align(dagcon::Alignment& aln) {
    // This alignment type defined in blasr code base
    blasr::Alignment initialAln, refinedAln;
    FASTQSequence query;
    query.seq = (Nucleotide*)aln.qstr.c_str();
    query.length = aln.qstr.length();

    DNASequence target;
    target.seq = (Nucleotide*)aln.tstr.c_str();
    target.length = aln.tstr.length();
    SDPAlign(query, target, distScoreFn_, tupleMetrics_.tupleSize,
             config_.sdpIndel, config_.sdpIndel, config_.indelRate*2,
             initialAln, Local);

    GuidedAlign(query, target, initialAln, distScoreFn_, 
        config_.bandSize, refinedAln);

    std::string queryStr, alignStr, targetStr;

    //StickPrintAlignment(initialAln, query, target, std::cout);
    //StickPrintAlignment(refinedAln, query, target, std::cout);

    CreateAlignmentStrings(refinedAln, query.seq, target.seq, 
            targetStr, alignStr, queryStr, query.length, target.length);

    // alignment coordinates may change, update alignment object
    aln.start += refinedAln.GenomicTBegin();
    aln.end = aln.start + refinedAln.GenomicTEnd();

    if (aln.strand == '-') {
        aln.start = aln.tlen - aln.end;
        aln.qstr = revComp(queryStr);
        aln.tstr = revComp(targetStr);
    } else {
        aln.qstr = queryStr;
        aln.tstr = targetStr;
    }
    aln.start++;
}

void SimpleAligner::operator() (dagcon::Alignment& aln) {
    align(aln);
}
