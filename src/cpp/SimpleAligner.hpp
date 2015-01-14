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


#ifndef __GCON_SIMPLE_ALIGNER__
#define __GCON_SIMPLE_ALIGNER__
#include "Types.h"
#include "Enumerations.h"
#include "DNASequence.hpp"
#include "datastructures/alignment/Alignment.hpp"
#include "algorithms/alignment/AlignmentUtils.hpp"
#include "algorithms/alignment/SDPAlign.hpp"
#include "algorithms/alignment/GuidedAlign.hpp"
#include "format/StickAlignmentPrinter.hpp"
#include "FASTQSequence.hpp"

namespace Aligner {
struct Config {
    float indelRate;
    int indel;
    int match;
    int sdpIndel;
    int sdpIns;
    int sdpDel;
    int kmer;
    int bandSize;
};
}

class SimpleAligner {
public:
    SimpleAligner();
    void align(dagcon::Alignment& aln); 
    void operator() (dagcon::Alignment& aln);
private:
    Aligner::Config config_;
    TupleMetrics tupleMetrics_;
    DistanceMatrixScoreFunction<DNASequence, FASTQSequence> distScoreFn_;
};

#endif // __GCON_SIMPLE_ALIGNER__
