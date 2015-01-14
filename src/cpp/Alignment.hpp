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


#ifndef __GCON_ALIGNMENT_HPP__
#define __GCON_ALIGNMENT_HPP__

/// 
/// Super-simple alignment representation.  Represents an alignment between two
/// PacBio reads, one of which we're trying to correct.  The read to correct
/// may be either the target or the query, depending on how the alignment was
/// done.
///
namespace dagcon {
class Alignment {
public:
    typedef void (*ParseFunc)(std::istream&, Alignment* aln);

    // May correct the target or the query, default is target
    static bool groupByTarget;

    // length of the sequence we are trying to correct
    uint32_t tlen;
    
    // conforming offsets are 1-based
    uint32_t start;

    uint32_t end;

    // ID of the read we're trying to correct (target)
    std::string id;

    // ID of the supporting read (query)
    std::string sid; 

    char strand;

    // query and target strings must be equal length
    std::string qstr;
    std::string tstr;

    Alignment();

    static ParseFunc parse;
};
}

std::istream& operator>>(std::istream& instrm, dagcon::Alignment& data);

void parseM5(std::istream& stream, dagcon::Alignment* aln);

void parsePre(std::istream& stream, dagcon::Alignment* aln);

/// Simplifies the alignment by normalizing gaps.  Converts mismatches into
/// indels ... 
///      query: CAC        query:  C-AC
///             | |  --->          |  |
///     target: CGC       target:  CG-C
///
/// Shifts equivalent gaps to the right in the reference ...
///      query: CAACAT        query: CAACAT
///             | | ||  --->         |||  |
///     target: C-A-AT       target: CAA--T
///
/// Shifts equivalent gaps to the right in the read ...
///      query: -C--CGT       query: CCG--T
///              |  | |  --->        |||  |
///     target: CCGAC-T      target: CCGACT  
dagcon::Alignment normalizeGaps(dagcon::Alignment& aln);

void trimAln(dagcon::Alignment& aln, int trimLen=50);

std::string revComp(std::string& seq); 

#endif // __GCON_ALIGNMENT_HPP__
