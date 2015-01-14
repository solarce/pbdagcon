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


#ifndef __GCON_BLASRM5_ALN_PROVIDER__
#define __GCON_BLASRM5_ALN_PROVIDER__

#include "AlnProvider.hpp"

///
/// Exceptions thrown by this class
///
namespace M5Exception {
    struct FormatError {
        std::string msg;
        FormatError(std::string m) { msg = m; }
    };

    struct SortError {};

    struct FileOpenError {};
}

///
/// Provides sets of alignments for a given target sequence from a blasr M5 
/// file.  File may be grouped by target or query.  The grouping determines
/// which set gets corrected.  Earlier, pre-assembly reads were corrected as
/// targets.  However, we can avoid the sort step if we can correct the reads 
/// as queries, since blasr groups alignments by query.
///
class BlasrM5AlnProvider : public AlnProvider {
public:
    /// Constructs a new alignment provider.  Checks the format of the file and
    /// throws an exception if it's malformed.
    /// \param fpath Path to the file containing alignments.
    BlasrM5AlnProvider(const std::string& fpath);

    /// Constructs a provider based on the given stream.  Note that no checks
    /// are actually made on the validity of the format, caveat emptor. This
    /// can be used to take a piped stream of alignments straight from blasr.
    BlasrM5AlnProvider(std::istream* stream);
    
    /// Cleans up some stuff.
    ~BlasrM5AlnProvider();

    /// Gets the set of alignments for the next target and puts them into the
    /// given vector.  Note this function will clear the contents of the vector
    /// prior to adding the next set of alignments.
    /// \param dest reference to a vector to hold the alignments.
    /// \return True if there are more targets, otherwise false.
    bool nextTarget(std::vector<dagcon::Alignment>& dest);
    
    /// Called during constructor, checks that the file is formatted correctly.
    /// Also determines if the input is grouped by query or target.
    void checkFormat();

private:
    /// Path to the input file
    const std::string fpath_;

    /// State variables 
    std::string currId_;
    dagcon::Alignment prevAln_;
    bool firstAln_;
    
    /// Represents an input stream to the alignments.
    std::ifstream fs_;
    std::istream* is_;
};

#endif //__GCON_BLASRM5_ALN_PROVIDER__
