// ==========================================================================
//                             primerIdIdentifyer
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Jochen Singer <jochen.singer@bsse.ethz.ch>
// ==========================================================================

namespace seqan
{

struct Ednafull_ {};
typedef Ednafull_ ScoreSpecEdnafull;

typedef Score<int, ScoreMatrix<Iupac, ScoreSpecEdnafull> > Ednafull;

template <typename TScore>
struct ScoringMatrixData_<TScore, Iupac, ScoreSpecEdnafull>;

template <>
struct ScoringMatrixData_<int, Iupac, ScoreSpecEdnafull> {
    enum {
        VALUE_SIZE = ValueSize<Iupac>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline int const * getData() {
        //SEQAN_CHECKPOINT;
        //Entropy =   0.1424, Expected =  -0.1074
        static int const _data[TAB_SIZE] = {
            5, 5, -4, 1, -4, 1, -4, -1, -4, 1, -4, -1, -4, -1, -4, -2,
            5, 5, -4, 1, -4, 1, -4, -1, -4, 1, -4, -1, -4, -1, -4, -2,
            -4, -4, 5, 1, -4, -4, 1, -1, -4, -4, 1, -1, -4, -4, -1, -2,
            1, 1, 1, -1, -4, -2, -2, -1, -4, -2, -2, -1, -4, -3, -3, -1,
            -4, -4, -4, -4, 5, 1, 1, -1, -4, -4, -4, -4, 1, -1, -1, -2,
            1, 1, -4, -2, 1, -1, -2, -1, -4, -2, -4, -3, -2, -1, -3, -1,
            -4, -4, 1, -2, 1, -2, -1, -1, -4, -4, -2, -3, -2, -3, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -4, -3, -3, -2, -3, -2, -2, -1,
            -4, -4, -4, -4, -4, -4, -4, -4, 5, 1, 1, -1, 1, -1, -1, -2,
            1, 1, -4, -2, -4, -2, -4, -3, 1, -1, -2, -1, -2, -1, -3, -1,
            -4, -4, 1, -2, -4, -4, -2, -3, 1, -2, -1, -1, -2, -3, -1, -1,
            -1, -1, -1, -1, -4, -3, -3, -2, -1, -1, -1, -1, -3, -2, -2, -1,
            -4, -4, -4, -4, 1, -2, -2, -3, 1, -2, -2, -3, -1, -1, -1, -1,
            -1, -1, -4, -3, -1, -1, -3, -2, -1, -1, -3, -2, -1, -1, -2, -1,
            -4, -4, -1, -3, -1, -3, -1, -2, -1, -3, -1, -2, -1, -2, -1, -1,
            -2, -2, -2, -1, -2, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1
        };
        return _data;
    }
};

}

