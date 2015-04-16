/**
 * Copyright (c) 2014-2015 Jochen Singer, David Seifert
 *
 * This file is part of PIDalign
 *
 * PIDalign is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or any later version.
 *
 * PIDalign is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * PIDalign. If not, see <http://www.gnu.org/licenses/>.
 */

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

