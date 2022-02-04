/* -------------------------------------------------------------

This file is a component of SDPA
Copyright (C) 2004-2012 SDPA Project

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

------------------------------------------------------------- */

/*-----------------------------------------
  sdpa_chordal.h
-----------------------------------------*/

#ifndef __sdpa_chordal_h__
#define __sdpa_chordal_h__

#include "sdpa_dataset.h"
#include <dmumps_c.h>
  
#define SELECT_MUMPS_BEST  7 // MUMPS selects automatically when 7
#define SELECT_DENSE      -1 // This value must be minus

namespace sdpa {

class Chordal {
public:

  // condition of sparse computation
  // m_threshold < mDim, 
  // b_threshold < nBlock, 
  // aggregate_threshold >= aggrigated sparsity ratio
  // extend_threshold    >= extended sparsity ratio
  int m_threshold;
  int b_threshold;
  double aggregate_threshold;
  double extend_threshold;

  int   best;
/* indicates the used ordering method */
  /* -1: dense  computation */
  /*  7: sparse computation by MUMPS */

  SparseMatrix* sparse_bMat_ptr;
  DMUMPS_STRUC_C mumps_id;
  bool mumps_usage;
  

  Chordal(void);
  ~Chordal();
  void initialize(SparseMatrix* sparse_bMat_ptr);
  void terminate();

  // merge array1 to array2
  void mergeArray(int na1, int* array1, int na2, int* array2);
  void catArray(int na1, int* array1, int na2, int* array2);
  void slimArray(int i, int length, int* array, int& slimedLength);
 
  void makeGraph(InputData& inputData, int m);
  
  void ordering_bMat(int m, int nBlock,
                     InputData& inputData, FILE* Display,
                     FILE* fpOut);
  double analysisAndcountLowerNonZero(int m);
  bool factorizeSchur(int m, int* diagonalIndex,
		      FILE* Display, FILE* fpOut);
  bool solveSchur(Vector& rhs);

};

} // end of namespace 'sdpa'

#endif // __sdpa_chordal_h__
