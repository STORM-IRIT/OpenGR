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

#ifndef __sdpa_detaset_h__
#define __sdpa_detaset_h__

#include "sdpa_include.h"
#include "sdpa_struct.h"

namespace sdpa {

class Newton;

class Solutions;
class InputData;
class Residuals;
class WorkVariables;

class ComputeTime;
class Parameter;
class StepLength;
class DirectionParameter;
class Switch;
class RatioInitResCurrentRes;
class SolveInfo;
class Phase;
class AverageComplementarity;


class Solutions
{
public:
  int nDim;
  int mDim;

  DenseLinearSpace xMat;
  DenseLinearSpace zMat;
  Vector           yVec;

  DenseLinearSpace invCholeskyX;
  DenseLinearSpace invCholeskyZ;
  DenseLinearSpace invzMat;

  double xzMinEigenValue;

  Solutions();
  Solutions(int m, BlockStruct& bs,
	    double lambda,ComputeTime& com);
  ~Solutions();
  void initialize(int m, BlockStruct& bs,
		  double lambda,ComputeTime& com);
  void terminate();

  void initializeZero(int m, BlockStruct& bs,
		      ComputeTime& com);
  
  void copyFrom(Solutions& other);
  bool update(StepLength& alpha, Newton& newton,
	      WorkVariables& work,
	      ComputeTime& com);
  bool computeInverse(WorkVariables& work,
		      ComputeTime& com);
  void display(FILE* fpout=stdout);
};

class InputData
{
public:
  Vector b;
  SparseLinearSpace C;
  SparseLinearSpace* A;

  // nBLock : number of block
  // nConstraint[k]: number of nonzero matrix in k-th block
  // When A[i].block[k] is nonzero matrix,  for t,
  //     i             <-> constraint[k][t]
  //     A[i].block[k] <-> A[i].sp_block[blockIndex[k][t]]
  int SDP_nBlock;  int* SDP_nConstraint;
  int** SDP_constraint;  int** SDP_blockIndex;
  int SOCP_nBlock;  int* SOCP_nConstraint;
  int** SOCP_constraint;  int** SOCP_blockIndex;
  int LP_nBlock;  int* LP_nConstraint;  
  int** LP_constraint;  int** LP_blockIndex;

  InputData();
  ~InputData();
  void initialize(BlockStruct& bs);
  void terminate();
  void initialize_bVec(int m);
  void initialize_index_SDP();
  void initialize_index_SOCP();
  void initialize_index_LP();
  void initialize_index();

  //   retVec_i := A_i bullet xMat (for i)
  void multi_InnerProductToA(DenseLinearSpace& xMat,Vector& retVec);
  //   retMat := \sum_{i} A_i xVec_i
  void multi_plusToA(Vector& xVec, DenseLinearSpace& retMat);
  void display(FILE* fpout=stdout);
  void display_index(FILE* fpout=stdout);
};

class Residuals
{
public:
  Vector           primalVec;
  DenseLinearSpace dualMat;
  double           normPrimalVec;
  double           normDualMat;
  double           centerNorm;

  Residuals();
  Residuals(int m, BlockStruct& bs,
	    InputData& inputData, Solutions& currentPt);
  ~Residuals();

  void initialize(int m, BlockStruct& bs,
		  InputData& inputData, Solutions& currentPt);
  void terminate();

  void copyFrom(Residuals& other);
  
  double computeMaxNorm(Vector& primalVec);
  double computeMaxNorm(DenseLinearSpace& dualMat);

  void update(int m,
	      InputData& inputData,
	      Solutions& currentPt,
	      ComputeTime& com);
  void compute(int m, 
	       InputData& inputData, 
	       Solutions& currentPt);
  void display(FILE* fpout = stdout);

};


class WorkVariables
{
public:
  DenseLinearSpace DLS1;
  DenseLinearSpace DLS2;

  // Vector DV1;
  // Vector DV2;

  BlockVector SDP_BV1;
  BlockVector SDP_BV2;
  BlockVector SDP_BV3;
  BlockVector SDP_BV4;
  BlockVector SDP_BV5;
  BlockVector SDP_BV6;
  BlockVector SDP_BV7;
  BlockVector SDP_BV8;
  BlockVector SDP_BV9;

  BlockVector SDP2_BV1;

  WorkVariables();
  WorkVariables(int m, BlockStruct& bs);
  ~WorkVariables();

  void initialize(int m, BlockStruct& bs);
  void terminate();

};

} // end of namespace 'sdpa'

#endif // __sdpa_dataset_h__
