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

#ifndef __sdpa_newton_h__
#define __sdpa_newton_h__

#include "sdpa_chordal.h"

#define SparseCholesky 1

#include <pthread.h> // for multiple threading

#ifdef GOTO_BLAS
extern "C" {
  void goto_set_num_threads(int);
};
#endif

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



class Newton
{
public:
  enum bMat_Sp_De {SPARSE, DENSE};
  bMat_Sp_De bMat_type;

  SparseMatrix sparse_bMat;
  DenseMatrix bMat; // the coefficent of Schur complement
  Vector      gVec; // the right hand side of Schur complement
  DenseLinearSpace DxMat;
  Vector           DyVec;
  DenseLinearSpace DzMat;

  DenseLinearSpace r_zinvMat;
  DenseLinearSpace x_rd_zinvMat;
  
  enum FormulaType {F1,F2,F3};
  FormulaType** useFormula;
  
  // Caution: 
  // if SDPA doesn't use sparse bMat, following variables are indefinite.
  //
  // nBLock : number of block
  // nConstraint[k]: number of combination of nonzero matrices in k-th block
  // when A[k].block[i] and A[k].block[j] are nonzero matrices, 
  //     i             <-> constraint1[k][t]
  //     j             <-> constraint2[k][t]
  //     A[k].block[i] <-> A[k].sp_block[blockIndex1[k][t]]
  //     A[k].block[j] <-> A[k].sp_block[blockIndex2[k][t]]
  //     B_{ij}        <-> sparse_bMat.sp_ele[location_sparse_bMat[k][t]]
  int SDP_nBlock;  int* SDP_number;  
  int** SDP_constraint1;  int** SDP_constraint2;
  int** SDP_blockIndex1;  int** SDP_blockIndex2;
  int** SDP_location_sparse_bMat;
  int SOCP_nBlock;  int* SOCP_number;  
  int** SOCP_constraint1;  int** SOCP_constraint2;
  int** SOCP_blockIndex1;  int** SOCP_blockIndex2;
  int** SOCP_location_sparse_bMat;
  int LP_nBlock;  int* LP_number;  
  int** LP_constraint1;  int** LP_constraint2;
  int** LP_blockIndex1;  int** LP_blockIndex2;
  int** LP_location_sparse_bMat;

  // from index of aggrigate sparsity pattern to index of sparse_bMat
  // B_{ii} <-> sparse_bMat[diagonalIndex[i]]
  int* diagonalIndex;
  // B_{ij} for all i is between diagonalIndex[j] and rowStartIndex[j+1]

  Newton();
  Newton(int m, BlockStruct& bs);
  ~Newton();
  
  void initialize(int m, BlockStruct& bs);

  void terminate();

  void initialize_dense_bMat(int m);
  // 2008/03/12 kazuhide nakata
  void initialize_sparse_bMat(int m);
  // 2008/03/12 kazuhide nakata
  void initialize_bMat(int m, Chordal& chordal, InputData& inputData, 
                       FILE* Display, FILE* fpOut);

  int binarySearchIndex(int i, int j);
  void make_aggrigateIndex_SDP(InputData& inputData);
  void make_aggrigateIndex_SOCP(InputData& inputData);
  void make_aggrigateIndex_LP(InputData& inputData);
  void make_aggrigateIndex(InputData& inputData);

  void computeFormula_SDP(InputData& inputData,
			  double DenseRatio,double Kappa);

  enum WHICH_DIRECTION {PREDICTOR, CORRECTOR};
  void compute_rMat(WHICH_DIRECTION direction,
		    AverageComplementarity& mu,
		    DirectionParameter& beta,
		    Solutions& cuurentPt,
		    WorkVariables& work);

  void Make_gVec(Newton::WHICH_DIRECTION direction,
	       InputData& inputData,
	       Solutions& currentPt,
	       Residuals& currentRes,
	       AverageComplementarity& mu,
	       DirectionParameter& beta,
	       Phase& phase,
	       WorkVariables& work,
	       ComputeTime& com);

  void calF1(double& ret, DenseMatrix& G,
	     SparseMatrix& Ai);
  void calF2(double& ret, DenseMatrix& F, DenseMatrix& G,
	     DenseMatrix& invZ, SparseMatrix& Ai, bool& hasF2Gcal);
  void calF3(double& ret,
	     DenseMatrix& X, DenseMatrix& invZ,
	     SparseMatrix& Ai, SparseMatrix& Aj);

  // B_{i,j} = (X A_i Z^{-1}) \bullet A_j
  void compute_bMat_dense_SDP(InputData& inputData,
			      Solutions& currentPt,
			      WorkVariables& work,
			      ComputeTime& com);

  void compute_bMat_dense_SDP_thread(InputData& inputData,
			      Solutions& currentPt,
			      WorkVariables& work,
			      ComputeTime& com);

  // B_{i,j} = (X A_i Z^{-1}) \bullet A_j 
  void compute_bMat_sparse_SDP(InputData& inputData,
			       Solutions& currentPt,
			       WorkVariables& work,
			       ComputeTime& com);

  void compute_bMat_sparse_SDP_thread(InputData& inputData,
			       Solutions& currentPt,
			       WorkVariables& work,
			       ComputeTime& com);

  void compute_bMat_dense_SOCP(InputData& inputData,
			       Solutions& currentPt,
			       WorkVariables& work,
			       ComputeTime& com);

  void compute_bMat_sparse_SOCP(InputData& inputData,
			       Solutions& currentPt,
			       WorkVariables& work,
			       ComputeTime& com);

  void compute_bMat_dense_LP(InputData& inputData,
			     Solutions& currentPt,
			     WorkVariables& work,
			     ComputeTime& com);

  void compute_bMat_sparse_LP(InputData& inputData,
			      Solutions& currentPt,
			      WorkVariables& work,
			      ComputeTime& com);

  void Make_bMat(InputData& inputData,
		 Solutions& currentPt,
		 WorkVariables& work,
		 ComputeTime& com);

  bool compute_DyVec(Newton::WHICH_DIRECTION direction,
		     int m,
		     InputData& inputData,
		     Chordal& chordal,
		     Solutions& currentPt,
		     WorkVariables& work,
		     ComputeTime& com,
		     FILE* Display, FILE* fpOut);

  void compute_DzMat(InputData& inputData,
		     Residuals& currentRes,
		     Phase& phase,
		     ComputeTime& com);
  
  void compute_DxMat(Solutions& currentPt,
		     WorkVariables& work,
		     ComputeTime& com);
  
  
  bool Mehrotra(WHICH_DIRECTION direction,
		int m,
		InputData& inputData,
		Chordal& chordal,
		Solutions& currentPt,
		Residuals& currentRes,
		AverageComplementarity& mu,
		DirectionParameter& beta,
		Switch& reduction,
		Phase& phase,
		WorkVariables& work,
		ComputeTime& com,
		FILE* Display, FILE* fpOut);
  
  void display(FILE* fpout=stdout);
  void display_index(FILE* fpout=stdout);
  void display_sparse_bMat(FILE* fpout=stdout);

  static pthread_mutex_t job_mutex;
  static pthread_cond_t  job_cond;
  static int  Column_Number;
  static bool mutex_flag;
  static int  Calc_F1;

  static void calF1_thread(double& ret, DenseMatrix& G,
			   SparseMatrix& Aj);
  static void calF2_thread(double& ret, DenseMatrix& F, DenseMatrix& G,
			   DenseMatrix& X, SparseMatrix& Aj,
			   bool& hasF2Gcal);
  static void calF3_thread(double& ret,
			   DenseMatrix& X, DenseMatrix& invZ,
			   SparseMatrix& Ai, SparseMatrix& Aj);
  static void calF3_thread_1x1(double& ret,
			       DenseMatrix& X, DenseMatrix& invZ,
			       SparseMatrix& Ai, SparseMatrix& Aj);
  static void calF3_thread_2(double& ret,
			     DenseMatrix& X, DenseMatrix& invZ,
			     SparseMatrix& Ai, SparseMatrix& Aj);

  static  void* compute_bMat_dense_SDP_thread_func(void *arg);
  static  void* compute_bMat_sparse_SDP_thread_func(void *arg);

  int NUM_THREADS;
  int NUM_GOTOBLAS;
  void setNumThreads(FILE* Display, FILE* fpOut, int NumThreads=0);
};

 typedef struct _thread_arg {
   int Block_Number;
   int thread_num;
   int mDIM;
   int SDP_nBlock;
   int *SDP_number;
   int **SDP_constraint1;
   int **SDP_constraint2;
   int **SDP_blockIndex1;
   int **SDP_blockIndex2;
   int **SDP_location_sparse_bMat;
   DenseMatrix* bMat;
   SparseMatrix* sparse_bMat;
   Newton::FormulaType** useFormula;
   InputData* inputData;
   Solutions* currentPt;
   WorkVariables* work;
   ComputeTime* com;
 } thread_arg_t;

} // end of namespace 'sdpa'

#endif // __sdpa_newton_h__
