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

#include <algorithm>

#define PrintSparsity 1
#define PLUS_ADJUST_DIAGONAL (1.0e-10)
  // #define PLUS_ADJUST_DIAGONAL (m*1.0e-10)

  // For debug
#define FORCE_SCHUR_DENSE  0
#define FORCE_SCHUR_SPARSE 0
  
#include "sdpa_chordal.h"
#define MUMPS_JOB_INIT      -1
#define MUMPS_JOB_END       -2
#define MUMPS_JOB_ANALYSIS   1
#define MUMPS_JOB_FACTORIZE  2
#define MUMPS_JOB_SOLVE      3
#define MUMPS_USE_COMM_WORLD (-987654)

namespace sdpa {

Chordal::Chordal()
{
  mumps_usage     = false;
  sparse_bMat_ptr = NULL;
  best = 0;
}

Chordal::~Chordal()
{
  terminate();
}

void Chordal::initialize(SparseMatrix* sparse_bMat_ptr)
{
  // condition of sparse computation
  // m_threshold < mDim, 
  // b_threshold < nBlock, 
  // aggregate_threshold >= aggrigated sparsity ratio
  // extend_threshold    >= extended sparsity ratio
  m_threshold = 100;
  b_threshold = 5;
  aggregate_threshold = 0.70;
  extend_threshold = 0.80;
  
#if FORCE_SCHUR_DENSE // DENSE computation for debugging
  m_threshold = 10000000;
  b_threshold = 1000000; 
  aggregate_threshold = 0.0; 
  extend_threshold = 0.0; 
#endif
#if FORCE_SCHUR_SPARSE // SPARSE computation for debugging
  m_threshold = 0;
  b_threshold = 0; 
  aggregate_threshold = 2.0; 
  extend_threshold = 2.0; 
#endif
  // initialize by assuming Schur would be DENSE
  best = SELECT_DENSE;

  this->sparse_bMat_ptr = sparse_bMat_ptr;

  // initialize MUMPS
  mumps_id.job = MUMPS_JOB_INIT;
  mumps_id.comm_fortran = MUMPS_USE_COMM_WORLD;
  // rank 0 process participates factorizations
  mumps_id.par = 1;
  // Only symmetric positive definite matricies
  mumps_id.sym = 1;
  // No OUTPUTS
  mumps_id.icntl[1-1] = -1;
  mumps_id.icntl[2-1] = -1;
  mumps_id.icntl[3-1] = -1;
  mumps_id.icntl[4-1] =  0;

  // MUMPS selects ordering automatically
  mumps_id.icntl[7-1] =  SELECT_MUMPS_BEST;
  // for Minumum Degree Ordering
  // mumps_id.icntl[7-1] =  0;
  
  dmumps_c(&mumps_id);
  mumps_usage = true;

}

void Chordal::terminate()
{
  if (mumps_usage == true) {
    mumps_id.job = MUMPS_JOB_END;
    dmumps_c(&mumps_id);
    mumps_usage = false;
  }
  if (sparse_bMat_ptr) {
    sparse_bMat_ptr->terminate();
  }
  sparse_bMat_ptr = NULL;
}

// merge array1 to array2
void Chordal::mergeArray(int na1, int* array1, int na2, int* array2)
{

  int ptr = na1 + na2 - 1;
  int ptr1 = na1-1;
  int ptr2 = na2-1;
  int idx1, idx2;
  
  while ((ptr1 >= 0) || (ptr2 >= 0)){

    if (ptr1 >= 0){
      idx1 = array1[ptr1];
    } else {
      idx1 = -1;
    }
    if (ptr2 >= 0 ){
      idx2 = array2[ptr2];
    } else {
      idx2 = -1;
    }
    if (idx1 > idx2){
      array2[ptr] = idx1;
      ptr1--;
    } else {
      array2[ptr] = idx2;
      ptr2--;
    }
    ptr--;

  }

  // error check
  if (ptr != -1){
    rMessage("Chordal::mergeArray:: program bug");
  }
}

void Chordal::catArray(int na1, int* array1, int na2, int* array2)
{
  int ind1 = 0;
  for (int index=0; index<na1; ++index) {
    array2[na2] = array1[ind1];
    ind1++;
    na2++;
  }
}

void Chordal::slimArray(int j, int length, int* array, int& slimedLength)
{
  if (length == 0) {
    return;
  }
  sort(&array[0],&array[length]);

  // We list up only lower triangular
  int index = 0;
  while (array[index] != j) {
    index++;
  }
  array[0] = j;
  slimedLength = 0;
  index++;
  for (; index<length; ++index) {
    if (array[slimedLength] == array[index]) {
      continue;
    }
    slimedLength++;
    array[slimedLength] = array[index];
  }
  slimedLength++;
}

  // make aggrigate sparsity pattern
void Chordal::makeGraph(InputData& inputData, int m)
{

  int j,k,l;
  int SDP_nBlock  = inputData.SDP_nBlock;
  int SOCP_nBlock = inputData.SOCP_nBlock;
  int LP_nBlock   = inputData.LP_nBlock;
  int* counter;
  NewArray(counter,int,m);
  for (j=0; j<m; j++){
    counter[j] = 0;
  }

  // count maximum mumber of index
  for (l = 0; l<SDP_nBlock; l++){
    int SDP_nConstraint = inputData.SDP_nConstraint[l];
    for (k=0; k<SDP_nConstraint; k++){
      j = inputData.SDP_constraint[l][k];
      counter[j] += SDP_nConstraint;
    }
  }
  for (l = 0; l<SOCP_nBlock; l++){
    int SOCP_nConstraint = inputData.SOCP_nConstraint[l];
    for (k=0; k<SOCP_nConstraint; k++){
      j = inputData.SOCP_constraint[l][k];
      counter[j] += SOCP_nConstraint;
    }
  }
  for (l = 0; l<LP_nBlock; l++){
    int LP_nConstraint = inputData.LP_nConstraint[l];
    for (k=0; k<LP_nConstraint; k++){
      j = inputData.LP_constraint[l][k];
      counter[j] += LP_nConstraint;
    }
  }

  // allocate temporaly workspace
  int** tmp;
  NewArray(tmp,int*,m);
  for (j=0; j<m; j++){
    // Anyway B(j,j) exists even if A(j) is empty
    counter[j]++;
  }
  for (j=0; j<m; j++){
    NewArray(tmp[j],int,counter[j]);
  }
  for (j=0; j<m; j++){
    // Anyway B(j,j) exists even if A(j) is empty
    tmp[j][0] = j;
  }

  // merge index
  for (j=0; j<m; j++){
    // Slide 1: Anyway B(j,j) exists even if A(j) is empty
    counter[j] = 1;
  }
  // merge index of for SDP
  for (l = 0; l<SDP_nBlock; l++){
    for (k=0; k<inputData.SDP_nConstraint[l]; k++){
      j = inputData.SDP_constraint[l][k];
      catArray(inputData.SDP_nConstraint[l],
		 inputData.SDP_constraint[l],
                 counter[j], tmp[j]);
      counter[j] += inputData.SDP_nConstraint[l];
    }
  }
  // merge index of for SOCP
  for (l = 0; l<SOCP_nBlock; l++){
    for (k=0; k<inputData.SOCP_nConstraint[l]; k++){
      j = inputData.SOCP_constraint[l][k];
      catArray(inputData.SOCP_nConstraint[l],
	       inputData.SOCP_constraint[l],
	       counter[j], tmp[j]);
      counter[j] += inputData.SOCP_nConstraint[l];
    }
  }
  // merge index of for LP
  for (l = 0; l<LP_nBlock; l++){
    for (k=0; k<inputData.LP_nConstraint[l]; k++){
      j = inputData.LP_constraint[l][k];
      catArray(inputData.LP_nConstraint[l], inputData.LP_constraint[l],
	       counter[j], tmp[j]);
      counter[j] += inputData.LP_nConstraint[l];
    }
  }

  for (j=0; j<m; j++){
    #if 0
    printf("BeforeSlimArray[%d] = ",j);
    for (int index=0; index<counter[j]; ++index) {
      printf(" %d", tmp[j][index]);
    }
    printf("\n");
    #endif

    int tmp2 = 0;
    slimArray(j,counter[j],tmp[j],tmp2);
    counter[j] = tmp2;

    #if 0
    printf("slimArray[%d] = ",j);
    for (int index=0; index<counter[j]; ++index) {
      printf(" %d", tmp[j][index]);
    }
    printf("\n");
    #endif
  }
  int nz = 0;
  for (j=0; j<m; j++){
    nz += counter[j];
  }

  sparse_bMat_ptr -> initialize(m,m,SparseMatrix::SPARSE,
				nz,SparseMatrix::DSarrays);
  sparse_bMat_ptr -> NonZeroCount = nz;
  int indexNZ = 0;
  for (j=0; j<m; ++j) {
    for (int index_i=0; index_i<counter[j]; ++index_i) {
      // Note that MUMPS is written in FORTRAN
      // So, we need to slide all indices by +1
      sparse_bMat_ptr ->    row_index[indexNZ] = tmp[j][index_i]+1;
      sparse_bMat_ptr -> column_index[indexNZ] = j+1;
      sparse_bMat_ptr ->       sp_ele[indexNZ] = 0.0;
      indexNZ++;
    }
  }

  DeleteArray(counter);
  for (j=0; j<m; j++){
    DeleteArray(tmp[j]);
  }
  DeleteArray(tmp);
}

double Chordal::analysisAndcountLowerNonZero(int m)
{
  mumps_id.job = MUMPS_JOB_ANALYSIS;
  mumps_id.n   = m;
  mumps_id.nz  = sparse_bMat_ptr->NonZeroCount;
  mumps_id.irn = sparse_bMat_ptr->row_index;
  mumps_id.jcn = sparse_bMat_ptr->column_index;
  mumps_id.a   = sparse_bMat_ptr->sp_ele;
  
  // sparse_bMat_ptr->display();
  // rMessage("m = " << m);
  // rMessage("NonZeroCount = " << sparse_bMat_ptr->NonZeroCount);
  
  // No OUTPUTS for analysis
  mumps_id.icntl[1-1] = -1;
  mumps_id.icntl[2-1] = -1;
  mumps_id.icntl[3-1] = -1;
  mumps_id.icntl[4-1] =  0;
  // strcpy(mumps_id.write_problem,"write_problem");
  dmumps_c(&mumps_id);
  double lower_nonzeros = (double)mumps_id.infog[20-1];
  // if lower_nonzeros  is greater than 1.0e+6,
  // the value infog[20-1] is lower_nonzeros*(-1)/(1.0e+6).
  // we need to adjust the value.
  if (lower_nonzeros < 0) {
    lower_nonzeros *= (-1.0e+6);
  }
  #if 0
  rMessage("lower_nonzeros = " << lower_nonzeros);
  rMessage("Schur density = " << lower_nonzeros/((m+1)*m/2)*100 << "%" );
  #endif

  if (mumps_id.infog[1-1] != 0) {
    rError("MUMPS ERROR " << mumps_id.infog[1-1]);
  }
  return lower_nonzeros;
}

void Chordal::ordering_bMat(int m, int nBlock,
                            InputData& inputData,
			    FILE* Display, FILE* fpOut)
{
  best = SELECT_MUMPS_BEST;
  #if 0
  if ((m <= m_threshold)||(nBlock <= b_threshold)) {
    best = SELECT_DENSE;
    return;
  }
  #else
  if (m <= m_threshold) {
    best = SELECT_DENSE;
    return;
  }
  #endif

  #if 1
  for (int b=0; b<inputData.SDP_nBlock; b++){
    if (inputData.SDP_nConstraint[b] > m * sqrt(aggregate_threshold)){
      best = SELECT_DENSE;
      return;
    }      
  }
  for (int b=0; b<inputData.SOCP_nBlock; b++){
    if (inputData.SOCP_nConstraint[b] > m * sqrt(aggregate_threshold)){
      best = SELECT_DENSE;
      return;
    }      
  }
  for (int b=0; b<inputData.LP_nBlock; b++){
    if (inputData.LP_nConstraint[b] > m * sqrt(aggregate_threshold)){
      best = SELECT_DENSE;
      return;
    }     
  }
  #endif

  makeGraph(inputData,m);
  // Here, we initialize sparse_bMat

  int NonZeroAggregate = sparse_bMat_ptr->NonZeroCount*2-m;
  if (NonZeroAggregate  > aggregate_threshold * (double)m * (double) m) {
    best = SELECT_DENSE;
    return;
  }

  double lowerExtended   = analysisAndcountLowerNonZero(m);
  double NonZeroExtended = lowerExtended*2 - m;
  double overM2 = 1.0/((double)m*(double)m)*100.0;
#if PrintSparsity
  /* print sparsity information */
  if (Display) {

    #if 0
    fprintf(Display,"dense matrix               :\t\t\t%14d elements\n", m*m);
    fprintf(Display,"aggregate sparsity pattern :\t\t\t%14d elements\n",
	    NonZeroAggregate);
    fprintf(Display,"extended  sparsity pattern :\t\t\t%14d elements\n",
	    (int)NonZeroExtended);
    fprintf(Display,"Schur density = %.8lf%%\n",
	    (double)NonZeroExtended*overM2);
    fprintf(Display,"Fill in       = %e%%\n",
	    (double)(NonZeroExtended-NonZeroAggregate)*overM2);
    fprintf(Display, "Estimated FLOPs for elimation process = %e\n",
	    mumps_id.rinfog[1-1]);
    fprintf(Display,
	    "Maximum Processor  Memory Requirement = %d MB = %.2lf GB\n",
	    mumps_id.infog[16-1],(double)mumps_id.infog[16-1]/1024);
    fprintf(Display,
	    "Total   Processors Memory Requirement = %d MB = %.2lf GB\n",
	    mumps_id.infog[17-1],(double)mumps_id.infog[17-1]/1024);
    #else
    fprintf(Display, "Full Schur Elements %ld, %.2e\n",
	    (long int)((double)m*m),(double)m*m);
    fprintf(Display, "Agg %d (%.2e%%)->Ext %d (%.2e%%)"
	    " [Fill %d (%.2e%%)]\n",
	    NonZeroAggregate,
	    (double)NonZeroAggregate*overM2,
	    (int)NonZeroExtended,
	    (double)NonZeroExtended*overM2,
	    (int)(NonZeroExtended-NonZeroAggregate),
	    (double)(NonZeroExtended-NonZeroAggregate)*overM2);
    fprintf(Display, "Est FLOPs Elim = %.2e:",
	    mumps_id.rinfog[1-1]);
    fprintf(Display,
	    "MaxMem = %dMB = %.2lfGB:",
	    mumps_id.infog[16-1],(double)mumps_id.infog[16-1]/1024);
    fprintf(Display,
	    "TotMem = %dMB = %.2lfGB\n",
	    mumps_id.infog[17-1],(double)mumps_id.infog[17-1]/1024);
    #endif
  }
  if (fpOut) {
    #if 0
    fprintf(fpOut,"dense matrix               :\t\t\t%14d elements\n", m*m);
    fprintf(fpOut,"aggregate sparsity pattern :\t\t\t%14d elements\n",
	    NonZeroAggregate);
    fprintf(fpOut,"extended  sparsity pattern :\t\t\t%14d elements\n",
	    (int)NonZeroExtended);
    fprintf(fpOut,"Schur density = %.8lf%%\n",
	    (double)NonZeroExtended*overM2);
    fprintf(fpOut,"Fill in       = %e%%\n",
	    (double)(NonZeroExtended-NonZeroAggregate)*overM2);
    fprintf(fpOut, "Estimated FLOPS for elimation process = %e\n",
	    mumps_id.rinfog[1-1]);
    fprintf(fpOut,
	    "Maximum Processor  Memory Requirement = %d MB = %.2lf GB\n",
	    mumps_id.infog[16-1],(double)mumps_id.infog[16-1]/1024);
    fprintf(fpOut,
	    "Total   Processors Memory Requirement = %d MB = %.2lf GB\n",
	    mumps_id.infog[17-1],(double)mumps_id.infog[17-1]/1024);
  #else
    fprintf(fpOut, "Full Schur Elements Number %ld, %.2e\n",
    	    (long int)((double)m*m),(double)m*m);
    fprintf(fpOut, "Agg %d (%.2e%%)->Ext %d (%.2e%%)"
	    " [Fill %d (%.2e%%)]\n",
	    NonZeroAggregate,
	    (double)NonZeroAggregate*overM2,
	    (int)NonZeroExtended,
	    (double)NonZeroExtended*overM2,
	    (int)(NonZeroExtended-NonZeroAggregate),
	    (double)(NonZeroExtended-NonZeroAggregate)*overM2);
    fprintf(fpOut, "Est FLOPs Elim = %.2e:",
	    mumps_id.rinfog[1-1]);
    fprintf(fpOut,
	    "MaxMem = %dMB = %.2lfGB:",
	    mumps_id.infog[16-1],(double)mumps_id.infog[16-1]/1024);
    fprintf(fpOut,
	    "TotMem = %dMB = %.2lfGB\n",
	    mumps_id.infog[17-1],(double)mumps_id.infog[17-1]/1024);
    #endif
  }
#endif


  
  if (NonZeroExtended > extend_threshold * m * m){
    best = SELECT_DENSE;
  }

  double sparse_cost =  mumps_id.rinfog[1-1] * 1.15;
  double dense_cost  =  1.0/3.0 * (double)m * (double)m * (double) m;
  double sd_ratio = 0.85;
  // The ratio of (dense/sparse)
  // estimated by BbRosenB10.dat-s
  #if 0
  rMessage("sparse_cost = " << sparse_cost
	   << " : dense_cost = " << dense_cost
	   << " : dense_cost * sd_ratio = "
	   << (dense_cost * sd_ratio));
  #endif

  #if !FORCE_SCHUR_SPARSE
  if (sparse_cost > dense_cost * sd_ratio) {
    best = SELECT_DENSE;
  }
  #endif
}

bool Chordal::factorizeSchur(int m, int* diagonalIndex,
			     FILE* Display, FILE* fpOut)
{
  // I need to adjust Schur before factorization
  // to loose Numerical Error Condition
  double adjustSize = PLUS_ADJUST_DIAGONAL;
  for (int i=0; i<m; ++i) {
    sparse_bMat_ptr->sp_ele[diagonalIndex[i]] += adjustSize;
  }
  
  mumps_id.job = MUMPS_JOB_FACTORIZE;
  mumps_id.a   = sparse_bMat_ptr->sp_ele;
  dmumps_c(&mumps_id);
  bool isSuccess = SDPA_SUCCESS;
  while (mumps_id.infog[1-1] == -9) {
    #if 0
    rMessage("mumps icntl(14) = " << mumps_id.icntl[14-1]);
    rMessage("mumps icntl(23) = " << mumps_id.icntl[23-1]);
    #endif
    if (Display) {
      fprintf(Display,"MUMPS needs more memory space. Trying ANALYSIS phase once more\n");
    }
    if (fpOut) {
      fprintf(fpOut,  "MUMPS needs more memory space. Trying ANALYSIS phase once more\n");
    }
    mumps_id.icntl[14-1] += 20; // More 20% working memory space
    analysisAndcountLowerNonZero(m);
    mumps_id.job = MUMPS_JOB_FACTORIZE;
    dmumps_c(&mumps_id);
  }
   if (mumps_id.infog[1-1] < 0) {
    isSuccess = SDPA_FAILURE;
    if (mumps_id.infog[1-1] == -10) {
      rMessage("Cholesky failed by NUMERICAL ERROR");
      rMessage("There are some possibilities.");
      rMessage("1. SDPA terminates due to inaccuracy of numerical error");
      rMessage("2. The input problem may not have (any) interior-points");
      rMessage("3. Input matrices are linearly dependent");
    }
    else {
      rMessage("Cholesky failed with Error Code "
	       << mumps_id.infog[1-1]);
    }
  }
  return isSuccess;
}

bool Chordal::solveSchur(Vector& rhs)
{
  mumps_id.job = MUMPS_JOB_SOLVE;
  mumps_id.rhs = rhs.ele;
  dmumps_c(&mumps_id);
  return SDPA_SUCCESS;
}

} // end of namespace 'sdpa'

