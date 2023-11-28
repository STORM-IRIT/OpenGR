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
/*---------------------------------------------------
function [objVal,x,X,Y,INFO] = read_output(filename,mDIM,nBLOCK,bLOCKsTRUCT);
---------------------------------------------------*/

#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <mex.h>

using namespace std;
#define lengthOfString 10240

#define MX_CALLOC 0
#define MX_DEBUG  0

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{

  mwSize mwsize;
  mxArray *field_ptr = NULL;
  mxArray *cell_ptr  = NULL;
  char str[lengthOfString];
  double *tmp_ptr = NULL;

  const mxArray* filename_ptr      = prhs[0];
  const mxArray* m_ptr             = prhs[1];
  const mxArray* nBlock_ptr        = prhs[2];
  const mxArray* blockStruct_ptr   = prhs[3];
  
  char* filename = NULL;
  mwsize = mxGetM(filename_ptr)*mxGetN(filename_ptr)+1;
  filename = (char*)mxCalloc(mwsize, sizeof(char));
  mxGetString(filename_ptr,filename,mwsize);

  mwSize m = 0;
  mwSize nBlock = 0;
  double* m_ele = mxGetPr(m_ptr);
  m = (mwSize) m_ele[0];
  double* nBlock_ele = mxGetPr(nBlock_ptr);
  nBlock = (mwSize) nBlock_ele[0];
  
  double* blockStruct_ele = mxGetPr(blockStruct_ptr);
  #if MX_DEBUG
  mexPrintf("blockStruct = [");
  for (mwSize l=0; l<nBlock; ++l) {
    mexPrintf("%d ", (int) blockStruct_ele[l]);
  }
  mexPrintf("\n");
  #endif



  FILE* fpData;
  if ((fpData = fopen(filename,"r")) == NULL) {
    mexPrintf("**** Cannot Open [%s] ****\n",filename);
    mexPrintf("Output [objVal,x,X,Y,INFO] will be set empty\n");
    plhs[0] = mxCreateSparse(0, 0, 0,mxREAL);
    plhs[1] = mxCreateSparse(0, 0, 0,mxREAL);
    plhs[2] = mxCreateSparse(0, 0, 0,mxREAL);
    plhs[3] = mxCreateSparse(0, 0, 0,mxREAL);
    plhs[4] = mxCreateSparse(0, 0, 0,mxREAL);
    return;
  }

  const char *fnames[] = {
    "phasevalue",
    "iteration",
    // "dimacs", // dimacs is not available when option is not assigned
    "cpusec"
  };  

  plhs[0] = mxCreateDoubleMatrix(1,2,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(m,1,mxREAL);
  plhs[2] = mxCreateCellMatrix(nBlock,1);
  plhs[3] = mxCreateCellMatrix(nBlock,1);
  plhs[4] = mxCreateStructMatrix(1,1,3,fnames);

  double* objVal_ele = mxGetPr(plhs[0]);
  double* x_ele      = mxGetPr(plhs[1]);
  mxArray* X_ptr     = plhs[2];
  mxArray* Y_ptr     = plhs[3];
  mxArray* INFO_ptr  = plhs[4];
  
  while (true) {
    volatile int dummy=0; dummy++;//for gcc-3.3 bug
    fgets(str,lengthOfString,fpData);
    if (strstr(str,"phase.value") != NULL) {
      break;
    }
  }
  char* phasevalue = strstr(str,"phase.value") 
    + strlen("phase.value  = ");
  for (mwSize strIndex = 0; strIndex < 30; ++strIndex) {
    if ( phasevalue[strIndex] == ' '
	 || phasevalue[strIndex] == '\n') { // remove redundant space
      phasevalue[strIndex] = '\0';
      break;
    }
  }
  field_ptr = mxCreateString(phasevalue);
  mxSetField(INFO_ptr, 0, "phasevalue", field_ptr);
  #if MX_DEBUG
  mexPrintf("phasevalue = %s\n", phasevalue);
  #endif
  
  while (true) {
    volatile int dummy=0; dummy++;//for gcc-3.3 bug
    fgets(str,lengthOfString,fpData);
    if (strstr(str,"Iteration") != NULL) {
      break;
    }
  }
  char* stringIteration = strstr(str,"Iteration") 
    + strlen("Iteration = ");
  field_ptr = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
  *mxGetPr(field_ptr) = (double)atoi(stringIteration);
  mxSetField(INFO_ptr, 0, "iteration", field_ptr);
  #if MX_DEBUG
  mexPrintf("Iteration = %d\n", atoi(stringIteration));
  #endif

  while (true) {
    volatile int dummy=0; dummy++;//for gcc-3.3 bug
    fgets(str,lengthOfString,fpData);
    if (strstr(str,"objValPrimal") != NULL) {
      break;
    }
  }
  char* stringPrimal = strstr(str,"objValPrimal") 
    + strlen("objValPrimal = ");
  field_ptr = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
  objVal_ele[0] = (double)atof(stringPrimal);
  #if MX_DEBUG
  mexPrintf("primal = %.16e\n", objVal_ele[0]);
  #endif

  while (true) {
    volatile int dummy=0; dummy++;//for gcc-3.3 bug
    fgets(str,lengthOfString,fpData);
    if (strstr(str,"objValDual") != NULL) {
      break;
    }
  }
  char* stringDual = strstr(str,"objValDual") 
    + strlen("objValDual   = ");
  field_ptr = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
  objVal_ele[1] = (double)atof(stringDual);
  #if MX_DEBUG
  mexPrintf("dual = %.16e\n", objVal_ele[1]);
  #endif

  while (true) {
    volatile int dummy=0; dummy++;//for gcc-3.3 bug
    fgets(str,lengthOfString,fpData);
    if (strstr(str,"total time") != NULL) {
      break;
    }
  }
  char* stringTime = strstr(str,"total time") 
    + strlen("total time   = ");
  field_ptr = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
  *mxGetPr(field_ptr) = (double)atof(stringTime);
  mxSetField(INFO_ptr, 0, "cpusec", field_ptr);
  #if MX_DEBUG
  mexPrintf("cpusec = %lf\n", atof(stringTime));
  #endif
  
  // result vectors and matrices
  while (true) {
    volatile int dummy=0; dummy++;//for gcc-3.3 bug
    fgets(str,lengthOfString,fpData);
    if (strstr(str,"xVec = ") != NULL) {
      break;
    }
  }

  // for xVec
  double tmpd;
  for (mwSize k=0; k<m; ++k) {
    fscanf(fpData,"%*[^0-9+-]%lf",&tmpd);
    x_ele[k] = tmpd;
    #if MX_DEBUG
    mexPrintf("x_ele[%d] = %e\n",k,x_ele[k]);
    #endif
  }

  mwIndex cell_index;
  // for xMat
  cell_index = 0;
  for (mwSize l=0; l<nBlock; ++l) {
    mwSize sizeM,sizeN;
    sizeM = sizeN = 0;
    long int blk = (long int)blockStruct_ele[l];
    sizeM = sizeN = blk;
    if (blk < 0) {
      sizeM = 1;
      sizeN = -blk;
    }
    mwSize lengthMN = sizeM * sizeN;
    #if MX_DEBUG
    mexPrintf("X{%d} size  = %d, sizeM = %d, sizeN = %d, lengthMN = %d\n",
	      cell_index,blk,sizeM,sizeN,lengthMN);
    #endif
    cell_ptr = mxCreateDoubleMatrix(sizeM, sizeN, mxREAL);
    tmp_ptr = mxGetPr(cell_ptr);
    for (mwSize index = 0; index < lengthMN; ++index) {
      fscanf(fpData,"%*[^0-9+-]%lf",&tmpd);
      tmp_ptr[index] = tmpd;
      #if MX_DEBUG
      mexPrintf("X[%d,%d] = %e\n",cell_index,index,tmp_ptr[index]);
      #endif
    }
    mxSetCell(X_ptr, cell_index++, mxDuplicateArray(cell_ptr));
  }
    
  
  // for yMat
  cell_index = 0;
  for (mwSize l=0; l<nBlock; ++l) {
    mwSize sizeM,sizeN;
    sizeM = sizeN = 0;
    long int blk = (long int)blockStruct_ele[l];
    sizeM = sizeN = blk;
    if (blk < 0) {
      sizeM = 1;
      sizeN = -blk;
    }
    mwSize lengthMN = sizeM * sizeN;
    #if MX_DEBUG
    mexPrintf("Y{%d} length  = %d\n",cell_index,lengthMN);
    #endif
    
    cell_ptr = mxCreateDoubleMatrix(sizeM, sizeN, mxREAL);
    tmp_ptr = mxGetPr(cell_ptr);
    for (mwSize index = 0; index < lengthMN; ++index) {
      fscanf(fpData,"%*[^0-9+-]%lf",&tmpd);
      tmp_ptr[index] = tmpd;
      #if MX_DEBUG
      mexPrintf("Y[%d,%d] = %e\n",cell_index,index,tmp_ptr[index]);
      #endif
    }
    mxSetCell(Y_ptr, cell_index++, mxDuplicateArray(cell_ptr));
  }
    
  fclose(fpData);
  return;
}
