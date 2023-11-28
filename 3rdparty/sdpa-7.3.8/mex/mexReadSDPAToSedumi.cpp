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
/*----------------------------------------------------------
  function [At,b,c,K,blockStruct] = SDPAToSedumi(filename);
----------------------------------------------------------*/

#include <iostream>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <mex.h>

using namespace std;
#define lengthOfString 10240

#define MX_CALLOC 0
#define MX_DEBUG  0

enum BlockType {btSDP, btLP};

class LIJV
{
public:
  int l,i,j;
  int index;
  double value;

  LIJV()
  {
    i = j = l = index = 0;
    value = 0.0;
  }
  static bool compare(LIJV* a, LIJV* b)
  {
    // if  (a <  b) return true;
    // if  (a >= b) return  false;

    if      ( a[0].index < b[0].index ) { return true; }
    else if ( a[0].index > b[0].index ) { return false; }
    return false; // a == b
  }
};

void printData(vector<LIJV*>* A,int m)
{
  mexPrintf("------------------------\n");
  for (int k=0; k<m+1; ++k) {
    int length = A[k].size();
    for (int index=0; index<length; ++index) {
      LIJV* ele = A[k].at(index);
      mexPrintf("%d %d %d %d %e [%d]\n",
		k, ele->l, ele->i, ele->j, ele->value, ele->index);
    }
  }
  mexPrintf("------------------------\n");
}

void readLIJV(vector<LIJV*>* A,FILE* fpData)
{
  int i,j,k,l;
  double value;
  while (true) {
    if (fscanf(fpData,"%*[^0-9+-]%d",&k)<=0) {
      break;
    }
    if (fscanf(fpData,"%*[^0-9+-]%d",&l)<=0) {
      break;
    }
    if (fscanf(fpData,"%*[^0-9+-]%d",&i)<=0) {
      break;
    }
    if (fscanf(fpData,"%*[^0-9+-]%d",&j)<=0) {
      break;
    }
    if (fscanf(fpData,"%*[^0-9+-]%lf",&value)<=0) {
      break;
    }

    if (value == 0.0 || value == -0.0) {
      continue;
    }

    LIJV* ele = NULL;
    #if MX_CALLOC
    ele = (LIJV*)mxCalloc(1, sizeof(LIJV));
    #else
    try {
      ele = new LIJV;
    }
    catch (bad_alloc){
        mexPrintf("Memory Exhausted (bad_alloc)");
    }
    catch(...){
        mexPrintf("Fatal Error (related memory allocation)");
    }
    #endif
    ele->i = i;
    ele->j = j;
    ele->l = l;
    ele->value = value;
    A[k].push_back(ele);
  }
}

void pushNonSymmetric(vector<LIJV*>* A, int m)
{
  for (int k=0; k<m+1; ++k) {
    const int length = A[k].size();
    for (int index=0; index<length; ++index) {
      LIJV* ele1 = A[k].at(index);
      if (ele1->i != ele1->j) {
	LIJV* ele2 = NULL;
        #if MX_CALLOC
	ele2 = (LIJV*)mxCalloc(1, sizeof(LIJV));
        #else
	try {
	  ele2 = new LIJV;
	}
	catch(bad_alloc){
	  mexPrintf("Memory Exhausted (bad_alloc)");
	}
	catch(...){
	  mexPrintf("Fatal Error (related memory allocation");
	}
        #endif

	ele2->i = ele1->j;
	ele2->j = ele1->i;
	ele2->l = ele1->l;
	ele2->value = ele1->value;
	A[k].push_back(ele2);
      }
    }
  }
}

void changeIndex(vector<LIJV*>* A, int m,
		 int nBlock, BlockType* blockType,
		 int* blockStart, int* blockStruct)
{
  for (int k=0; k<m+1; ++k) {
    const int length = A[k].size();
    for (int index=0; index<length; ++index) {
      LIJV* ele = A[k].at(index);
      const int l = ele->l-1;
      const int i = ele->i-1;
      const int j = ele->j-1;
      int index2 = i;
      if (blockType[l] == btSDP) {
	index2 = i*blockStruct[l]+j;
      }
      ele->index = index2 + blockStart[l];
    }
  }
  for (int k=0; k<m+1; ++k) {
    sort(A[k].begin(), A[k].end(), LIJV::compare);
  }
}

void mexFunction(int nlhs, mxArray* plhs[],
		 int nrhs, const mxArray* prhs[])
{
  mwSize mwsize;
  mxArray *field_ptr = NULL;
  mxArray* filename_ptr = (mxArray*) prhs[0];
  // Get filename
  char* filename = NULL;
  mwsize = mxGetM(filename_ptr)*mxGetN(filename_ptr)+1;
  filename = (char*)mxCalloc(mwsize, sizeof(char));
  mxGetString(filename_ptr,filename,mwsize);

  FILE* fpData;
  if ((fpData = fopen(filename,"r")) == NULL) {
    mexPrintf("**** Cannot Open [%s] ****\n",filename);
    mexPrintf("Output [At,b,c,K,blockStruct] will be set empty\n");
    plhs[0] = mxCreateSparse(0, 0, 0,mxREAL);
    plhs[1] = mxCreateSparse(0, 0, 0,mxREAL);
    plhs[2] = mxCreateSparse(0, 0, 0,mxREAL);
    plhs[3] = mxCreateSparse(0, 0, 0,mxREAL);
    plhs[4] = mxCreateSparse(0, 0, 0,mxREAL);
    return;
  }
  // mexPrintf("Opened [%s]\n", filename);

  int m      = 0;
  int nBlock = 0;
  int* blockStruct = NULL;
  char str[lengthOfString];

  // Read Comment & m
  while (true) {
    volatile int dummy=0; dummy++;//for gcc-3.3 bug
    fgets(str,lengthOfString,fpData);
    if (str[0]=='*' || str[0]=='"') {
      mexPrintf("%s",str);
    } else {
      sscanf(str,"%d",&m);
      break;
    }
  }
  #if MX_DEBUG
  mexPrintf("m = %d\n",m);
  #endif
  fscanf(fpData,"%d",&nBlock);
  #if MX_DEBUG
  mexPrintf("nBlock = %d\n",nBlock);
  #endif
  blockStruct = (int*)mxCalloc(nBlock, sizeof(int));
  for (int l=0; l<nBlock; ++l) {
    fscanf(fpData,"%*[^0-9+-]%d",&blockStruct[l]);
  }
  for (int l=0; l<nBlock; ++l) {
    if (blockStruct[l] == 1) {
      blockStruct[l] = -1;
    }
  }

  #if MX_DEBUG
  mexPrintf("blockStruct = ");
  for (int l=0; l<nBlock; ++l) {
    mexPrintf("%d ",blockStruct[l]);
  }
  mexPrintf("\n");
  #endif

  BlockType* blockType;
  blockType = (BlockType*)mxCalloc(nBlock, sizeof(BlockType));
  int* blockStart;
  blockStart = (int*)mxCalloc(nBlock, sizeof(int));

  int nLP = 0;
  for (int l=0; l<nBlock; ++l) {
    if (blockStruct[l] < 0) {
      nLP -= blockStruct[l];
    }
  }

  int startLP  = 0;
  int startSDP = nLP;
  for (int l=0; l<nBlock; ++l) {
    if (blockStruct[l] < 0) {
      blockType [l] = btLP;
      blockStart[l] = startLP;
      startLP -= blockStruct[l];
    }
    else {
      blockType [l] = btSDP;
      blockStart[l] = startSDP;
      startSDP += blockStruct[l]*blockStruct[l];
    }
  }

  int totalLength = startSDP;

  double* b;
  b = (double*) mxCalloc(m, sizeof(double));

  for (int k=0; k<m; ++k) {
    fscanf(fpData,"%*[^0-9+-]%lf",&b[k]);
  }
  #if MX_DEBUG
  mexPrintf("b = ");
  for (int k=0; k<m; ++k) {
    mexPrintf("%e ",b[k]);
  }
  mexPrintf("\n");
  #endif

  vector<LIJV*>* A;
  A = (vector<LIJV*>*) mxCalloc(m+1, sizeof(vector<LIJV*>));

  readLIJV(A,fpData);
  #if MX_DEBUG
  printData(A,m);
  #endif
  pushNonSymmetric(A,m);
  #if MX_DEBUG
  printData(A,m);
  #endif
  changeIndex(A, m, nBlock, blockType, blockStart, blockStruct);
  #if MX_DEBUG
  printData(A,m);
  #endif

  mwSize Annz = 0;
  for (int k=1; k<m+1; ++k) {
    Annz += A[k].size();
  }

  const char *Knames[] = {
    "l","s"
  };
  
  /* Create cellarrays for the output variables */
  plhs[0] = mxCreateSparse(totalLength, m, Annz,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(m,1,mxREAL);
  plhs[2] = mxCreateSparse(totalLength,1, A[0].size(),mxREAL);
  plhs[3] = mxCreateStructMatrix(1,1,2,Knames);
  plhs[4] = mxCreateDoubleMatrix(nBlock,1,mxREAL);

  mxArray* At_ptr = plhs[0];
  mxArray* b_ptr  = plhs[1];
  mxArray* c_ptr  = plhs[2];
  mxArray* K_ptr  = plhs[3];
  mxArray* blockStruct_ptr = plhs[4];

  
  mwIndex* At_row    = mxGetIr(At_ptr);
  mwIndex* At_column = mxGetJc(At_ptr);
  double*  At_ele    = mxGetPr(At_ptr);
  At_column[0] = 0;
  for (int k=0; k<m; ++k) {
    At_column[k+1] = At_column[k] + A[k+1].size();
  }
  int cIndex = 0;
  for (int k=0; k<m; ++k) {
    int length = A[k+1].size();
    for (int index2 = 0; index2 < length; ++index2) {
      LIJV* ele = A[k+1].at(index2);
      At_row[cIndex] = ele->index;
      At_ele[cIndex] = ele->value;
      cIndex++;
    }
  }
  
  mwIndex* c_row    = mxGetIr(c_ptr);
  mwIndex* c_column = mxGetJc(c_ptr);
  double*  c_ele    = mxGetPr(c_ptr);
  c_column[0] = 0;
  c_column[1] = A[0].size();
  cIndex = 0;
  for (mwSize index2 = 0; index2 < (mwSize) A[0].size(); ++index2) {
    LIJV* ele = A[0].at(index2);
    c_row[cIndex] = ele->index;
    // Note that 'c' must be negative
    c_ele[cIndex] = -(ele->value);
    cIndex++;
  }

  double* b_ele = mxGetPr(b_ptr);
  for (int k=0; k<m; ++k) {
    b_ele[k] = b[k];
  }

  field_ptr = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
  *mxGetPr(field_ptr) = nLP;
  mxSetField(K_ptr, 0, "l", field_ptr);
  int noSDP = 0;
  for (int l=0; l<nBlock; ++l) {
    if (blockType[l] == btSDP) {
      noSDP++;
    }
  }
  field_ptr = mxCreateDoubleMatrix(noSDP,1,mxREAL);
  double* Ks = mxGetPr(field_ptr);
  int indexSDP = 0;
  for (int l=0; l<nBlock; ++l) {
    if (blockType[l] == btSDP) {
      Ks[indexSDP] = blockStruct[l];
      indexSDP++;
    }
  }
  mxSetField(K_ptr, 0, "s", field_ptr);

  double* blockStruct_ele = mxGetPr(blockStruct_ptr);
  for (int l=0; l<nBlock; ++l) {
    blockStruct_ele[l] = blockStruct[l];
  }

  for (int k=0; k<m+1; ++k) {
    int length = A[k].size();
    for (int index=0; index<length; ++index) {
      #if MX_CALLOC
      mxFree(A[k].at(index));
      #else
      delete A[k].at(index);
      #endif
    }
  }    
  
  fclose(fpData);
  mxFree(A);
  mxFree(b);
  mxFree(blockStruct);
  mxFree(blockType);
  mxFree(blockStart);
  mxFree(filename);
  return;
}

