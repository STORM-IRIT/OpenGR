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
  function mexWriteSedumiToSDPA(filename,At,b,c,K,accuracy)

  Note:: c must be set as (-c) and K.f must be processed
  before entering this routine.
----------------------------------------------------------*/

#include <iostream>
#include <cstdio>
#include <mex.h>
using namespace std;


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  mwSize mwsize;
  mxArray *field_ptr = NULL;
  
  mxArray* filename_ptr = (mxArray*) prhs[0];
  mxArray* At_ptr       = (mxArray*) prhs[1];
  mxArray* b_ptr        = (mxArray*) prhs[2];
  mxArray* c_ptr        = (mxArray*) prhs[3];
  mxArray* K_ptr        = (mxArray*) prhs[4];
  mxArray* accuracy_ptr = (mxArray*) prhs[5];

  // Get filename
  char* filename = NULL;
  mwsize = mxGetM(filename_ptr)*mxGetN(filename_ptr)+1;
  filename = (char*)mxCalloc(mwsize, sizeof(char));
  mxGetString(filename_ptr,filename,mwsize);

  FILE* fp;
  if ((fp = fopen(filename,"w")) == NULL) {
    mexPrintf("Cannot Open %s\n",filename);
    return;
  }
  
  char* accuracy = NULL;
  mwsize = mxGetM(accuracy_ptr)*mxGetN(accuracy_ptr)+1;
  accuracy = (char*)mxCalloc(mwsize, sizeof(char));
  mxGetString(accuracy_ptr,accuracy,mwsize);
  // mexPrintf("accuracy = %s\n", accuracy);
  
  mwSize mDIM = mxGetN(At_ptr);
  fprintf(fp, PRINTF_INT_STYLE "\n", mDIM);
  
  
  /* nBLOCK */
  mwSize nBLOCK = 0;
  mwSize K_l    = 0;
  int    isK_l  = 0; // 1 (K_l > 0) or 0 (K_l == 0)
  field_ptr = mxGetField(K_ptr, 0, "l");
  if (field_ptr != NULL) {
    K_l   = (mwSize)((mxGetPr(field_ptr))[0]);
    #if 0
    mexPrintf("K_l = " PRINT_INT_STYLE "\n", K_l);
    #endif
    if (K_l > 0) {
      isK_l = 1;
      nBLOCK++;
    }
  }
  mwSize* K_s             = NULL;
  mwSize* K_sdpConeStart  = NULL;
  int     K_sdpNoCones    = 0;
  field_ptr = mxGetField(K_ptr, 0, "s");
  if (field_ptr != NULL) {
    K_sdpNoCones = (int) mxGetM(field_ptr);
    K_s = (mwSize*) mxCalloc(K_sdpNoCones, sizeof(mwSize));
    K_sdpConeStart = (mwSize*) mxCalloc(K_sdpNoCones+1, sizeof(mwSize));
    K_sdpConeStart[0] = 0;
    #if 0
    mexPrintf("K_sdpNoCones = " PRINTF_INT_STYLE "\n", K_sdpNoCones);
    #endif
    for (int l=0; l<K_sdpNoCones; ++l) {
      K_s[l] = (mwSize)((mxGetPr(field_ptr))[l]);
      #if 0
      mexPrintf("K_s[" PRINTF_INT_STYLE "] = " PRINTF_INT_STYLE "\n", l+isK_l+1, K_s[l]);
      #endif
      K_sdpConeStart[l+1] = K_sdpConeStart[l] + K_s[l]*K_s[l];
      nBLOCK++;
    }
  }

  fprintf(fp, PRINTF_INT_STYLE "\n", nBLOCK);
  if (isK_l) {
    fprintf(fp, PRINTF_INT_STYLE " ", -K_l);
  }
    
  for (int l=0; l<K_sdpNoCones-1; ++l) {
    fprintf(fp, PRINTF_INT_STYLE " ", K_s[l]);
  }
  if (K_sdpNoCones >= 1) {
    fprintf(fp, PRINTF_INT_STYLE "\n", K_s[K_sdpNoCones-1]);
  }
  else {
    fprintf(fp,"\n");
  }
    

  double* b = mxGetPr(b_ptr);
  for(mwSize k = 0; k < mDIM-1; ++k){
    fprintf(fp, accuracy, b[k]);
    fprintf(fp, " ");
  }
  fprintf(fp, accuracy, b[mDIM-1]);
  fprintf(fp, "\n");

  /* F_0 = - C */ 
  if (mxIsEmpty(c_ptr) || mxGetNzmax(c_ptr) == 0
      || (mxGetJc(c_ptr))[1] == 0) {
    mexPrintf("c = empty\n");
  }
  else {
    mwIndex* C_row     = mxGetIr(c_ptr);
    mwIndex* C_column  = mxGetJc(c_ptr);
    double*  C_ele     = mxGetPr(c_ptr);
    mwIndex  C_length  = C_column[1];
    int currentSdpCone = 0;
    for (mwSize index = 0; index<C_length; ++index) {
      mwSize C_j =  C_row[index];
      double ele =  C_ele[index];
      // mexPrintf("C_j = " PRINTF_INT_STYLE ", ele = %e\n", C_j, ele);
      if (C_j < K_l) {
	fprintf(fp, "0 1 " PRINTF_INT_STYLE " " PRINTF_INT_STYLE " ", C_j+1, C_j+1);
	fprintf(fp, accuracy, ele);
	fprintf(fp, "\n");
      }
      else {
	C_j -= K_l;
	while (K_sdpConeStart[currentSdpCone+1] <= C_j) {
	  currentSdpCone++;
	}
	mwSize index2, i,j;
	index2 = C_j - K_sdpConeStart[currentSdpCone];
	i = index2 / K_s[currentSdpCone];
	j = index2 % K_s[currentSdpCone];
	if (i <= j) {
	  // Only upper triangular is input
	  #if 0
	  mexPrintf("input k=%d, l=%d, i=%d, j=%d, v=%e\n",
		    0, isK_l + currentSdpCone+1,i+1, j+1, ele);
	  #endif
	  fprintf(fp, "0 " PRINTF_INT_STYLE " " PRINTF_INT_STYLE " " PRINTF_INT_STYLE " ",
		  (mwSize)(isK_l + currentSdpCone + 1), i+1, j+1);
	  fprintf(fp, accuracy, ele);
	  fprintf(fp, "\n");
	}
      }
    }
  }

  /* F_{k+1} = A_k */
  mwIndex* At_row     = mxGetIr(At_ptr);
  mwIndex* At_column  = mxGetJc(At_ptr);
  double*  At_ele     = mxGetPr(At_ptr);
  for (mwSize k=0; k<mDIM; ++k) {
    mwIndex  Ak_start  = At_column[k  ];
    mwIndex  Ak_end    = At_column[k+1];
    int currentSdpCone = 0;
    for (mwSize index = Ak_start; index<Ak_end; ++index) {
      mwSize Ak_j =  At_row[index];
      /* F_{k+1} = A_k */
      double ele  = At_ele[index];
      if (Ak_j < K_l) {
	fprintf(fp, PRINTF_INT_STYLE " 1 " PRINTF_INT_STYLE " " PRINTF_INT_STYLE " ", k+1, Ak_j+1, Ak_j+1);
	fprintf(fp, accuracy, ele);
	fprintf(fp, "\n");
      }
      else {
	Ak_j -= K_l;
	while (K_sdpConeStart[currentSdpCone+1] <= Ak_j) {
	  currentSdpCone++;
	}
	mwSize index2, i,j;
	index2 = Ak_j - K_sdpConeStart[currentSdpCone];
	i = index2 / K_s[currentSdpCone];
	j = index2 % K_s[currentSdpCone];
	if (i <= j) {
	  // Only upper triangular is input
	  #if 0
	  mexPrintf("input k=%d, l=%d, i=%d, j=%d, v=%e\n",
		    k+1, isK_l + currentSdpCone+1,i+1, j+1, ele);
	  #endif
	  fprintf(fp, PRINTF_INT_STYLE " " PRINTF_INT_STYLE " " PRINTF_INT_STYLE " " PRINTF_INT_STYLE " ", k+1,
		  (mwSize)(isK_l + currentSdpCone+1), i+1, j+1);
	  fprintf(fp, accuracy, ele);
	  fprintf(fp, "\n");
	}
      }
    }
  }
    
  fclose(fp);
  mxFree(filename);
  mxFree(accuracy);
  if (K_sdpNoCones > 0) {
    mxFree(K_s);
    mxFree(K_sdpConeStart);
  }

  return;
}
