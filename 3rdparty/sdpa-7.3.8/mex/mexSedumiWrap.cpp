/*
 * This file is a component of SDPA
 * Copyright (C) 2004-2012 SDPA Project
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
 *
 * SDPA-M: 7.3
 * mexsdpa.cpp
 */

#include <mex.h>
/*
 * SDPA header files
 */
#include <sdpa_tool.h>
#include <sdpa_call.h>
using namespace sdpa;

extern void _main();

void sdpasolver(mxArray* At_ptr, mxArray* b_ptr, mxArray* c_ptr,
		mxArray* K_ptr, mxArray* OPTION_ptr,
		mxArray* x_ptr, mxArray* y_ptr,
		mxArray* info_ptr)
{
  time_t ltime;
  time(&ltime);
  char string_time[1024];
  strcpy(string_time,ctime(&ltime));
  string_time[strlen(string_time)-1]='\0';

  SDPA sdpa;
  
  int    maxIteration = 0;
  double param        = 0.0;
  /* mxArray pointer */
  mxArray *field_ptr = NULL;
  int nSymmChk = 0;
  int nDimacs = 0;
  
  /* strings for phase value */
  const char *szPhase[] = {
    "noINFO", "pFEAS", "dFEAS", "pdFEAS", "pdINF",
    "pFEAS_dINF", "pINF_dFEAS", "pdOPT", "pUNBD", "dUNBD"};
  
  /* output file */
  char *outfile = NULL;
  FILE *fp = NULL;
  FILE *fpResult = NULL;
  int nOutfile = 0;

  mwSize mDIM;
  mwSize nBLOCK;
  
  /* temporary variables */
  mwIndex k;
  int size;
  mwSize mwsize;
  double *tmp_ptr = NULL;
  char* tmpPrint = NULL;

  TimeStart(SDPA_START);
  TimeStart(SDPA_CONVERT_START);
  
  /*** Set SDPA parameters by OPTIONS ***/
  /* Max Iteration */
  field_ptr = mxGetField(OPTION_ptr, 0, "maxIteration");
  if( field_ptr != NULL ){
    maxIteration = (int)mxGetScalar(field_ptr);
    // mexPrintf("maxIteration = %d\n",maxIteration);
    sdpa.setParameterMaxIteration(maxIteration);
  }
  /* epsilonStar */
  field_ptr = mxGetField(OPTION_ptr, 0, "epsilonStar");
  if( field_ptr != NULL ){
    param = *mxGetPr(field_ptr);
    sdpa.setParameterEpsilonStar(param);
  }
  /* lambdaStar */
  field_ptr = mxGetField(OPTION_ptr, 0, "lambdaStar");
  if( field_ptr != NULL ){
    param = *mxGetPr(field_ptr);
    sdpa.setParameterLambdaStar(param);
  }
  /* omegaStar */
  field_ptr = mxGetField(OPTION_ptr, 0, "omegaStar");
  if( field_ptr != NULL ){
    param = *mxGetPr(field_ptr);
    sdpa.setParameterOmegaStar(param);
  }
  /* lowerBound */
  field_ptr = mxGetField(OPTION_ptr, 0, "lowerBound");
  if( field_ptr != NULL ){
    param = *mxGetPr(field_ptr);
    sdpa.setParameterLowerBound(param);
  }
  /* upperBound */
  field_ptr = mxGetField(OPTION_ptr, 0, "upperBound");
  if( field_ptr != NULL ){
    param = *mxGetPr(field_ptr);
    sdpa.setParameterUpperBound(param);
  }
  /* betaStar */
  field_ptr = mxGetField(OPTION_ptr, 0, "betaStar");
  if( field_ptr != NULL ){
    param = *mxGetPr(field_ptr);
    sdpa.setParameterBetaStar(param);
  }
  /* betaBar */
  field_ptr = mxGetField(OPTION_ptr, 0, "betaBar");
  if( field_ptr != NULL ){
    param = *mxGetPr(field_ptr);
    sdpa.setParameterBetaBar(param);
  }
  /* gammaStar */
  field_ptr = mxGetField(OPTION_ptr, 0, "gammaStar");
  if( field_ptr != NULL ){
    param = *mxGetPr(field_ptr);
    sdpa.setParameterGammaStar(param);
  }
  /* epsilonDash */
  field_ptr = mxGetField(OPTION_ptr, 0, "epsilonDash");
  if( field_ptr != NULL ){
    param = *mxGetPr(field_ptr);
    sdpa.setParameterEpsilonDash(param);
  }
  /* xPrint */
  field_ptr = mxGetField(OPTION_ptr, 0, "xPrint");
  if( field_ptr != NULL ){
    mwsize = mxGetM(field_ptr) * mxGetN(field_ptr) + 1;
    tmpPrint = (char*)mxCalloc(mwsize, sizeof(char));
    mxGetString(field_ptr, tmpPrint, mwsize);
    sdpa.setParameterPrintXVec(tmpPrint);
    mxFree(tmpPrint);
  }
  /* XPrint */
  field_ptr = mxGetField(OPTION_ptr, 0, "XPrint");
  if( field_ptr != NULL ){
    mwsize = mxGetM(field_ptr) * mxGetN(field_ptr) + 1;
    tmpPrint = (char*)mxCalloc(mwsize, sizeof(char));
    mxGetString(field_ptr, tmpPrint, mwsize);
    sdpa.setParameterPrintXMat(tmpPrint);
    mxFree(tmpPrint);
  }
  /* YPrint */
  field_ptr = mxGetField(OPTION_ptr, 0, "YPrint");
  if( field_ptr != NULL ){
    mwsize = mxGetM(field_ptr) * mxGetN(field_ptr) + 1;
    tmpPrint = (char*)mxCalloc(mwsize, sizeof(char));
    mxGetString(field_ptr, tmpPrint, mwsize);
    sdpa.setParameterPrintYMat(tmpPrint);
    mxFree(tmpPrint);
  }
  /* infPrint */
  field_ptr = mxGetField(OPTION_ptr, 0, "infPrint");
  if( field_ptr != NULL ){
    mwsize = mxGetM(field_ptr) * mxGetN(field_ptr) + 1;
    tmpPrint = (char*)mxCalloc(mwsize, sizeof(char));
    mxGetString(field_ptr, tmpPrint, mwsize);
    sdpa.setParameterPrintInformation(tmpPrint);
    mxFree(tmpPrint);
  }


  /* isSymmetric */
  field_ptr = mxGetField(OPTION_ptr, 0, "isSymmetric");
  if( field_ptr != NULL ){
    nSymmChk = (int)mxGetScalar(field_ptr);
  }
  /* isDimacs */
  field_ptr = mxGetField(OPTION_ptr, 0, "isDimacs");
  if( field_ptr != NULL ){
    nDimacs = (int)mxGetScalar(field_ptr);
  }
  /* print */
  field_ptr = mxGetField(OPTION_ptr, 0, "print");
  if( field_ptr != NULL ){
    mwsize = mxGetM(field_ptr) * mxGetN(field_ptr) + 1;
    if (mwsize == 1) {
      // mexPrintf("display is NULL\n");
      fp = NULL;
    }
    else {
      outfile = (char*)mxCalloc(mwsize, sizeof(char));
      mxGetString(field_ptr, outfile, mwsize);
      if( strncmp("display", outfile, mwsize - 1) == 0 ){
	fp = stdout;
      } else if( strncmp("no", outfile, mwsize - 1) == 0 ){
	fp = NULL;
      } else {
	fp = fopen(outfile, "at");
	if( fp == NULL ){
	  mexPrintf("Failed to open %s\n", outfile);
	  fp = stdout;
	} else {
	  nOutfile = 1;
	}
      }
      mxFree(outfile);
    }
  } else {
    /* default setting is displaying information to stdout */
    fp = stdout;
  }
  sdpa.setDisplay(fp);
  
  /* resultFile */
  field_ptr = mxGetField(OPTION_ptr, 0, "resultFile");
  if( field_ptr != NULL ){
    mwsize = mxGetM(field_ptr) * mxGetN(field_ptr) + 1;
    if (mwsize == 1) {
      // mexPrintf("resultFile is NULL\n");
    }
    else {
      outfile = (char*)mxCalloc(mwsize, sizeof(char));
      mxGetString(field_ptr, outfile, mwsize);
      if ( strncmp("no", outfile, mwsize - 1) == 0 ) {
	mexPrintf("resultFile is NULL\n");
      }
      else {
	fpResult = fopen(outfile, "w");
	if ( fpResult == NULL ) {
	  mexPrintf("Failed to open %s\n", outfile);
	  mexPrintf("Skip the detail file\n");
	} else {
	  sdpa.setResultFile(fpResult);
	}
      }
      mxFree(outfile);
    }
  }
  
  if (fp) {
    fprintf(fp,"SDPA start at [%s]\n",string_time);
  }
  if (fpResult) {
    fprintf(fpResult,"SDPA start at [%s]\n",string_time);
  }

  /* NumThreads */
  field_ptr = mxGetField(OPTION_ptr, 0, "NumThreads");
  if( field_ptr != NULL ){
    sdpa.setNumThreads((int)mxGetScalar(field_ptr));
  }
  
  /*** initialize SDPA class members ***/
  /* mDIM */
  mDIM = mxGetN(At_ptr);
  sdpa.inputConstraintNumber(mDIM);
  /* nBLOCK */
  nBLOCK = 0;
  mwSize K_l   = 0;
  int    isK_l = 0; // 1 (K_l > 0) or 0 (K_l == 0)
  field_ptr = mxGetField(K_ptr, 0, "l");
  if (field_ptr != NULL) {
    K_l   = (mwSize)((mxGetPr(field_ptr))[0]);
    #if 0
    mexPrintf("K_l = %zd\n", K_l);
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
    for (int l=0; l<K_sdpNoCones; ++l) {
      K_s[l] = (mwSize)((mxGetPr(field_ptr))[l]);
      #if 0
      mexPrintf("K_s[%zd] = %zd\n", l+isK_l+1, K_s[l]);
      #endif
      K_sdpConeStart[l+1] = K_sdpConeStart[l] + K_s[l]*K_s[l];
      nBLOCK++;
    }
  }

  #if 0
  mexPrintf("isK_l = %d, nBlock = %d \n", isK_l, nBLOCK);
  #endif
  sdpa.inputBlockNumber(nBLOCK);
  if (isK_l > 0) {
    sdpa.inputBlockSize(1, K_l);
    sdpa.inputBlockType(1, SDPA::LP);
  }
  for (int l=0; l<K_sdpNoCones; ++l) {
    #if 0
    mexPrintf("l+isK_l+1 = %d, K_s = %d \n", l+isK_l+1, K_s[l]);
    #endif
    sdpa.inputBlockSize(l+isK_l+1, K_s[l]);
    sdpa.inputBlockType(l+isK_l+1, SDPA::SDP);
  }

  /* Execute initializeUpperTriangleSpace() */
  sdpa.initializeUpperTriangleSpace();

  /* cVECT = -b*/
  double* b = mxGetPr(b_ptr);
  for(mwSize i = 0; i < mDIM; i++){
    sdpa.inputCVec((int)i+1, -b[i]);
  }
  
  /*** Count NonZeroNumber in coefficience matrices ***/
  // Do nothing for SDPA 7 
  
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
      /* A_0 = - C */
      double ele = -C_ele[index];
      // mexPrintf("C_j = %zd, ele = %e\n", C_j, ele);
      if (C_j < K_l) {
	sdpa.inputElement(0, 1, C_j+1, C_j+1, ele);
      }
      else {
	C_j -= K_l;
	while (K_sdpConeStart[currentSdpCone+1] <= C_j) {
	  currentSdpCone++;
	}
	int index2, i,j;
	index2 = C_j - K_sdpConeStart[currentSdpCone];
	i = index2 / K_s[currentSdpCone];
	j = index2 % K_s[currentSdpCone];
	if (i <= j) {
	  // Only upper triangular is input
	  #if 0
	  mexPrintf("input k=%d, l=%d, i=%d, j=%d, v=%e\n",
		    0, isK_l + currentSdpCone+1,i+1, j+1, ele);
	  #endif
	  sdpa.inputElement(0, isK_l + currentSdpCone+1,i+1, j+1, ele);
	}
      }
    }
  }
  
  /* F_{k+1} = - A_k */
  mwIndex* At_row     = mxGetIr(At_ptr);
  mwIndex* At_column  = mxGetJc(At_ptr);
  double*  At_ele     = mxGetPr(At_ptr);
  for (k=0; k<mDIM; ++k) {
    mwIndex  Ak_start  = At_column[k  ];
    mwIndex  Ak_end    = At_column[k+1];
    int currentSdpCone = 0;
    for (mwSize index = Ak_start; index<Ak_end; ++index) {
      mwSize Ak_j =  At_row[index];
      /* F_{k+1} = - A_k */
      double ele  = -At_ele[index];
      if (Ak_j < K_l) {
	sdpa.inputElement(k+1, 1, Ak_j+1, Ak_j+1, ele);
      }
      else {
	Ak_j -= K_l;
	while (K_sdpConeStart[currentSdpCone+1] <= Ak_j) {
	  currentSdpCone++;
	}
	int index2, i,j;
	index2 = Ak_j - K_sdpConeStart[currentSdpCone];
	i = index2 / K_s[currentSdpCone];
	j = index2 % K_s[currentSdpCone];
	if (i <= j) {
	  // Only upper triangular is input
	  #if 0
	  mexPrintf("input k=%d, l=%d, i=%d, j=%d, v=%e\n",
		    k+1, isK_l + currentSdpCone+1,i+1, j+1, ele);
	  #endif
	  sdpa.inputElement(k+1, isK_l + currentSdpCone+1,i+1, j+1, ele);
	}
      }
    }
  }

  /*** Check the consistence of F, c ***/
  if( nSymmChk ){
    sdpa.initializeUpperTriangle(true);
  }
  else {
    sdpa.initializeUpperTriangle(false);
  }

  // sdpa.writeInputSparse((char*)"b.dat-s",(char*)"%e");

  /*** Solve SDP ***/
  sdpa.initializeSolve();
  mexPrintf("Converted to SDPA internal data / ");
  mexPrintf("Starting SDPA main loop\n");
  TimeEnd(SDPA_CONVERT_END);
  TimeStart(SDPA_SOLVE_START);
  sdpa.solve();
  TimeEnd(SDPA_SOLVE_END);
  TimeStart(SDPA_RETRIEVE_START);
  mexPrintf("Converting optimal solution to Sedumi format\n");
  /**** Set output values to arguments ****/
  
  /* Optimal value for xVec */
  double* y = mxGetPr(y_ptr);
  tmp_ptr = sdpa.getResultXVec();
  if( tmp_ptr != NULL ){
    for(k = 0; k < mDIM; k++){
      y[k] = tmp_ptr[k];
    }
  }

  /* Optimal value for YMat */
  double* x = mxGetPr(x_ptr);
  if (isK_l > 0) {
    size = sdpa.getBlockSize(1);
    tmp_ptr = sdpa.getResultYMat(1);
    for (int index = 0; index < size; ++index) {
      x[index] = tmp_ptr[index];
    }
  }
  for (int l=0; l<K_sdpNoCones; ++l) {
    size = sdpa.getBlockSize(l+isK_l+1);
    tmp_ptr = sdpa.getResultYMat(l+isK_l+1);
    for (int index = 0; index < size*size; ++index) {
      x[K_l + K_sdpConeStart[l] + index] = tmp_ptr[index];
    }
  }
  TimeEnd(SDPA_RETRIEVE_END);

  /* Dimacs Error Information */
  if (nDimacs != 0) {
    mexPrintf("Computing Dimacs Error\n");
    field_ptr = mxCreateNumericMatrix(6,1,mxDOUBLE_CLASS,mxREAL);
    double dimacs_error[7];
    sdpa.getDimacsError(dimacs_error);
    double* dimacs_store = mxGetPr(field_ptr);
    for (int i=1; i<=6; i++) {
      dimacs_store[i-1] = dimacs_error[i];
    }
    mxSetField(info_ptr, 0, "dimacs", field_ptr);
  }

  
  /* Phase information */
  field_ptr = mxCreateString(szPhase[sdpa.getPhaseValue()]);
  mxSetField(info_ptr, 0, "phasevalue", field_ptr);
  /* Iteration */
  field_ptr = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
  *mxGetPr(field_ptr) = (double)sdpa.getIteration();
  mxSetField(info_ptr, 0, "iteration", field_ptr);

  /* primalObj */
  field_ptr = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
  *mxGetPr(field_ptr) = -sdpa.getDualObj();
  mxSetField(info_ptr, 0, "primalObj", field_ptr);
  /* dualObj */
  field_ptr = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
  *mxGetPr(field_ptr) = -sdpa.getPrimalObj();
  mxSetField(info_ptr, 0, "dualObj", field_ptr);
  /* primalError */
  field_ptr = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
  *mxGetPr(field_ptr) = sdpa.getDualError();
  mxSetField(info_ptr, 0, "primalError", field_ptr);
  /* dualError */
  field_ptr = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
  *mxGetPr(field_ptr) = sdpa.getPrimalError();
  mxSetField(info_ptr, 0, "dualError", field_ptr);
  /* digits */
  field_ptr = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
  *mxGetPr(field_ptr) = sdpa.getDigits();
  mxSetField(info_ptr, 0, "digits", field_ptr);
  /* dualityGap */
  field_ptr = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
  *mxGetPr(field_ptr) = sdpa.getDualityGap();
  mxSetField(info_ptr, 0, "dualityGap", field_ptr);
  /* mu */
  field_ptr = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
  *mxGetPr(field_ptr) = sdpa.getMu();
  mxSetField(info_ptr, 0, "mu", field_ptr);

  /* solveTime */
  field_ptr = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
  *mxGetPr(field_ptr) = TimeCal(SDPA_SOLVE_START,SDPA_SOLVE_END);
  mxSetField(info_ptr, 0, "solveTime", field_ptr);
  /* convertingTime */
  field_ptr = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
  *mxGetPr(field_ptr) = TimeCal(SDPA_CONVERT_START,SDPA_CONVERT_END);
  mxSetField(info_ptr, 0, "convertingTime", field_ptr);
  /* retrivingTime */
  field_ptr = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
  *mxGetPr(field_ptr) = TimeCal(SDPA_RETRIEVE_START,SDPA_RETRIEVE_END);
  mxSetField(info_ptr, 0, "retrievingTime", field_ptr);
  /* retrivingTime */
  field_ptr = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
  TimeEnd(SDPA_END);
  *mxGetPr(field_ptr) = TimeCal(SDPA_START,SDPA_END);
  mxSetField(info_ptr, 0, "sdpaTime", field_ptr);

  time(&ltime);
  strcpy(string_time,ctime(&ltime));
  string_time[strlen(string_time)-1]='\0';
  
  if (fp) {
    fprintf(fp,"SDPA end at [%s]\n",string_time);
  }
  if (fpResult) {
    fprintf(fpResult,"SDPA end at [%s]\n",string_time);
  }
  /* close output file */
  if( nOutfile ){
    fclose(fp);
  }
  if (fpResult != NULL) {
    fclose(fpResult);
  }
  
  /*** Free allocated memory ****/
  mxFree(K_s);
  mxFree(K_sdpConeStart);
  sdpa.terminate();
  
  return;
}

/*
 * Matlab gateway function
 */
void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  /* Decleration of variables */

  mxArray* At_ptr;
  mxArray* b_ptr;
  mxArray* c_ptr;
  mxArray* K_ptr;
  mxArray* OPTION_ptr;
  
  mxArray* x_ptr;
  mxArray* y_ptr;
  mxArray* info_ptr;

  const char *fnames[] = {
    "phasevalue",
    "iteration",
    "cpusec",
    "primalObj",
    "dualObj",
    "primalError",
    "dualError",
    "digits",
    "dualityGap",
    "mu",
    "dimacs",
    "solveTime",
    "convertingTime",
    "retrievingTime",
    "sdpaTime"
  };
  
  /* Get the pointer of input variables */
  At_ptr     = (mxArray*)prhs[0];
  b_ptr      = (mxArray*)prhs[1];
  c_ptr      = (mxArray*)prhs[2];
  K_ptr      = (mxArray*)prhs[3];
  OPTION_ptr = (mxArray*)prhs[4];

  mwSize m = mxGetM(b_ptr);
  mwSize n = mxGetM(c_ptr);
  
  /* Create cellarrays for the output variables */
  plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(m,1,mxREAL);
  plhs[2] = mxCreateStructMatrix(1,1,15,fnames);
  
  //Get the pointer of output variables
  x_ptr    = plhs[0];
  y_ptr    = plhs[1];
  info_ptr = plhs[2];
  
  /* Call sdpasolver here */
  sdpasolver(At_ptr,b_ptr,c_ptr,K_ptr,OPTION_ptr,
	     x_ptr,y_ptr,info_ptr);   
  return;
}

/*
 * End of File
 */
