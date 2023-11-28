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
#include <sdpa_call.h>

extern void _main();

void sdpasolver(double  *mDIM_ptr,
		double  *nBLOCK_ptr,
		double  *bLOCKsTRUCT_ptr,
		double  *c_ptr,
		mxArray *F_ptr,
		double  *x0_ptr,
		mxArray *X0_ptr,
		mxArray *Y0_ptr,
		mxArray *OPTION_ptr,
		int     IniPt,
		double  *objVal_ptr,
		double  *x_ptr,
		mxArray *X_ptr,
		mxArray *Y_ptr,
		mxArray *INFO_ptr)
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
  mxArray *cell_ptr  = NULL;
  
  int nSymmChk = 0;
  int nDimacs  = 0;
  
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
  mwIndex i,j,k,l;
  int size;
  mwSize mwsize;
  mwIndex idx, startidx, endidx;
  mwSize sizeM, sizeN;
  mwIndex *subscript = NULL;
  mwIndex *dims = NULL;
  int cell_index = 0;
  mwIndex *ir_ptr = NULL;
  mwIndex *jc_ptr = NULL;
  double *tmp_ptr = NULL;
  double *result_ptr = NULL;
  char* tmpPrint = NULL;
 
  /*** Set SDPA parameters by OPTIONS ***/
  /* Max Iteration */
  field_ptr = mxGetField(OPTION_ptr, 0, "maxIteration");
  if( field_ptr != NULL ){
    maxIteration = (int)mxGetScalar(field_ptr);
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
      fpResult = fopen(outfile, "w");
      if( fpResult == NULL ){
	mexPrintf("Failed to open %s\n", outfile);
	mexPrintf("Skip the detail file\n");
      } else {
	sdpa.setResultFile(fpResult);
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
  mDIM = (mwSize)(*mDIM_ptr);
  sdpa.inputConstraintNumber(mDIM);
  /* nBLOCK */
  nBLOCK = (mwSize)(*nBLOCK_ptr);
  sdpa.inputBlockNumber(nBLOCK);
  /* bLOCKsTRUCT */
  for(i = 0; i < nBLOCK; i++){
    int bs = (int)(bLOCKsTRUCT_ptr[i]);
    sdpa.inputBlockSize(i+1,bs);
    if (bs < 0 || bs == 1) {
      sdpa.inputBlockType(i+1,SDPA::LP);
    }
    else {
      sdpa.inputBlockType(i+1,SDPA::SDP);
    }
  }
  /* Execute initializeUpperTriangleSpace() */
  sdpa.initializeUpperTriangleSpace();
  /* cVECT */
  for(i = 0; i < mDIM; i++){
    sdpa.inputCVec(i+1, c_ptr[i]);
  }
  
  /*** Count NonZeroNumber in coefficience matrices ***/
  // Do nothing for SDPA 7 
  
  /*** Set coefficience matrices value ***/
  subscript = (mwIndex*)mxCalloc(2,sizeof(mwIndex));
  dims = (mwIndex*)mxGetDimensions(F_ptr);
  for(l = 0; l < dims[0]; l++){
    for(k = 0; k < dims[1]; k++){
      subscript[0]=l; subscript[1]=k;
      cell_index = mxCalcSingleSubscript(F_ptr, 2, subscript);
      cell_ptr = mxGetCell(F_ptr, cell_index);
      // mxGetDimensions at the next line is redundant,
      // but this is necessary to octave
      dims = (mwIndex*)mxGetDimensions(F_ptr);
      if( cell_ptr == NULL || mxIsEmpty(cell_ptr) ){
	/* If the cell is empty, we assume this cell as zero matrix */
	continue;
      }
      sizeM = mxGetM(cell_ptr);
      sizeN = mxGetN(cell_ptr);
      tmp_ptr = mxGetPr(cell_ptr);
      if( mxIsSparse(cell_ptr) ){
	/* Sparse Matrix */
	ir_ptr = mxGetIr(cell_ptr);
	jc_ptr = mxGetJc(cell_ptr);
	if( sizeM == 1 && sizeN != 1 ){
	  /* Row Vector Case */
	  for(j = 0; j < sizeN; j++){
	    startidx = jc_ptr[j];
	    endidx   = jc_ptr[j+1];
	    if( startidx == endidx ){ continue; }
	    sdpa.inputElement( k, l+1, j+1, j+1, tmp_ptr[startidx]);
	  }
	} else if( sizeM != 1 && sizeN == 1 ){
	  /* Column Vector Case */
	  endidx = jc_ptr[sizeN];
	  for(idx = 0; idx < endidx; idx++){
	    i = ir_ptr[idx];
	    sdpa.inputElement( k, l+1, i+1, i+1, tmp_ptr[idx]);
	  }
	} else {
	  /* Matrix Case */
	  for(j = 0; j < sizeN; j++){
	    startidx = jc_ptr[j];
	    endidx   = jc_ptr[j+1];
	    if( startidx == endidx ){ continue; }
	    for(idx = startidx; idx < endidx; idx++){
	      i = ir_ptr[idx];
	      if( i <= j ){ 
		sdpa.inputElement( k, l+1, i+1, j+1, tmp_ptr[idx]);
	      }
	    }
	  }
	}
      } else {
	/* Dense Matrix */
	if( sdpa.getBlockType(l+1)==SDPA::SDP ){
	  /* Full Matrix */
	  for(j = 0; j < sizeN; j++){
	    for(i = 0; i < sizeM; i++){
	      if( i <= j ){
		if( tmp_ptr[j * sizeN + i] != 0 ){
		  sdpa.inputElement(k, l+1, i+1, j+1,
				    tmp_ptr[j * sizeN + i]);
		}
	      }
	    }
	  }
	} else {
	  /* Diagonal Matrix */
	  if( sizeM == 1 && sizeN != 1 ){
	    /* Row Vector Case */
	    for(j = 0; j < sizeN; j++){
	      if( tmp_ptr[j] != 0 ){
		sdpa.inputElement( k, l+1, j+1, j+1, tmp_ptr[j]);
	      }
	    }
	  } else if( sizeM != 1 && sizeN == 1 ){
	    /* Column Vector Case */
	    for(i = 0; i < sizeM; i++){
	      if( tmp_ptr[i] != 0 ){
		sdpa.inputElement( k, l+1, i+1, i+1, tmp_ptr[i]);
	      }
	    }
	  } else {
	    /* Matrix Case */
	    for(j = 0; j < sizeN; j++){
	      if( tmp_ptr[j * sizeN + j] != 0 ){
		sdpa.inputElement( k, l+1, j+1, j+1, 
				  tmp_ptr[j * sizeN + j]);
	      }
	    }
	  }
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
  
  /*** Check initial point ***/
  if( IniPt == 1 ){
    sdpa.setInitPoint(true);
    /* initial value for xVec */
    for(k = 0; k < mDIM; k++){
      sdpa.inputInitXVec( k+1, x0_ptr[k]);
    }
    /* initial value for XMat */
    cell_index = 0;
    for(l = 0; l < nBLOCK; l++){
      cell_ptr = mxGetCell(X0_ptr, cell_index);
      if( cell_ptr == NULL || mxIsEmpty(cell_ptr) ){
	/* If the cell is empty, we assume this cell as zero matrix */
	continue; 
      }
      sizeM = mxGetM(cell_ptr);
      sizeN = mxGetN(cell_ptr);
      tmp_ptr = mxGetPr(cell_ptr);
      if( mxIsSparse(cell_ptr) ){
	/* Sparse Matrix */
	ir_ptr = mxGetIr(cell_ptr);
	jc_ptr = mxGetJc(cell_ptr);
	if( sizeM == 1 && sizeN != 1 ){
	  /* Row Vector Case */
	  for(j = 0; j < sizeN; j++){
	    startidx = jc_ptr[j];
	    endidx   = jc_ptr[j+1];
	    if( startidx == endidx ){ continue; }
	    sdpa.inputInitXMat(l+1, j+1, j+1, tmp_ptr[startidx]);
	  }
	} else if( sizeM != 1 && sizeN == 1 ){
	  /* Column Vector Case */
	  endidx = jc_ptr[sizeN];
	  for(idx = 0; idx < endidx; idx++){
	    i = ir_ptr[idx];
	    sdpa.inputInitXMat(l+1, i+1, i+1, tmp_ptr[idx]);
	  }
	} else {
	  /* Matrix Case */
	  for(j = 0; j < sizeN; j++){
	    startidx = jc_ptr[j];
	    endidx   = jc_ptr[j+1];
	    if( startidx == endidx ){ continue; }
	    for(idx = startidx; idx < endidx; idx++){
	      i = ir_ptr[idx];
	      if( i <= j ){
		sdpa.inputInitXMat( l+1, i+1, j+1, tmp_ptr[idx]);
	      }
	    }
	  }
	}
      } else {
	/* Dense Matrix */
	if( sdpa.getBlockType(l+1) == SDPA::SDP ){
	  /* Full Matrix */
	  for(j = 0; j < sizeN; j++){
	    for(i = 0; i < sizeM; i++){
	      if( i <= j ){
		if( tmp_ptr[j * sizeN + i] != 0 ){
		  sdpa.inputInitXMat( l+1, i+1, j+1,
				      tmp_ptr[j * sizeN + i]);
		}
	      }
	    }
	  }
	} else {
	  /* Diagonal Matrix */
	  if( sizeM == 1 && sizeN != 1 ){
	    /* Row Vector Case */
	    for(j = 0; j < sizeN; j++){
	      if( tmp_ptr[j] != 0 ){
		sdpa.inputInitXMat( l+1, j+1, j+1, tmp_ptr[j]);
	      }
	    }
	  } else if( sizeM != 1 && sizeN == 1 ){
	    /* Column Vector Case */
	    for(i = 0; i < sizeM; i++){
	      if( tmp_ptr[i] != 0 ){
		sdpa.inputInitXMat( l+1, i+1, i+1, tmp_ptr[i]);
	      }
	    }
	  } else {
	    /* Matrix Case */
	    for(j = 0; j < sizeN; j++){
	      if( tmp_ptr[j * sizeN + j] != 0 ){
		sdpa.inputInitXMat( l+1, j+1, j+1, 
				   tmp_ptr[j * sizeN + j]);
	      }
	    }
	  }
	}
      }
      cell_index++;
    }
    /* initial value for YMat */
    cell_index = 0;
    for(l = 0; l < nBLOCK; l++){
      cell_ptr = mxGetCell(Y0_ptr, cell_index);
      if( cell_ptr == NULL || mxIsEmpty(cell_ptr) ){
	/* If the cell is empty, we assume this cell as zero matrix */
	continue; 
      }
      sizeM = mxGetM(cell_ptr);
      sizeN = mxGetN(cell_ptr);
      tmp_ptr = mxGetPr(cell_ptr);
      if( mxIsSparse(cell_ptr) ){
	/* Sparse Matrix */
	ir_ptr = mxGetIr(cell_ptr);
	jc_ptr = mxGetJc(cell_ptr);
	if( sizeM == 1 && sizeN != 1 ){
	  /* Row Vector Case */
	  for(j = 0; j < sizeN; j++){
	    startidx = jc_ptr[j];
	    endidx   = jc_ptr[j+1];
	    if( startidx == endidx ){ continue; }
	    sdpa.inputInitYMat( l+1, j+1, j+1, tmp_ptr[startidx]);
	  }
	} else if( sizeM != 1 && sizeN == 1 ){
	  /* Column Vector Case */
	  endidx = jc_ptr[sizeN];
	  for(idx = 0; idx < endidx; idx++){
	    i = ir_ptr[idx];
	    sdpa.inputInitYMat( l+1, i+1, i+1, tmp_ptr[idx]);
	  }
	} else {
	  /* Matrix Case */
	  for(j = 0; j < sizeN; j++){
	    startidx = jc_ptr[j];
	    endidx   = jc_ptr[j+1];
	    if( startidx == endidx ){ continue; }
	    for(idx = startidx; idx < endidx; idx++){
	      i = ir_ptr[idx];
	      if( i <= j ){
		sdpa.inputInitYMat( l+1, i+1, j+1, tmp_ptr[idx]);
	      }
	    }
	  }
	}
      } else {
	/* Dense Matrix */
	if( sdpa.getBlockType(l+1) == SDPA::SDP ){
	  /* Full Matrix */
	  for(j = 0; j < sizeN; j++){
	    for(i = 0; i < sizeM; i++){
	      if( i <= j ){
		if( tmp_ptr[j * sizeN + i] != 0 ){
		  sdpa.inputInitYMat( l+1, i+1, j+1,
				     tmp_ptr[j * sizeN + i]);
		}
	      }
	    }
	  }
	} else {
	  /* Diagonal Matrix */
	  if( sizeM == 1 && sizeN != 1 ){
	    /* Row Vector Case */
	    for(j = 0; j < sizeN; j++){
	      if( tmp_ptr[j] != 0 ){
		sdpa.inputInitYMat( l+1, j+1, j+1, tmp_ptr[j]);
	      }
	    }
	  } else if( sizeM != 1 && sizeN == 1 ){
	    /* Column Vector Case */
	    for(i = 0; i < sizeM; i++){
	      if( tmp_ptr[i] != 0 ){
		sdpa.inputInitYMat( l+1, i+1, i+1, tmp_ptr[i]);
	      }
	    }
	  } else {
	    /* Matrix Case */
	    for(j = 0; j < sizeN; j++){
	      if( tmp_ptr[j * sizeN + j] != 0 ){
		sdpa.inputInitYMat( l+1, j+1, j+1, 
				    tmp_ptr[j * sizeN + j]);
	      }
	    }
	  }
	}
      }
      cell_index++;
    }
  }
  /*** Solve SDP ***/
  sdpa.initializeSolve();
  sdpa.solve();
  /*** Dimacs Error Information ****/
  if (nDimacs != 0) {
    field_ptr = mxCreateNumericMatrix(6,1,mxDOUBLE_CLASS,mxREAL);
    double dimacs_error[7];
    sdpa.getDimacsError(dimacs_error);
    double* dimacs_store = mxGetPr(field_ptr);
    for (int i=1; i<=6; i++) {
      dimacs_store[i-1] = dimacs_error[i];
    }
    mxSetField(INFO_ptr, 0, "dimacs", field_ptr);
  }
  /**** Set output values to arguments ****/
  /* Optimal value of Primal objective */
  objVal_ptr[0] = sdpa.getPrimalObj();
  /* Optimal value of Dual objective */
  objVal_ptr[1] = sdpa.getDualObj();
  /* Optimal value for xVec */
  tmp_ptr = sdpa.getResultXVec();
  if( tmp_ptr != NULL ){
    for(k = 0; k < mDIM; k++){
      x_ptr[k] = tmp_ptr[k];
    }
  }
  /* Optimal value for XMat */
  cell_index = 0;
  for(l = 0; l < nBLOCK; l++){
    size = sdpa.getBlockSize(l+1);
    if( sdpa.getBlockType(l+1) == SDPA:: SDP){
      sizeM = size;
      sizeN = size;
    } else {
      sizeM = 1;
      sizeN = abs(size);
    }
    cell_ptr = mxCreateDoubleMatrix(sizeM, sizeN, mxREAL);
    tmp_ptr = mxGetPr(cell_ptr);
    idx = 0;
    result_ptr = sdpa.getResultXMat(l+1);
    if( size >= 0 ){
      for(j = 0; j < sizeN; j++){
	for(i = 0; i < sizeM; i++){
	  tmp_ptr[idx++] = result_ptr[j + sizeN * i];
	}
      }
    } else {
      for(idx = 0; idx < sizeN; idx++){
	tmp_ptr[idx] = result_ptr[idx];
      }
    }
    mxSetCell(X_ptr, cell_index++, mxDuplicateArray(cell_ptr));
  }
  /* Optimal value for YMat */
  cell_index = 0;
  for(l = 0; l < nBLOCK; l++){
    size = sdpa.getBlockSize(l+1);
    if( sdpa.getBlockType(l+1) == SDPA:: SDP){
      sizeM = size;
      sizeN = size;
    } else {
      sizeM = 1;
      sizeN = abs(size);
    }
    cell_ptr = mxCreateDoubleMatrix(sizeM, sizeN, mxREAL);
    tmp_ptr = mxGetPr(cell_ptr);
    idx = 0;
    result_ptr = sdpa.getResultYMat(l+1);
    if( size >= 0 ){
      for(j = 0; j < sizeN; j++){
	for(i = 0; i < sizeM; i++){
	  tmp_ptr[idx++] = result_ptr[j + sizeN * i];
	}
      }
    } else {
      for(idx = 0; idx < sizeN; idx++){
	tmp_ptr[idx] = result_ptr[idx];
      }
    }
    mxSetCell(Y_ptr, cell_index++, mxDuplicateArray(cell_ptr));
  }
  /* Phase information */
  field_ptr = mxCreateString(szPhase[sdpa.getPhaseValue()]);
  mxSetField(INFO_ptr, 0, "phasevalue", field_ptr);
  /* Iteration */
  field_ptr = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
  *mxGetPr(field_ptr) = (double)sdpa.getIteration();
  mxSetField(INFO_ptr, 0, "iteration", field_ptr);
  /* close output file */
  time(&ltime);
  strcpy(string_time,ctime(&ltime));
  string_time[strlen(string_time)-1]='\0';
  if (fp) {
    fprintf(fp,"SDPA end at [%s]\n",string_time);
  }
  if (fpResult) {
    fprintf(fpResult,"SDPA end at [%s]\n",string_time);
  }

  if( nOutfile ){
    fclose(fp);
  }
  if (fpResult != NULL) {
    fclose(fpResult);
  }
  /*** Free allocated memory ****/
  mxFree(subscript);
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
  
  double
    *mDIM_ptr,
    *nBLOCK_ptr,
    *bLOCKsTRUCT_ptr,
    *c_ptr,
    *x_ptr,
    *objVal_ptr,
    *x0_ptr;

  mxArray
    *F_ptr,
    *X_ptr,
    *Y_ptr,
    *X0_ptr,
    *Y0_ptr,
    *OPTION_ptr,
    *INFO_ptr;
  
  int IniPt;

  const char *fnames[] = {
    "phasevalue",
    "iteration",
    "dimacs",
    "cpusec"
  };
  
  /* 
   * check arguments
   */
  if( nrhs != 9 ){
    mexErrMsgTxt("Input arguments must be 9.");
  }
  
  /* Get the pointer of input variables */
  mDIM_ptr        = mxGetPr(prhs[0]);
  nBLOCK_ptr      = mxGetPr(prhs[1]);
  bLOCKsTRUCT_ptr = mxGetPr(prhs[2]);
  c_ptr           = mxGetPr(prhs[3]);
  F_ptr           = (mxArray*)prhs[4];
  x0_ptr          = mxGetPr(prhs[5]);
  X0_ptr          = (mxArray*)prhs[6];
  Y0_ptr          = (mxArray*)prhs[7];
  OPTION_ptr      = (mxArray*)prhs[8];
  
  if( mxIsEmpty(prhs[5])    /* x0 */
      || mxIsEmpty(prhs[6]) /* X0 */
      || mxIsEmpty(prhs[7]) /* Y0 */
      ){
    /* nouse of Initial Point*/
    IniPt = 0;
  } else {
    /* use of Initial Point */
    IniPt = 1;
    /* check of argument dimensions */
    /* if(*mDIM_ptr != mxGetM(prhs[5]) || !mxIsDouble(prhs[5]))
       mexErrMsgTxt("x0 must be (mDIM x 1) column vector of double");
       if( *nBLOCK_ptr != sizeM_X0 || 1 != sizeN_X0 || !mxIsCell(prhs[6]))
       mexErrMsgTxt("X0 must be (nBLOCK x 1) cell array");
       if( *nBLOCK_ptr != sizeM_Y0 || 1 != sizeN_Y0 || !mxIsCell(prhs[7]))
       
       mexErrMsgTxt("Y0 must be (nBLOCK x 1) cell array");*/
  }
  /* Create cellarrays for the output variables */
  plhs[0] = mxCreateDoubleMatrix(1,2,mxREAL);
  plhs[1] = mxCreateDoubleMatrix((int) *mDIM_ptr,1,mxREAL);
  plhs[2] = mxCreateCellMatrix((int) *nBLOCK_ptr,1);
  plhs[3] = mxCreateCellMatrix((int) *nBLOCK_ptr,1);
  plhs[4] = mxCreateStructMatrix(1,1,4,fnames);
  
  //Get the pointer of output variables
  objVal_ptr = mxGetPr(plhs[0]);
  x_ptr      = mxGetPr(plhs[1]);
  X_ptr      = plhs[2];
  Y_ptr      = plhs[3];
  INFO_ptr   = plhs[4];

  /* Call sdpasolver here */
  sdpasolver(mDIM_ptr,
	     nBLOCK_ptr,
	     bLOCKsTRUCT_ptr,
	     c_ptr,
	     F_ptr,
	     x0_ptr,
	     X0_ptr,
	     Y0_ptr,
	     OPTION_ptr,
	     IniPt,
	     objVal_ptr,
	     x_ptr,
	     X_ptr,
	     Y_ptr,
	     INFO_ptr);
   
  return;
}

/*
 * End of File
 */
