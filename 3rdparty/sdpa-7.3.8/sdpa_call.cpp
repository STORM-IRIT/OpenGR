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
/*--------------------------------------------------
  sdpa_call.cpp
--------------------------------------------------*/

#include "sdpa_call.h"
#include "sdpa_io.h"
#include "sdpa_linear.h"

using namespace sdpa;
#define LengthOfBuffer 1024

IndexLIJv::IndexLIJv()
{
  // nothing is needed
}
		     
IndexLIJv::~IndexLIJv()
{
  // nothing is needed
}
		     

bool IndexLIJv::compare(IndexLIJv* a, IndexLIJv* b)
{
  // if  (a <  b) return true;
  // if  (a >= b) return  false;
  
  if      ( a[0].l < b[0].l ) { return true; }
  else if ( a[0].l > b[0].l ) { return false; }
  if      ( a[0].i < b[0].i ) { return true; }
  else if ( a[0].i > b[0].i ) { return false; }
  if      ( a[0].j < b[0].j ) { return true; }
  else if ( a[0].j > b[0].j ) { return false; }
  return false; // a == b
}
  


SDPA::SDPA()
{
  KAPPA      = 1.5;
  m          = 0;
  nBlock     = 0;
  fpout      = NULL;
  Display    = NULL;
  isInitPoint = false;

  typeParameter = PARAMETER_DEFAULT;
  param.setDefaultParameter(Parameter::PARAMETER_DEFAULT);
}

SDPA::~SDPA()
{
  terminate();
}

void SDPA::setParameterType(ParameterType type)
{
  if (type == PARAMETER_DEFAULT) {
    param.setDefaultParameter(Parameter::PARAMETER_DEFAULT);
  } else if (type == PARAMETER_UNSTABLE_BUT_FAST) {
    param.setDefaultParameter(Parameter::PARAMETER_UNSTABLE_BUT_FAST);
  } else if (type == PARAMETER_STABLE_BUT_SLOW) {
    param.setDefaultParameter(Parameter::PARAMETER_STABLE_BUT_SLOW);
  }
  typeParameter = type;
}

void SDPA::setParameterMaxIteration(int maxIteration)
{
  param.maxIteration = maxIteration;
}

void SDPA::setParameterEpsilonStar (double epsilonStar)
{
  param.epsilonStar = epsilonStar;
}

void SDPA::setParameterLambdaStar  (double lambdaStar)
{
  param.lambdaStar = lambdaStar;
}

void SDPA::setParameterOmegaStar   (double omegaStar)
{
  param.omegaStar = omegaStar;
}

void SDPA::setParameterLowerBound  (double lowerBound)
{
  param.lowerBound = lowerBound;
}

void SDPA::setParameterUpperBound  (double upperBound)
{
  param.upperBound = upperBound;
}
void SDPA::setParameterBetaStar    (double betaStar)
{
  param.betaStar = betaStar;
}

void SDPA::setParameterBetaBar     (double betaBar)
{
  param.betaBar = betaBar;
}

void SDPA::setParameterGammaStar   (double gammaStar)
{
  param.gammaStar = gammaStar;
}

void SDPA::setParameterEpsilonDash (double epsilonDash)
{
  param.epsilonDash = epsilonDash;
}

void SDPA::setParameterPrintXVec(char* xPrint)
{
  strncpy(param.xPrint,xPrint,PRINT_DEFAULT_LENGTH);
}

void SDPA::setParameterPrintXMat(char* XPrint)
{
  strncpy(param.XPrint,XPrint,PRINT_DEFAULT_LENGTH);
}

void SDPA::setParameterPrintYMat(char* YPrint)
{
  strncpy(param.YPrint,YPrint,PRINT_DEFAULT_LENGTH);
}

void SDPA::setParameterPrintInformation(char* infPrint)
{
  strncpy(param.infPrint,infPrint,PRINT_DEFAULT_LENGTH);
}

void SDPA::setDisplay(FILE* Display)
{
  this->Display = Display;
}

void SDPA::setResultFile(FILE* fpout)
{
  this->fpout = fpout;
}

void SDPA::setInitPoint(bool isInitPoint)
{
  if (this->isInitPoint == false && isInitPoint == true) {
    initPt_xMat.initialize(bs);
    initPt_zMat.initialize(bs);
  }
  this->isInitPoint = isInitPoint;
  if (isInitPoint == false) {
    mu.initialize(param.lambdaStar);
    currentPt.initialize(m, bs, param.lambdaStar, com);
  }
}

void SDPA::setNumThreads(int NumThreads)
{
  newton.setNumThreads(Display,fpout,NumThreads);
}

SDPA::ParameterType SDPA::getParameterType()
{
  return typeParameter;
}

int SDPA::getParameterMaxIteration()
{
  return param.maxIteration;
}
double SDPA::getParameterEpsilonStar()
{
  return param.epsilonStar;
}
double SDPA::getParameterLambdaStar()
{
  return param.lambdaStar;
}
double SDPA::getParameterOmegaStar()
{
  return param.omegaStar;
}
double SDPA::getParameterLowerBound()
{
  return param.lowerBound;
}
double SDPA::getParameterUpperBound()
{
  return param.upperBound;
}
double SDPA::getParameterBetaStar()
{
  return param.betaStar;
}
double SDPA::getParameterBetaBar()
{
  return param.betaBar;
}
double SDPA::getParameterGammaStar()
{
  return param.gammaStar;
}
double SDPA::getParameterEpsilonDash()
{
  return param.epsilonDash;
}
char* SDPA::getParameterPrintXVec()
{
  return param.xPrint;
}
char* SDPA::getParameterPrintXMat()
{
  return param.XPrint;
}
char* SDPA::getParameterPrintYMat()
{
  return param.YPrint;
}
char* SDPA::getParameterPrintInformation()
{
  return param.infPrint;
}
FILE* SDPA::getDisplay()
{
  return Display;
}
FILE* SDPA::getResultFile()
{
  return fpout;
}
bool SDPA::getInitPoint()
{
  return isInitPoint;
}

int SDPA::getNumThreads()
{
  return newton.NUM_THREADS;
}

void SDPA::inputConstraintNumber(int m)
{
  this->m = m;
}

void SDPA::inputBlockNumber(int nBlock)
{
  this->nBlock = nBlock;
  bs.initialize(nBlock);
}

void SDPA::inputBlockSize(int l, int size)
{
  bs.blockStruct[l-1] = size;
}

void SDPA::inputBlockType(int l, ConeType coneType)
{
  if (coneType == SDPA::SDP) {
    bs.blockType[l-1] = BlockStruct::btSDP;
  }
  if (coneType == SDPA::SOCP) {
    bs.blockType[l-1] = BlockStruct::btSOCP;
  }
  if (coneType == SDPA::LP) {
    bs.blockType[l-1] = BlockStruct::btLP;
  }
}

void SDPA::inputCVec(int k, double value)
{
  if (k > m || k <= 0) {
    rError("k exceeds ConstraintNumber or "
	   "k is less than or equal to zero :: m= "
	   << m << " : k= " << k);

  }
  inputData.b.ele[k-1] = value;
}

void SDPA::inputElement(int k, int l, int i, int j, double value,
			bool inputCheck)
{
  if (inputCheck) {
    if (k > m || k < 0) {
      rError ("k exceeds ConstraintNumber or "
	      "k is less than zero :: m= "
	      << m << " : k= " << k << " : l= " << l
	      << " : i= " << i << " : j= " << j);
    }
    if (l > nBlock || l <= 0) {
      rError ("l exceeds nBlock or "
	      "l is less than or equal to zero :: nBlock= "
	      << nBlock << " : k= " << k << " : l= " << l
	      << " : i= " << i << " : j= " << j);
    }
    int dim = bs.blockStruct[l-1];
    if (i > dim || i <= 0) {
      rError ("i exceeds dimension of the block or "
	      "i is less than or equal to zero :: dim= "
	      << dim << " : k= " << k << " : l= " << l
	      << " : i= " << i << " : j= " << j);
    }
    if (j > dim || j <= 0) {
      rError ("j exceeds dimension of the block or "
	      "j is less than or equal to zero :: dim= "
	      << dim << " : k= " << k << " : l= " << l
	      << " : i= " << i << " : j= " << j);
    }
    if (bs.blockType[l-1] == BlockStruct::btSDP) {
      if (i > j) {
	rMessage("Swap i and j [Only Upper Triangle]"
		 " : k= " << k << " : l= " << l
		 << " : i= " << i << " : j= " << j);
      }
    }
    if (bs.blockType[l-1] == BlockStruct::btLP) {
      if (i!=j) {
	rError("i should be j in LP block"
	       " : k= " << k << " : l= " << l
	       << " : i= " << i << " : j= " << j);
      }
    }
  }

  if (i > j) {
    int tmp = i; i = j; j = tmp;
  }
  
  IndexLIJv* indexLIJv;
  NewArray(indexLIJv,IndexLIJv,1);
  indexLIJv[0].l     = l;
  indexLIJv[0].i     = i;
  indexLIJv[0].j     = j;
  indexLIJv[0].value = value;
  NonZeroElements[k].push_back(indexLIJv);
}

void SDPA::inputInitXVec(int k, double value)
{
  if (k > m || k <= 0) {
    rError("k exceeds ConstraintNumber or "
	   "k is less than or equal to zero :: m= "
	   << m << " : k= " << k);
  }
  // Note reverse primal-dual
  currentPt.yVec.ele[k-1] = -value;
}

void SDPA::inputInitXMat(int l, int i, int j, double value)
{
  if (l > nBlock || l <= 0) {
    rError ("l exceeds nBlock or "
	    "l is less than or equal to zero :: nBlock= "
	    << nBlock << " : l= " << l
	    << " : i= " << i << " : j= " << j);
  }
  int dim = bs.blockStruct[l-1];
  if (i > dim || i <= 0) {
    rError ("i exceeds dimension of the block or "
	    "i is less than or equal to zero :: dim= "
	    << dim << " : l= " << l
	    << " : i= " << i << " : j= " << j);
  }
  if (j > dim || j <= 0) {
    rError ("j exceeds dimension of the block or "
	    "j is less than or equal to zero :: dim= "
	    << dim << " : l= " << l
	    << " : i= " << i << " : j= " << j);
  }
  if (bs.blockType[l-1] == BlockStruct::btLP) {
    if (i!=j) {
      rError("i should be j in LP block"
	     " : l= " << l
	     << " : i= " << i << " : j= " << j);
    }
  }

  if (bs.blockType[l-1] == BlockStruct::btSDP) {
    int l2 = bs.blockNumber[l-1];
    currentPt.zMat.setElement_SDP(l2,i-1,j-1,value);
  } else if (bs.blockType[l-1] == BlockStruct::btSOCP) {
    rError("io:: current version does not support SOCP");
    int l2 = bs.blockNumber[l-1];
    currentPt.zMat.setElement_SOCP(l2,i-1,j-1,value);
  } else if (bs.blockType[l-1] == BlockStruct::btLP) {
    currentPt.zMat.setElement_LP(bs.blockNumber[l-1]+i-1,value);
  }
}

void SDPA::inputInitYMat(int l, int i, int j, double value)
{
  if (l > nBlock || l <= 0) {
    rError ("l exceeds nBlock or "
	    "l is less than or equal to zero :: nBlock= "
	    << nBlock << " : l= " << l
	    << " : i= " << i << " : j= " << j);
  }
  int dim = bs.blockStruct[l-1];
  if (i > dim || i <= 0) {
    rError ("i exceeds dimension of the block or "
	    "i is less than or equal to zero :: dim= "
	    << dim << " : l= " << l
	    << " : i= " << i << " : j= " << j);
  }
  if (j > dim || j <= 0) {
    rError ("j exceeds dimension of the block or "
	    "j is less than or equal to zero :: dim= "
	    << dim << " : l= " << l
	    << " : i= " << i << " : j= " << j);
  }
  if (bs.blockType[l-1] == BlockStruct::btLP) {
    if (i!=j) {
      rError("i should be j in LP block"
	     " : l= " << l
	     << " : i= " << i << " : j= " << j);
    }
  }

  if (bs.blockType[l-1] == BlockStruct::btSDP) {
    int l2 = bs.blockNumber[l-1];
    currentPt.xMat.setElement_SDP(l2,i-1,j-1,value);
  } else if (bs.blockType[l-1] == BlockStruct::btSOCP) {
    rError("io:: current version does not support SOCP");
    int l2 = bs.blockNumber[l-1];
    currentPt.xMat.setElement_SOCP(l2,i-1,j-1,value);
  } else if (bs.blockType[l-1] == BlockStruct::btLP) {
    currentPt.xMat.setElement_LP(bs.blockNumber[l-1]+i-1,value);
  }
}


void SDPA::initializeUpperTriangleSpace()
{
  bs.makeInternalStructure();
  NewArray(NonZeroElements,vector<IndexLIJv*>,m+1);
  currentPt.initialize(m, bs, param.lambdaStar, com);
  inputData.initialize(bs);
  inputData.initialize_bVec(m);
}


void SDPA::printNonZeroElements(FILE* fp)
{
  for (int k=0; k<=m; ++k) {
    int size = NonZeroElements[k].size();
    for (int index = 0; index<size; ++index) {
      IndexLIJv* a = NonZeroElements[k][index];
      int l        = a[0].l;
      int i        = a[0].i;
      int j        = a[0].j;
      double value = a[0].value;
      fprintf(fp,"%d, %d, %d, %d, ",k,l,i,j);
      fprintf(fp,param.infPrint,value);
      fprintf(fp,"\n");
    }
  }
}

void SDPA::sortNonZeroElements()
{
  for (int k=0; k<=m; ++k) {
    sort(NonZeroElements[k].begin(), NonZeroElements[k].end(),
	 IndexLIJv::compare);
  }
}

void SDPA::checkNonZeroElements()
{
  TimeStart(FILE_CHECK_START1);
  for (int k=0; k<=m; ++k) {
    int size = NonZeroElements[k].size();
    for (int index = 0; index<size-1; ++index) {
      IndexLIJv* a = NonZeroElements[k][index];
      IndexLIJv* b = NonZeroElements[k][index+1];
      if (a[0].l == b[0].l && a[0].i == b[0].i
	  && a[0].j == b[0].j) {
	int l        = a[0].l;
	int i        = a[0].i;
	int j        = a[0].j;
	rError("Twice input to the same index. "
	       ": k = " << k << ": l = " << l
	       << ": i = " << i << ": j = " << j);
      }
    }
  }
  TimeEnd(FILE_CHECK_END1);
  com.FileChange += TimeCal(FILE_CHECK_START1,
			    FILE_CHECK_END1);
  com.TotalTime += TimeCal(FILE_CHECK_START1,
			    FILE_CHECK_END1);
}

void SDPA::setNonZeroBlockStruct()
{
  // almost equivalent to IO::setBlockStruct
  NewArray(inputData.A,SparseLinearSpace,m);
  
  int  SDP_sp_nBlock;
  int* SDP_sp_NonZeroNumber;
  int* SDP_sp_index;
  int* SDP_sp_blockStruct;
  int  SOCP_sp_nBlock;
  int* SOCP_sp_NonZeroNumber;
  int* SOCP_sp_index;
  int* SOCP_sp_blockStruct;
  int  LP_sp_nBlock;
  int* LP_sp_index;

  NewArray(SDP_sp_index,         int,bs.SDP_nBlock);
  NewArray(SDP_sp_blockStruct,   int,bs.SDP_nBlock);
  NewArray(SDP_sp_NonZeroNumber, int,bs.SDP_nBlock);
  NewArray(SOCP_sp_index,        int,bs.SOCP_nBlock);
  NewArray(SOCP_sp_blockStruct,  int,bs.SOCP_nBlock);
  NewArray(SOCP_sp_NonZeroNumber,int,bs.SOCP_nBlock);
  NewArray(LP_sp_index,          int,bs.LP_nBlock);

  for (int k=0; k<=m; ++k) {
    SDP_sp_nBlock  = 0;
    SOCP_sp_nBlock = 0;
    LP_sp_nBlock   = 0;

    int previous_l = -1;
    int size = NonZeroElements[k].size();
    // rMessage("NonZeroElements[" << k << "].size = " << size);
    for (int index = 0; index < size; ++index) {
      IndexLIJv* a = NonZeroElements[k][index];
      int l = a[0].l;
      int i = a[0].i;
      if (bs.blockType[l-1] == BlockStruct::btSDP) {
	if (l!=previous_l) {
	  int l2 = bs.blockNumber[l-1];
	  SDP_sp_index        [SDP_sp_nBlock] = l2;
	  SDP_sp_blockStruct  [SDP_sp_nBlock] = bs.SDP_blockStruct[l2];
	  SDP_sp_NonZeroNumber[SDP_sp_nBlock] = 1;
	  previous_l = l;
	  SDP_sp_nBlock++;
	}
	else {
	  SDP_sp_NonZeroNumber[SDP_sp_nBlock-1]++;
	}
	#if 0
	rMessage("k = " << k <<
		 " : SDP_sp_NonZeroNumber["<< (SDP_sp_nBlock-1)
		 <<"] = " << SDP_sp_NonZeroNumber[SDP_sp_nBlock-1]);
	#endif

      }
      else if (bs.blockType[l-1] == BlockStruct::btSOCP) {
	rError("current version does not support SOCP");
      }
      else if (bs.blockType[l-1] == BlockStruct::btLP) {
	int l2 = bs.blockNumber[l-1];
	LP_sp_index[LP_sp_nBlock] = l2 + i - 1;
	previous_l = l;
	LP_sp_nBlock++;
      }
    }
    
    if (k==0) { // coefficient matrix
      inputData.C.initialize(SDP_sp_nBlock,
			     SDP_sp_index,
			     SDP_sp_blockStruct,
			     SDP_sp_NonZeroNumber,
			     SOCP_sp_nBlock,
			     SOCP_sp_index,
			     SOCP_sp_blockStruct,
			     SOCP_sp_NonZeroNumber,
			     LP_sp_nBlock,
			     LP_sp_index);
    }
    else { // input data matrix
      inputData.A[k-1].initialize(SDP_sp_nBlock,
			     SDP_sp_index,
			     SDP_sp_blockStruct,
			     SDP_sp_NonZeroNumber,
			     SOCP_sp_nBlock,
			     SOCP_sp_index,
			     SOCP_sp_blockStruct,
			     SOCP_sp_NonZeroNumber,
			     LP_sp_nBlock,
			     LP_sp_index);
    }
  }
  DeleteArray(SDP_sp_index);
  DeleteArray(SDP_sp_blockStruct);
  DeleteArray(SDP_sp_NonZeroNumber);
  DeleteArray(SOCP_sp_index);
  DeleteArray(SOCP_sp_blockStruct);
  DeleteArray(SOCP_sp_NonZeroNumber);
  DeleteArray(LP_sp_index);
}

void SDPA::setNonZeroElements()
{
  // almost equivalent to IO::setElement
  for (int k=0; k<=m; ++k) {
    int size = NonZeroElements[k].size();
    for (int index = 0; index < size; ++index) {
      IndexLIJv* a = NonZeroElements[k][index];
      int l        = a[0].l;
      int i        = a[0].i;
      int j        = a[0].j;
      double value = a[0].value;
      
      if (bs.blockType[l-1] == BlockStruct::btSDP) {
        int l2 = bs.blockNumber[l-1];
	#if 0
	rMessage("k = " << k << ": l = " << l
		 << " : l2 = " << l2 
		 << " : i = " << i << " : j = " << j);
	#endif
        if (k==0) {
          inputData.C.setElement_SDP(l2,i-1,j-1,-value);
        } else {
          inputData.A[k-1].setElement_SDP(l2,i-1,j-1,value);
        }
      }
      else if (bs.blockType[l-1] == BlockStruct::btSOCP) {
        rError("io:: current version does not support SOCP");
        int l2 = bs.blockNumber[l-1];
        if (k==0) {
          inputData.C.setElement_SOCP(l2,i-1,j-1,-value);
        } else {
          inputData.A[k-1].setElement_SOCP(l2,i-1,j-1,value);
        }
      }
      else if (bs.blockType[l-1] == BlockStruct::btLP) {
        if (i != j){
          rError("io:: LP part  3rd element != 4th element\n"
                 "column should be same as row in LP part.");
        }
        if (k==0) {
          inputData.C.setElement_LP(bs.blockNumber[l-1]+i-1,-value);
        } else {
          inputData.A[k-1].setElement_LP(bs.blockNumber[l-1]+i-1,value);
        }
      } else {
        rError("io::read not valid blockType");
      }
    }
  }

}

void SDPA::initializeUpperTriangle(bool checkTwiceInput)
{
  sortNonZeroElements();
  if (checkTwiceInput) {
    checkNonZeroElements();
  }
  // printNonZeroElements();
  setNonZeroBlockStruct();
  setNonZeroElements();
  for (int k=0; k<=m; ++k) {
    int size = NonZeroElements[k].size();
    for (int index = 0; index < size; ++index) {
      DeleteArray(NonZeroElements[k][index]);
    }
  }
  DeleteArray(NonZeroElements);
}
  
double* SDPA::getResultXVec()
{
  return currentPt.yVec.ele;
}

double* SDPA::getResultXMat(int l)
{
  if (l > nBlock || l <= 0) {
      rError ("l exceeds nBlock or "
	      "l is less than or equal to zero :: nBlock= "
	      << nBlock << " : l= " << l);
  }
  if (bs.blockType[l-1] == BlockStruct::btSDP) {
    int l2 = bs.blockNumber[l-1];
    return currentPt.zMat.SDP_block[l2].de_ele;
  }
  else if (bs.blockType[l-1] == BlockStruct::btSOCP) {
    rError("io:: current version does not support SOCP");
    int l2 = bs.blockNumber[l-1];
    return currentPt.zMat.SOCP_block[l2].de_ele;
  }
  else if (bs.blockType[l-1] == BlockStruct::btLP) {
    int start = bs.blockNumber[l-1];
    return &currentPt.zMat.LP_block[start];
  }
  return NULL;
}

double* SDPA::getResultYMat(int l)
{
  if (l > nBlock || l <= 0) {
      rError ("l exceeds nBlock or "
	      "l is less than or equal to zero :: nBlock= "
	      << nBlock << " : l= " << l);
  }
  if (bs.blockType[l-1] == BlockStruct::btSDP) {
    int l2 = bs.blockNumber[l-1];
    return currentPt.xMat.SDP_block[l2].de_ele;
  }
  else if (bs.blockType[l-1] == BlockStruct::btSOCP) {
    rError("io:: current version does not support SOCP");
    int l2 = bs.blockNumber[l-1];
    return currentPt.xMat.SOCP_block[l2].de_ele;
  }
  else if (bs.blockType[l-1] == BlockStruct::btLP) {
    int start = bs.blockNumber[l-1];
    return &currentPt.xMat.LP_block[start];
  }
  return NULL;
}

double SDPA::getPrimalObj()
{
  // Note reverse primal-dual
  return -solveInfo.objValDual;
}

double SDPA::getDualObj()
{
  // Note reverse primal-dual
  return -solveInfo.objValPrimal;
}

double SDPA::getPrimalError()
{
  // Note reverse primal-dual
  return currentRes.normDualMat;
}

double SDPA::getDualError()
{
  // Note reverse primal-dual
  return currentRes.normPrimalVec;
}

double SDPA::getDigits()
{
  double mean =  (fabs(solveInfo.objValPrimal)
		  + fabs(solveInfo.objValDual)) / 2.0;
  double PDgap = getDualityGap();
  double digits = -log10(fabs(PDgap/mean));
  return digits;
}

int SDPA::getIteration()
{
  return pIteration;
}

double SDPA::getMu()
{
  return mu.current;
}

double SDPA::getDualityGap()
{
  double PDgap = fabs(solveInfo.objValPrimal
		      - solveInfo.objValDual);
  return PDgap;
}

SDPA::PhaseType SDPA::getPhaseValue()
{
  // Note reverse primal-dual
  switch (phase.value) {
  case SolveInfo::noINFO    : return noINFO    ; break;
  case SolveInfo::pFEAS     : return pFEAS     ; break;
  case SolveInfo::dFEAS     : return dFEAS     ; break;
  case SolveInfo::pdFEAS    : return pdFEAS    ; break;
  case SolveInfo::pdINF     : return pdINF     ; break;
  case SolveInfo::pFEAS_dINF: return pINF_dFEAS; break;
  case SolveInfo::pINF_dFEAS: return pFEAS_dINF; break;
  case SolveInfo::pdOPT     : return pdOPT     ; break;
  case SolveInfo::pUNBD     : return dUNBD     ; break;
  case SolveInfo::dUNBD     : return pUNBD     ; break;
  default: break;
  }
  return noINFO;
}

void SDPA::getPhaseString(char* str)
{
  switch (phase.value) {
  case SolveInfo::noINFO    : strcpy(str,(char *)"noINFO    "); break;
  case SolveInfo::pFEAS     : strcpy(str,(char *)"pFEAS     "); break;
  case SolveInfo::dFEAS     : strcpy(str,(char *)"dFEAS     "); break;
  case SolveInfo::pdFEAS    : strcpy(str,(char *)"pdFEAS    "); break;
  case SolveInfo::pdINF     : strcpy(str,(char *)"pdINF     "); break;
  case SolveInfo::pFEAS_dINF: strcpy(str,(char *)"pFEAS_dINF"); break;
  case SolveInfo::pINF_dFEAS: strcpy(str,(char *)"pINF_dFEAS"); break;
  case SolveInfo::pdOPT     : strcpy(str,(char *)"pdOPT     "); break;
  case SolveInfo::pUNBD     : strcpy(str,(char *)"pUNBD     "); break;
  case SolveInfo::dUNBD     : strcpy(str,(char *)"dUNBD     "); break;
  default:
    strcpy(str,(char *)"phase error");
    break;
  }
  return;
}

double SDPA::getSolveTime()
{
  return com.TotalTime;
}

int SDPA::getConstraintNumber()
{
  return m;
}

int SDPA::getBlockNumber()
{
  return nBlock;
}

int SDPA::getBlockSize(int l)
{
  if (l<=0 || l>nBlock) {
    rMessage("out of range : getBlockSize "
	     ": l = " << l
	     << " should be between 1 and nBlock " << nBlock);
  }
  return bs.blockStruct[l-1];
}

SDPA::ConeType SDPA::getBlockType(int l)
{
  if (l<=0 || l>nBlock) {
    rMessage("out of range : getBlockSize "
	     ": l = " << l
	     << " should be between 1 and nBlock " << nBlock);
  }
  switch (bs.blockType[l-1]) {
  case BlockStruct::btSDP  : return SDPA::SDP ;
  case BlockStruct::btSOCP : return SDPA::SOCP;
  case BlockStruct::btLP   : return SDPA::LP  ;
  }
  rError("Type Error in getBlockType ");
  return SDPA::SDP; // dummy return
}

void SDPA::getDimacsError(double* DimacsError)
{
  IO::computeDimacs(DimacsError, solveInfo, currentRes,
		    currentPt, inputData, work);
}

void SDPA::printDimacsError(double* DimacsError, char* printFormat,
			    FILE* fpout)
{
  IO::printDimacs(DimacsError,printFormat,fpout);
}

  
void SDPA::printResultXVec(FILE* fp)
{
  // Note reverse primal-dual
  currentPt.yVec.display(fp,1.0,param.xPrint);
}

void SDPA::printResultXMat(FILE* fp)
{
  // Note reverse primal-dual
  currentPt.zMat.displaySolution(bs,fp,param.XPrint);
}

void SDPA::printResultYMat(FILE* fp)
{
  // Note reverse primal-dual
  currentPt.xMat.displaySolution(bs,fp,param.YPrint);
}

void SDPA::printComputationTime(FILE* fp)
{
  com.display(fp);
}

void SDPA::printParameters(FILE* fp)
{
  param.display(fp);
}

void SDPA::printSDPAVersion(FILE* fp)
{
  if (fp) {
    fprintf(fp,"%s\n",(char*)sdpa_right);
  }
}

void SDPA::readInput(char* filename, FILE* fpout, SparseType type)
{
  if (type == AUTO) {
    int len = strlen(filename);
    if (filename[len-1] =='s' && filename[len-2] == '-') {
      type = SPARSE;
    }
    else {
      type = DENSE;
    }
  }
  bool isDataSparse = true;
  if (type == DENSE) {
    isDataSparse = false;
  }
  TimeStart(FILE_READ_START1);
  FILE* fpinput = NULL;
  if ((fpinput = fopen(filename,"r")) == NULL) {
    rError("Cannot Open Data File " << filename);
  }
  if (fpout){ 
    fprintf(fpout,"data   is %s ", filename);
    if (isDataSparse) {
      fprintf(fpout," : sparse\n");
    }
    else {
      fprintf(fpout," : dense\n");
    }
  }
  
  char titleAndComment[LengthOfBuffer];
  IO::read(fpinput,fpout,m,titleAndComment);
  IO::read(fpinput,nBlock);
  bs.initialize(nBlock);
  IO::read(fpinput,bs);
  bs.makeInternalStructure();
  inputData.initialize(bs);
  IO::read(fpinput,m,bs,inputData,isDataSparse);
  // inputData.initialize_index();
  fclose(fpinput);
  currentPt.initialize(m, bs, param.lambdaStar,com);
  TimeEnd(FILE_READ_END1);
  com.FileRead += TimeCal(FILE_READ_START1,
			  FILE_READ_END1);
  com.TotalTime += TimeCal(FILE_READ_START1,
			  FILE_READ_END1);
  return;
}

void SDPA::readInit(char* filename,  FILE* fpout, SparseType type)
{
  TimeEnd(FILE_READ_START2);
  if (type == AUTO) {
    int len = strlen(filename);
    if (filename[len-1] =='s' && filename[len-2] == '-') {
      type = SPARSE;
    }
    else {
      type = DENSE;
    }
  }
  bool isInitSparse = true;
  if (type == DENSE) {
    isInitSparse = false;
  }
  FILE* fpinit = NULL;
  if ((fpinit = fopen(filename,"r")) == NULL) {
    rError("Cannot Open Init File " << filename);
  }
  if (fpout) {
    fprintf(fpout,"initial is %s ", filename);
    if (isInitSparse) {
      fprintf(fpout," : sparse\n");
    }
    else {
      fprintf(fpout," : dense\n");
    }
  }
  IO::read(fpinit,currentPt.xMat,currentPt.yVec,currentPt.zMat,
	   bs, isInitSparse);
#if 0
  rMessage("intial X = ");
  currentPt.xMat.display();
  rMessage("intial Z = ");
  currentPt.zMat.display();
#endif
  fclose(fpinit);

  TimeEnd(FILE_READ_END2);
  com.FileRead += TimeCal(FILE_READ_START2,
			   FILE_READ_END2);
  com.TotalTime += TimeCal(FILE_READ_START2,
			   FILE_READ_END2);
  return;
}

void SDPA::readParameter(char* filename, FILE* fpout)
{
  FILE* fp = NULL;
  if ((fp=fopen(filename,"r"))==NULL) {
    rError("Cannot Open parameter File " << filename);
  }
  if (fpout) {
    fprintf(fpout,"param  is %s \n", filename);
  }
  param.readFile(fp);
  fclose(fp);
  return;
}

void SDPA::writeSparseLinearSpace(FILE* fp, char* printFormat,
				  sdpa::SparseLinearSpace& A, int k)
{
  // bs.display();
  int  SDP_sp_nBlock  = A.SDP_sp_nBlock;
  // int  SOCP_sp_nBlock = A.SOCP_sp_nBlock;
  int  LP_sp_nBlock   = A.LP_sp_nBlock;
  int* SDP_sp_index   = A.SDP_sp_index;
  // int* SOCP_sp_index  = A.SOCP_sp_index;
  int* LP_sp_index    = A.LP_sp_index;

  for (int l=0; l<SDP_sp_nBlock; ++l) {
    int l2 = SDP_sp_index[l];
    int original_l = 0;
    for (int l3=0,countSDP=0; l3<bs.nBlock; ++l3) {
      if (bs.blockType[l3] == BlockStruct::btSDP) {
	if (countSDP == l2) {
	  original_l = l3;
	  break;
	}
	countSDP++;
      }
    }
    // rMessage("l2 = " << l2);
    SparseMatrix& Al = A.SDP_sp_block[l];
    if (Al.type == SparseMatrix::SPARSE) {
      // rMessage("sparse: Al.NonZeroCount = " << Al.NonZeroCount);
      for (int index=0; index<Al.NonZeroCount; ++index) {
	int i, j;
	double value;
	if (Al.DataStruct == SparseMatrix::DSarrays) {
	  i     = Al.row_index[index];
	  j     = Al.column_index[index];
	  value = Al.sp_ele[index];
	}
	else {
	  i     = Al.DataS[index].vRow;
	  j     = Al.DataS[index].vCol;
	  value = Al.DataS[index].vEle;
	}
	if (value == 0.0) {
	  continue;
	}
	if (k==0) {
	  value = -value;
	}
	fprintf(fp,"%d %d %d %d ", k, original_l+1,i+1,j+1);
	fprintf(fp,printFormat,value);
	fprintf(fp,"\n");
      }
    }
    else {
      for (int i=0; i<Al.nRow; ++i) {
	for (int j=i; j<Al.nCol; ++j) {
	  double value = Al.de_ele[i+Al.nRow*j];
	  if (value == 0.0) {
	    continue;
	  }
	  if (k==0) {
	    value = -value;
	  }
	  fprintf(fp,"%d %d %d %d ", k, original_l+1,i+1,j+1);
	  fprintf(fp,printFormat,value);
	  fprintf(fp,"\n");
	}
      }
    }
  }

  // for SOCP
  // rError("io:: current version does not support SOCP");

  // LP is always SPARSE format
  for (int l=0; l<LP_sp_nBlock; ++l) {
    double value = A.LP_sp_block[l];
    if (k==0) {
      value = -value;
    }
    int ik = LP_sp_index[l];
    int l2 = 0;
    for (l2=0; l2<nBlock; ++l2) {
      if (bs.blockType[l2] == BlockStruct::btLP
	  && bs.blockNumber[l2] <= ik
	  && ik < bs.blockNumber[l2] + bs.blockStruct[l2]) {
	break;
      }
    }
    int i = ik - bs.blockNumber[l2];
    fprintf(fp, "%d %d %d %d ", k, l2+1, i+1, i+1);
    fprintf(fp,printFormat,value);
    fprintf(fp,"\n");
  }
}

void SDPA::writeInputSparse(char* filename, char* printFormat)
{
  FILE* fp = NULL;
  if ((fp=fopen(filename,"w"))==NULL) {
    rError("Cannot Open Data File to Write" << filename);
  }
  fprintf(fp,"%d\n",m);
  fprintf(fp,"%d\n",nBlock);
  for (int l=0; l<nBlock-1; ++l) {
    if (bs.blockType[l] == BlockStruct::btSDP) {
      fprintf(fp,"%d,",bs.blockStruct[l]);
    }
    else if (bs.blockType[l] == BlockStruct::btSOCP) {
      rError("io:: current version does not support SOCP");
      fprintf(fp,"%d,",bs.blockStruct[l]);
    }
    else if (bs.blockType[l] == BlockStruct::btLP) {
      fprintf(fp,"%d,",-bs.blockStruct[l]);
    }
  }
  if (bs.blockType[nBlock-1] == BlockStruct::btSDP) {
    fprintf(fp,"%d\n",bs.blockStruct[nBlock-1]);
  }
  else if (bs.blockType[nBlock-1] == BlockStruct::btSOCP) {
    rError("io:: current version does not support SOCP");
    fprintf(fp,"%d\n",bs.blockStruct[nBlock-1]);
  }
  else if (bs.blockType[nBlock-1] == BlockStruct::btLP) {
    fprintf(fp,"%d\n",-bs.blockStruct[nBlock-1]);
  }

  if (strcmp(printFormat,NO_P_FORMAT) == 0) {
    fprintf(fp,"%s\n",NO_P_FORMAT);
  }
  else {
    for (int k=0; k<m; ++k) {
      fprintf(fp,printFormat,inputData.b.ele[k]);
      fprintf(fp," ");
    }
    fprintf(fp,"\n");
  
    writeSparseLinearSpace(fp, printFormat,inputData.C,0);
    for (int k=0; k<m; ++k) {
      // rMessage("k=" << k);
      writeSparseLinearSpace(fp, printFormat,inputData.A[k],k+1);
    }
  }
  fclose(fp);
  return;
}

void SDPA::writeDenseLinearSpace(FILE* fp, char* printFormat,
				 sdpa::DenseLinearSpace& X,int k)
{
  int  SDP_nBlock  = X.SDP_nBlock;
  // int  SOCP_nBlock = X.SOCP_nBlock;
  int  LP_nBlock   = X.LP_nBlock;
  
  for (int l=0; l<SDP_nBlock; ++l) {
    int l2 = 0;
    for (l2 = 0; l2 < nBlock; ++l2) {
      if (bs.blockNumber[l2] == l) {
	break;
      }
    }
    DenseMatrix& Xl = X.SDP_block[l];
    for (int i=0; i<Xl.nRow; ++i) {
      for (int j=i; j<Xl.nCol; ++j) {
	double value = Xl.de_ele[i+Xl.nRow*j];
	if (value == 0.0) {
	  continue;
	}
	fprintf(fp,"%d %d %d %d ", k,l2+1,i+1,j+1);
	fprintf(fp,printFormat,value);
	fprintf(fp,"\n");
      }
    }
  }
  
  // for SOCP
  // rError("io:: current version does not support SOCP");

  // LP 
  int l2 = 0;
  for (int l=0; l<LP_nBlock; ++l) {
    double value = X.LP_block[l];
    for (l2=0; l2<nBlock; ++l2) {
      if (bs.blockType[l2] == BlockStruct::btLP
	  && bs.blockNumber[l2] <= l
	  && l < bs.blockNumber[l2] + bs.blockStruct[l2]) {
	break;
      }
    }
    int i = l - bs.blockNumber[l2];
    fprintf(fp, "%d %d %d %d ", k, l2+1, i+1, i+1);
    fprintf(fp,printFormat,value);
    fprintf(fp,"\n");
  }
}

void SDPA::writeInitSparse(char* filename, char* printFormat)
{
  FILE* fp = NULL;
  if ((fp=fopen(filename,"w"))==NULL) {
    rError("Cannot Open Init File to Write" << filename);
  }
  if (strcmp(printFormat,NO_P_FORMAT) == 0) {
    fprintf(fp,"%s\n",NO_P_FORMAT);
    fclose(fp);
    return;
  }
  // Note reverse primal-dual
  for (int k=0; k<m; ++k) {
    fprintf(fp,printFormat,-currentPt.yVec.ele[k]);
    fprintf(fp," ");
  }
  fprintf(fp,"\n");
  writeDenseLinearSpace(fp,printFormat,currentPt.zMat,1);
  writeDenseLinearSpace(fp,printFormat,currentPt.xMat,2);
  fclose(fp);
  return;
}
  
void SDPA::terminate()
{
  bs.terminate();
  inputData.terminate();
  chordal.terminate();
  newton.terminate();
  currentPt.terminate();
  work.terminate();
  initPt_xMat.terminate();
  initPt_zMat.terminate();
  initRes.terminate();
  currentRes.terminate();
  alpha.terminate();
}


void SDPA::copyCurrentToInit()
{
  // initPt_xMat.copyFrom(currentPt.xMat);
  // initPt_zMat.copyFrom(currentPt.zMat);
  // Note reverse primal-dual
  Lal::let(currentPt.yVec,'=',currentPt.yVec,'*',&DMONE);
  return;
}

void SDPA::setKappa(double KAPPA)
{
  this->KAPPA = KAPPA;
}

