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

#include "sdpa_jordan.h"
#include "sdpa_dataset.h"
#include "sdpa_linear.h"
#include "sdpa_algebra.h"

namespace sdpa {

double Jal::trace(DenseLinearSpace& aMat)
{
  double ret = 0.0;

  for (int l=0; l<aMat.SDP_nBlock; ++l) {
    double* target = aMat.SDP_block[l].de_ele;
    int size = aMat.SDP_block[l].nRow;
    for (int j=0; j<size; ++j) {
	  ret += target[j+size*j];
	}
  }

  for (int l=0; l<aMat.SOCP_nBlock; ++l) {
    rError("dataset:: current version do not support SOCP");
  }

  for (int l=0; l<aMat.LP_nBlock; ++l) {
    ret += aMat.LP_block[l];
  }

  return ret;
}

// calculate the minimum eigen value of lMat*xMat*(lMat^T)
// lMat is lower triangular, xMat is symmetric
// block size > 64   : Lanczos method
// block size <= 64  : QR method
// QR method: workVec is temporary space and needs
//            3*xMat.nRow-1 length memory.
double Jal::getMinEigen(DenseLinearSpace& lMat,
		       DenseLinearSpace& xMat,
		       WorkVariables& work)
{
  double min = 1.0E50;
  double value;

  // for SDP
  for (int l=0; l<xMat.SDP_nBlock; ++l) {
    if (xMat.SDP_block[l].nRow > 64){ // use Lanczos method
      value = Lal::getMinEigen(lMat.SDP_block[l],
			       xMat.SDP_block[l],
			       work.DLS1.SDP_block[l],
			       work.SDP_BV1.ele[l],
			       work.SDP_BV2.ele[l],
			       work.SDP_BV3.ele[l],
			       work.SDP_BV4.ele[l],
			       work.SDP_BV5.ele[l],
			       work.SDP_BV6.ele[l],
			       work.SDP_BV7.ele[l],
			       work.SDP_BV8.ele[l],
			       work.SDP_BV9.ele[l],
			       work.SDP2_BV1.ele[l]);
    } else { // use QR method
      Lal::let(work.DLS2.SDP_block[l],
	       '=',xMat.SDP_block[l],'T',lMat.SDP_block[l]);
      Lal::let(work.DLS1.SDP_block[l],
	       '=',lMat.SDP_block[l],'*',work.DLS2.SDP_block[l]);
      Lal::getMinEigenValue(work.DLS1.SDP_block[l],
			    work.SDP_BV1.ele[l],work.SDP2_BV1.ele[l]);
      value = work.SDP_BV1.ele[l].ele[0];
    }
    if (value < min) {
      min = value;
    }
  } // end of 'for (int l)'

  // for SOCP
  for (int l=0; l<xMat.SOCP_nBlock; ++l) {
    rError("getMinEigen:: current version does not support SOCP");
  }

// for  LP
  for (int l=0; l<xMat.LP_nBlock; ++l) {
    value = xMat.LP_block[l] * lMat.LP_block[l] * lMat.LP_block[l];
    if (value < min) {
      min = value;
    }
  } // end of 'for (int l)'

  return min;
};

  // calculate the minimum eigen value of xMat by QR method.
double Jal::getMinEigen(DenseLinearSpace& xMat,
                        WorkVariables& work)
{
  double min = 1.0E50;
  double value;

  work.DLS1.copyFrom(xMat);

  // for SDP
  for (int l=0; l<xMat.SDP_nBlock; ++l) {
    Lal::getMinEigenValue(work.DLS1.SDP_block[l],
			  work.SDP_BV1.ele[l],work.SDP2_BV1.ele[l]);
    value = work.SDP_BV1.ele[l].ele[0];

    if (value < min) {
      min = value;
    }
  } // end of 'for (int l)'

  // for SOCP
  for (int l=0; l<xMat.SOCP_nBlock; ++l) {
    rError("getMinEigen:: current version does not support SOCP");
  }

// for  LP
  for (int l=0; l<xMat.LP_nBlock; ++l) {
    value = xMat.LP_block[l];
    if (value < min) {
      min = value;
    }
  } // end of 'for (int l)'

  return min;
};


bool Jal::getInvChol(DenseLinearSpace& invCholMat,
		     DenseLinearSpace& aMat,
		     DenseLinearSpace& workMat)
{
  // for SDP
  if (workMat.SDP_nBlock!=aMat.SDP_nBlock
      || invCholMat.SDP_nBlock!=aMat.SDP_nBlock) {
    rError("getInvChol:: different memory size");
  }
  for (int l=0; l<aMat.SDP_nBlock; ++l) {
    if (Lal::getCholesky(workMat.SDP_block[l],
			 aMat.SDP_block[l]) == false) {
      return false;
    }
    Lal::getInvLowTriangularMatrix(invCholMat.SDP_block[l],
				   workMat.SDP_block[l]);
  }

  // for SOCP
  for (int l=0; l<aMat.SOCP_nBlock; ++l) {
    rError("no support for SOCP");
  }

  // fo LP
  if (invCholMat.LP_nBlock!=aMat.LP_nBlock) {
    rError("getInvChol:: different memory size");
  }
  for (int l=0; l<aMat.LP_nBlock; ++l) {
    if (aMat.LP_block[l] < 0.0){
      return false;
    }
    invCholMat.LP_block[l] = 1.0 / sqrt(aMat.LP_block[l]);
  }

  return SDPA_SUCCESS;
}

bool Jal::getInvCholAndInv(DenseLinearSpace& invCholMat,
			   DenseLinearSpace& inverseMat,
			   DenseLinearSpace& aMat,
			   DenseLinearSpace& workMat)
{
  if (getInvChol(invCholMat, aMat, workMat)== false) {
    return SDPA_FAILURE;
  }

  for (int l=0; l<aMat.SDP_nBlock; ++l) {
    inverseMat.SDP_block[l].copyFrom(invCholMat.SDP_block[l]);
    dtrmm_fc ((char *)"Left",(char *)"Lower",
	       (char *)"Transpose",(char *)"NonUnitDiag",
	       &invCholMat.SDP_block[l].nRow,
	       &invCholMat.SDP_block[l].nCol,
	       &DONE,
	       invCholMat.SDP_block[l].de_ele,
	       &invCholMat.SDP_block[l].nRow,
	       inverseMat.SDP_block[l].de_ele,
	       &inverseMat.SDP_block[l].nRow,
	       strlen("Left"),strlen("Lower"),
	       strlen("Transpose"),strlen("NonUnitDiag") );
  }
  for (int l=0; l<aMat.SOCP_nBlock; ++l) {
    rError("rNewton:: we don't make this ruoutin");
    
  }
  for (int l=0; l<aMat.LP_nBlock; ++l) {
    inverseMat.LP_block[l] = 1.0 / aMat.LP_block[l];
  }
  return SDPA_SUCCESS;
}


bool Jal::multiply(DenseLinearSpace& retMat,
		   DenseLinearSpace& aMat,
		   DenseLinearSpace& bMat,
		   double* scalar)
{
  bool total_judge = SDPA_SUCCESS;

  // for SDP
  if (retMat.SDP_nBlock!=aMat.SDP_nBlock 
      || retMat.SDP_nBlock!=bMat.SDP_nBlock) {
    rError("multiply:: different nBlock size");
  }
  for (int l=0; l<retMat.SDP_nBlock; ++l) {
    bool judge = Lal::multiply(retMat.SDP_block[l],aMat.SDP_block[l],
			       bMat.SDP_block[l],scalar);
    if (judge == SDPA_FAILURE) {
      total_judge = SDPA_FAILURE;
    }
  }

  // for SOCP
#if 0
  if (retMat.SOCP_nBlock!=aMat.SOCP_nBlock 
      || retMat.SOCP_nBlock!=bMat.SOCP_nBlock) {
    rError("multiply:: different nBlock size");
  }
  for (int l=0; l<retMat.SOCP_nBlock; ++l) {
    bool judge = Lal::multiply(retMat.SOCP_block[l],aMat.SOCP_block[l],
			       bMat.SOCP_block[l],scalar);
    if (judge == SDPA_FAILURE) {
      total_judge = SDPA_FAILURE;
    }
  }
#endif

  // for LP
  if (retMat.LP_nBlock!=aMat.LP_nBlock 
      || retMat.LP_nBlock!=bMat.LP_nBlock) {
    rError("multiply:: different nBlock size");
  }
  if (scalar == NULL) {
    for (int l=0; l<retMat.LP_nBlock; ++l) {
      retMat.LP_block[l] = aMat.LP_block[l] * bMat.LP_block[l];
    }
  }
  else {
    for (int l=0; l<retMat.LP_nBlock; ++l) {
      retMat.LP_block[l] = aMat.LP_block[l]
	* bMat.LP_block[l] * (*scalar);
    }
  }
  return total_judge;
}

#if 0
// CAUTION!!! We don't initialize retMat to zero matrix for efficiently.
bool Jal::multiply(DenseLinearSpace& retMat,
		   SparseLinearSpace& aMat,
		   DenseLinearSpace& bMat,
		   double* scalar)
{
  bool total_judge = SDPA_SUCCESS;

  // for SDP
  for (int l=0; l<aMat.SDP_sp_nBlock; ++l) {
    int index = aMat.SDP_sp_index[l];
    bool judge = Lal::multiply(retMat.SDP_block[index],
			       aMat.SDP_sp_block[l],
			       bMat.SDP_block[index],scalar);
    if (judge == SDPA_FAILURE) {
      total_judge = SDPA_FAILURE;
    }
  }

  // for SOCP
  for (int l=0; l<aMat.SOCP_sp_nBlock; ++l) {
    int index = aMat.SOCP_sp_index[l];
    bool judge = Lal::multiply(retMat.SOCP_block[index],
			       aMat.SOCP_sp_block[l],
			       bMat.SOCP_block[index],scalar);
    if (judge == SDPA_FAILURE) {
      total_judge = SDPA_FAILURE;
    }
  }

  // for LP
  if (scalar == NULL) {
    for (int l=0; l<aMat.LP_sp_nBlock; ++l) {
      int index = aMat.LP_sp_index[l];
      retMat.LP_block[index] = 
	aMat.LP_sp_block[l] *  bMat.LP_block[index];
    }
  }
  else {
    for (int l=0; l<aMat.LP_sp_nBlock; ++l) {
      int index = aMat.LP_sp_index[l];
      retMat.LP_block[index] = 
	aMat.LP_sp_block[l] *  bMat.LP_block[index] * (*scalar);
    }
  }

  return total_judge;
}

// CAUTION!!! We don't initialize retMat to zero matrix for efficiently.
bool Jal::multiply(DenseLinearSpace& retMat,
		   DenseLinearSpace& aMat,
		   SparseLinearSpace& bMat,
		   double* scalar )
{
  bool total_judge = SDPA_SUCCESS;

  // for SDP
  for (int l=0; l<bMat.SDP_sp_nBlock; ++l) {
    int index = bMat.SDP_sp_index[l];
    bool judge = Lal::multiply(retMat.SDP_block[index],
			       aMat.SDP_block[index],
			       bMat.SDP_sp_block[l],scalar);
    if (judge == SDPA_FAILURE) {
      total_judge = SDPA_FAILURE;
    }
  }

  // for SOCP
#if 0
  for (int l=0; l<bMat.SOCP_sp_nBlock; ++l) {
    int index = bMat.SOCP_sp_index[l];
    bool judge = Lal::multiply(retMat.SOCP_block[index],
			       aMat.SOCP_block[index],
			       bMat.SOCP_sp_block[l],scalar);
    if (judge == SDPA_FAILURE) {
      total_judge = SDPA_FAILURE;
    }
  }
#endif

  // for LP
  for (int l=0; l<bMat.LP_sp_nBlock; ++l) {
    int index = bMat.LP_sp_index[l];
    if (scalar == NULL) {
      retMat.LP_block[index] = 
	aMat.LP_block[index] * bMat.LP_sp_block[l];
    } else {
      retMat.LP_block[index] = 
	aMat.LP_block[index] * bMat.LP_sp_block[l] * (*scalar);
    }
  }

  return total_judge;
}
#endif  

  //  retMat = (A * B + B * A)/2
bool Jal::jordan_product(DenseLinearSpace& retMat,
			 DenseLinearSpace& aMat,
			 DenseLinearSpace& bMat)
{
  multiply(retMat,aMat,bMat);
  Lal::getSymmetrize(retMat);
  return SDPA_SUCCESS;
}

  //  retMat = A * B
bool Jal::ns_jordan_product(DenseLinearSpace& retMat,
			    DenseLinearSpace& aMat,
			    DenseLinearSpace& bMat)
{
  multiply(retMat,aMat,bMat);
  return SDPA_SUCCESS;
}
  
  //  retMat = A * B * A 
bool Jal::jordan_quadratic_product(DenseLinearSpace& retMat,
				   DenseLinearSpace& aMat,
				   DenseLinearSpace& bMat,
				   DenseLinearSpace& workMat)
{
  multiply(workMat,aMat,bMat);
  multiply(retMat,workMat,aMat);
  return SDPA_SUCCESS;
}
  
  //  retMat = (A * B * C + C * B * A)/2
bool Jal::jordan_triple_product(DenseLinearSpace& retMat,
				DenseLinearSpace& aMat,
				DenseLinearSpace& bMat,
				DenseLinearSpace& cMat,
				DenseLinearSpace& workMat)
{
  multiply(workMat,aMat,bMat);
  multiply(retMat,workMat,cMat);
  Lal::getSymmetrize(retMat);
  return SDPA_SUCCESS;

}

//  retMat = A * B * C
bool Jal::ns_jordan_triple_product(DenseLinearSpace& retMat,
				   DenseLinearSpace& aMat,
				   DenseLinearSpace& bMat,
				   DenseLinearSpace& cMat,
				   DenseLinearSpace& workMat)
{
  multiply(workMat,aMat,bMat);
  multiply(retMat,workMat,cMat);
  return SDPA_SUCCESS;
}

} // end of namespace 'sdpa'
