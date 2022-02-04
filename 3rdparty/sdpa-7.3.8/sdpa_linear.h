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

#ifndef __sdpa_linear_h__
#define __sdpa_linear_h__

#include "sdpa_struct.h"

namespace sdpa {

class Lal
{
public:
// calculate the minimum eigenvalue of lMat*xMat*(lMat^T)
// by Lanczos metnod
  static double getMinEigen(DenseMatrix& lMat, DenseMatrix& xMat,
			    DenseMatrix& Q,
			    Vector& out, Vector& b,  Vector& r,
			    Vector& q, Vector& qold, 
			    Vector& w, Vector& tmp,
			    Vector& diagVec, Vector& diagVec2,
			    Vector& workVec);

// caluculate all eigenvalues of aMat  by QR method 
  static double getMinEigenValue(DenseMatrix& aMat,
				 Vector& eigenVec,
				 Vector& workVec);

  static double getOneNorm(Vector& b);
  static double getOneNorm(SparseMatrix& C);
  static double getOneNorm(SparseLinearSpace& C);

  static double getTwoNorm(Vector& b);
  static double getTwoNorm(DenseMatrix& X);
  static double getTwoNorm(DenseLinearSpace& X);

  
  static bool getInnerProduct(double& ret,
			      Vector& aVec, Vector& bVec);
  static bool getInnerProduct(double& ret,
			      BlockVector& aVec,
			      BlockVector& bVec);
  static bool getInnerProduct(double& ret,
			      DenseMatrix& aMat,
			      DenseMatrix& bMat);
  static bool getInnerProduct(double& ret,
			      SparseMatrix& aMat,
			      DenseMatrix&  bMat);

  static bool getCholesky(DenseMatrix& retMat, DenseMatrix& aMat);

// nakata 2004/12/01 
// diagonal part of Cholesky matrix is set these inverse.
  static bool getCholesky(SparseMatrix& aMat, int* diagonalIndex);

  static bool getInvLowTriangularMatrix(DenseMatrix& retMat,
					DenseMatrix& aMat);
  static bool getSymmetrize(DenseMatrix& aMat);

  static bool getTranspose(DenseMatrix& retMat,
			   DenseMatrix& aMat);

  static int rdpotf2_(char*uplo, int *n, double *a, int *lda, int *info);
  static int rdpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
  static bool choleskyFactorWithAdjust(DenseMatrix& aMat);
  
  static bool solveSystems(Vector& xVec,
			   DenseMatrix& aMat, Vector& bVec);
  // solve aMat * xVec = bVec
  // aMat must be Cholesky Factorized.

// nakata 2004/12/01 
  static bool solveSystems(Vector& xVec,
			   SparseMatrix& aMat, Vector& bVec);
  // solve aMat * xVec = bVec
  // aMat must be Cholesky Factorized.

  static bool getSymmetrize(DenseLinearSpace& aMat);

  static bool getTranspose(DenseLinearSpace& retMat,
			   DenseLinearSpace& aMat);

  static bool multiply(DenseMatrix& retMat,
		       DenseMatrix& aMat, DenseMatrix& bMat,
		       double* scalar = NULL);
  static bool multiply(DenseMatrix& retMat,
		       SparseMatrix& aMat, DenseMatrix& bMat,
		       double* scalar = NULL);
  static bool multiply(DenseMatrix& retMat,
		       DenseMatrix& aMat, SparseMatrix& bMat,
		       double* scalar = NULL);
  static bool multiply(DenseMatrix& retMat,
		       DenseMatrix& aMat, double* scalar = NULL);
  static bool multiply(Vector& retVec,
		       Vector& aVec, double* scalar = NULL);
  static bool multiply(BlockVector& retVec,
		       BlockVector& aVec,
		       double* scalar = NULL);
  static bool multiply(Vector& retVec,
		       DenseMatrix& aMat, Vector& bVec,
		       double* scalar = NULL);
  // ret = aMat**T * bMat
  static bool tran_multiply(DenseMatrix& retMat,
			    DenseMatrix& aMat, DenseMatrix& bMat,
			    double* scalar = NULL);
  // ret = aMat * bMat**T
  static bool multiply_tran(DenseMatrix& retMat,
			    DenseMatrix& aMat, DenseMatrix& bMat,
			    double* scalar = NULL);
  // ret = a + (*scalar)*b
  static bool plus(Vector& retVec, Vector& aVec,
		   Vector& bVec, double* scalar = NULL);
  static bool plus(DenseMatrix& retMat,
		   DenseMatrix& aMat, DenseMatrix& bMat,
		   double* scalar = NULL);
  static bool plus(DenseMatrix& retMat,
		   SparseMatrix& aMat, DenseMatrix& bMat,
		   double* scalar = NULL);
  static bool plus(DenseMatrix& retMat,
		   DenseMatrix& aMat, SparseMatrix& bMat,
		   double* scalar = NULL);
  
  static bool plus(BlockVector& retVec,
		   BlockVector& aVec,
		   BlockVector& bVec, double* scalar = NULL);

  // ret = a '*' (*scalar)
  static bool let(Vector& retVec, const char eq,
		  Vector& aVec, const char op,
		  double* scalar = NULL);

  // ret = a '*' (*scalar)
  static bool let(BlockVector& retVec, const char eq,
		  BlockVector& aVec, const char op,
		  double* scalar = NULL);

  // ret = a '*' (*scalar)
  static bool let(DenseMatrix& retMat, const char eq,
		  DenseMatrix& aMat, const char op,
		  double* scalar = NULL);

  // ret = a '+' '-' b*(*scalar)
  static bool let(Vector& retVec, const char eq,
		  Vector& aVec, const char op,
		  Vector& bVec, double* scalar = NULL);

  // ret = a '+' '-' '*' 't' 'T' b*(*scalar)
  static bool let(DenseMatrix& retMat, const char eq,
		  DenseMatrix& aMat, const char op,
		  DenseMatrix& bMat, double* scalar = NULL);

  // ret = a '+' '-' '*' b*(*scalar)
  static bool let(DenseMatrix& retMat, const char eq,
		  SparseMatrix& aMat, const char op,
		  DenseMatrix& bMat, double* scalar = NULL);

  // ret = a '+' '-' '*' b*(*scalar)
  static bool let(DenseMatrix& retMat, const char eq,
		  DenseMatrix& aMat, const char op,
		  SparseMatrix& bMat, double* scalar = NULL);

  // ret = aMat '*' '/' bVec
  static bool let(Vector& rVec, const char eq,
		  DenseMatrix& aMat, const char op,
		  Vector& bVec);

  // nakata 2004/12/01
  // ret = aMat '/' bVec
  static bool let(Vector& rVec, const char eq,
		  SparseMatrix& aMat, const char op,
		  Vector& bVec);

  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  Vector& aVec, const char op,
		  Vector& bVec);
  
  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  DenseMatrix& aMat, const char op,
		  DenseMatrix& bMat);
  
  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  DenseMatrix& aMat, const char op,
		  SparseMatrix& bMat);
  
  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  SparseMatrix& aMat, const char op,
		  DenseMatrix& bMat);

  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  BlockVector& aVec, const char op,
		  BlockVector& bVec);
  
  /////////////////////////////////////////////////////////////////////

  static bool getInnerProduct(double& ret,
			      DenseLinearSpace& aMat,
			      DenseLinearSpace&  bMat);

  static bool getInnerProduct(double& ret,
			      SparseLinearSpace& aMat,
			      DenseLinearSpace&  bMat);

  // ret = a (*scalar)*b
  static bool multiply(DenseLinearSpace& retMat,
		       DenseLinearSpace& aMat,
		       double* scalar = NULL);
  // ret = a + (*scalar)*b
  static bool plus(DenseLinearSpace& retMat,
		   DenseLinearSpace& aMat,
		   DenseLinearSpace& bMat,
		   double* scalar = NULL);
// CAUTION!!! We don't initialize retMat to zero matrix for efficiently.
  static bool plus(DenseLinearSpace& retMat,
		   SparseLinearSpace& aMat,
		   DenseLinearSpace& bMat,
		   double* scalar = NULL);
// CAUTION!!! We don't initialize retMat to zero matrix for efficiently.
  static bool plus(DenseLinearSpace& retMat,
		   DenseLinearSpace& aMat,
		   SparseLinearSpace& bMat,
		   double* scalar = NULL);

  // ret = a '*' (*scalar)
  static bool let(DenseLinearSpace& retMat, const char eq,
		  DenseLinearSpace& aMat, const char op,
		  double* scalar = NULL);

  // ret = a '+' '-' b*(*scalar)
  static bool let(DenseLinearSpace& retMat, const char eq,
		  DenseLinearSpace& aMat, const char op,
		  DenseLinearSpace& bMat, double* scalar = NULL);

  // ret = a '+' '-' b*(*scalar)
  static bool let(DenseLinearSpace& retMat, const char eq,
		  SparseLinearSpace& aMat, const char op,
		  DenseLinearSpace& bMat, double* scalar = NULL);

  // ret = a '+' '-' '*' b*(*scalar)
  static bool let(DenseLinearSpace& retMat, const char eq,
		  DenseLinearSpace& aMat, const char op,
		  SparseLinearSpace& bMat, double* scalar = NULL);

  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  DenseLinearSpace& aMat, const char op,
		  DenseLinearSpace& bMat);

  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  SparseLinearSpace& aMat, const char op,
		  DenseLinearSpace& bMat);

  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  DenseLinearSpace& aMat, const char op,
		  SparseLinearSpace& bMat);

};

} // end of namespace 'sdpa'

#endif // __sdpa_linear_h__
