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
/*  sdpa_algebra.h

LAPACK+BLAS definitions wrapper

Define macros to mangle the given C identifier (in lower and upper
case), which must not contain underscores, for linking with Fortran.

*/

#ifndef __sdpa_algebra_h__
#define __sdpa_algebra_h__

#define FC_RET_I int
#define FC_RET_D double

#if defined(__APPLE__) // Dirty...
#define FC_FUNC(name,NAME) name ## _
#endif

#define dtrsm_fc  FC_FUNC (dtrsm, DTRSM)
#define dsyrk_fc  FC_FUNC (dsyrk, DSYRK)
#define dcopy_fc  FC_FUNC (dcopy, DCOPY)
#define daxpy_fc  FC_FUNC (daxpy, DAXPY)
#define dgemm_fc  FC_FUNC (dgemm, DGEMM)
#define dgemv_fc  FC_FUNC (dgemv, DGEMV)
#define dscal_fc  FC_FUNC (dscal, DSCAL)
#define dtrsv_fc  FC_FUNC (dtrsv, DTRSV)
#define dtrmv_fc  FC_FUNC (dtrmv, DTRMV)
#define ddot_fc   FC_FUNC (ddot, DDOT)
#define dtrmm_fc  FC_FUNC (dtrmm, DTRMM)
#define ilaenv_fc FC_FUNC (ilaenv, ILAENV)
#define dsteqr_fc FC_FUNC (dsteqr, DSTEQR)
#define dsyev_fc  FC_FUNC (dsyev, DSYEV)
#define dpotrf_fc FC_FUNC (dpotrf, DPORTRF)


extern "C"
{
// BLAS
  FC_RET_I  dtrsm_fc
      (char* side, char* uplo, char* trans, char* diag,
       int* M, int* N,
       double* alpha,
       double* A, int* lda,
       double* B, int* ldb, int side_len,
       int uplo_len, int trans_len, int diag_len);

  FC_RET_I  dsyrk_fc
      (char* uplo, char* trans, int* N, int* K,
       double* alpha,
       double* A, int* lda,
       double* beta,
       double* C, int* ldc, int uplo_len, int trans_len);

  FC_RET_I  dcopy_fc
      (int* N,
       double* X, int* incX,
       double* Y, int* incY);

  FC_RET_I  daxpy_fc
      (int* N,
       double* alpha,
       double* X, int* incX,
       double* Y, int* incY);

  FC_RET_I  dgemm_fc
      (char* transA, char* transB, int* M, int* N, int* K,
       double* alpha,
       double* A, int* lda,
       double* B, int* ldb,
       double* beta,
       double* C, int* ldc, int transA_len, int transB_len);

  FC_RET_I  dgemv_fc
      (char* trans, int* M, int* N,
       double* alpha,
       double* A, int* lda,
       double* X, int* incX,
       double* beta,
       double* Y, int* incY, int trans_len);

  FC_RET_I  dscal_fc
      (int* N,
       double* alpha,
       double* X, int* incX);

  FC_RET_I  dtrsv_fc
      (char* uplo, char* trans, char* diag, int* N,
       double* A, int* lda,
       double* X, int* incX, int uplo_len,
       int trans_len, int diag_len);

  FC_RET_I  dtrmv_fc
      (char* uplo, char *trans, char* diag, int *N,  
       double *A, int *lda, 
       double *X, int *incX, int uplo_len, int trans_len, int diag_len);

  FC_RET_D  ddot_fc
      (int* N, double* X, int* incX, double* Y, int* incY);

  FC_RET_I  dtrmm_fc
      (char* side, char* uplo, char* trans, char* diag, 
       int* M, int* N,
       double* alpha,
       double* A, int* lda,
       double* B, int* ldb, int side_len, int uplo_len,
       int trans_len, int diag_len);

// LAPACK

  FC_RET_I  ilaenv_fc
      (int *ispec, char *name, char *opts, int *n1, 
	int *n2, int *n3, int *n4, int name_len, int opts_len);

  FC_RET_I  dsteqr_fc
      (char *compz, int *n, double *d, 
	double *e, double *z, int *ldz, double *work, 
	int *info, int compz_len);

  FC_RET_I  dsyev_fc
      (char *jobz, char *uplo, int *n, double *a,
        int *lda, double *w, double *work, int *lwork, 
	int *info, int jobz_len, int uplo_len);

  FC_RET_I  dpotrf_fc
     (char *uplo, int *n, double *a, int *lda,
      int *info, int uplo_len);
}

#endif // __sdpa_algebra_h__
