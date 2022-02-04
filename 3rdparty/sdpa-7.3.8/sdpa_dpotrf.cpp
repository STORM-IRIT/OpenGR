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
/*-----------------------------------------------
  rsdpa_dpotrf.cpp
  modification of ATL_dpotrfL
  for dealing with numerical error
  in diagonal elements.
  
  int rATL_dpotrfL(int N, double *A,int lda)

  modified by Makoto Yamshita 2002.07.11
-----------------------------------------------*/
#define POTRF_NONZERO (1.0e-14)
#define POTRF_ASSIGN  (1.0e+100)
#define POTRF_LIMIT   (-1.0e-6)

/*
 *             Automatically Tuned Linear Algebra Software v3.4.0
 *                    (C) Copyright 1999 R. Clint Whaley
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions, and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *   3. The name of the ATLAS group or the names of its contributers may
 *      not be used to endorse or promote products derived from this
 *      software without specific written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE ATLAS GROUP OR ITS CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */



#include "sdpa_include.h"
#include "sdpa_algebra.h"
#if 0
#define CHOLESKY_ADJUST(val) rMessage("Choleksy adjust from " << val << " to " << POTRF_NONZERO);
#else
#define CHOLESKY_ADJUST(val)       ;
#endif

namespace sdpa {

extern "C" {
static int potrf4(double* A,const int n)
{
  double* A1 = A+n+1;
  double* A2 = A1+n+1;
  double* A3 = A2+n+1;
  double L11 = *A;
  double L21 = A[1], L22 = *A1;
  double L31 = A[2], L32 = A1[1], L33 = *A2;
  double L41 = A[3], L42 = A1[2], L43 = A2[1], L44 = *A3;

  if (L11 < POTRF_LIMIT) {
    return 1;
  }
  if (L11 < POTRF_NONZERO) {
    CHOLESKY_ADJUST(L11);
    L11 = POTRF_ASSIGN;
  }

  *A = L11 = sqrt(L11);
  L11 = 1.0/L11;
  L21 *= L11;
  L31 *= L11;
  L41 *= L11;

  L22 -= L21*L21;
  if (L22 < POTRF_LIMIT) {
    return 2;
  }
  if (L22 < POTRF_NONZERO) {
    CHOLESKY_ADJUST(L22);
    L22 = POTRF_ASSIGN;
  }
  *A1 = L22 = sqrt(L22);
  L22 = 1.0/L22;
  L32 = (L32 - L31*L21)*L22;
  L42 = (L42 - L41*L21)*L22;

  L33 -= L31*L31 + L32*L32;
  if (L33 < POTRF_LIMIT) {
    return 3;
  }
  if (L33 < POTRF_NONZERO) {
    CHOLESKY_ADJUST(L33);
    L33 = POTRF_ASSIGN;
  }
  *A2 = L33 = sqrt(L33);
  L43 = (L43-L41*L31-L42*L32)/L33;
  L44 -= L41*L41 + L42*L42 + L43*L43;
  if (L44 < POTRF_LIMIT) {
    return 4;
  }
  if (L44 < POTRF_NONZERO) {
    CHOLESKY_ADJUST(L44);
    L44 = POTRF_ASSIGN;
  }
  *A3 = sqrt(L44);
  
  A[1] = L21;
  A[2] = L31; A1[1] = L32;
  A[3] = L41; A1[2] = L42; A2[1] = L43;
  return 0;
}

static int potrf3(double* A,const int n)
{
  double* A1 = A+n+1;
  double* A2 = A1+n+1;
  double L11 = *A;
  double L21 = A[1], L22 = *A1;
  double L31 = A[2], L32 = A1[1], L33 = *A2;

  if (L11 < POTRF_LIMIT) {
    return 1;
  }
  if (L11 < POTRF_NONZERO) {
    CHOLESKY_ADJUST(L11);
    L11 = POTRF_ASSIGN;
  }

  *A = L11 = sqrt(L11);
  L11 = 1.0/L11;
  L21 *= L11;
  L31 *= L11;

  L22 -= L21*L21;
  if (L22 < POTRF_LIMIT) {
    return 2;
  }
  if (L22 < POTRF_NONZERO) {
    CHOLESKY_ADJUST(L22);
    L22 = POTRF_ASSIGN;
  }
  L22 = sqrt(L22);
  L32 = (L32 - L31*L21)/L22;

  L33 -= L31*L31 + L32*L32;
  if (L33 < POTRF_LIMIT) {
    return 3;
  }
  if (L33 < POTRF_NONZERO) {
    CHOLESKY_ADJUST(L33);
    L33 = POTRF_ASSIGN;
  }
  *A2 = sqrt(L33);
  
  A[1] = L21; *A1   = L22;
  A[2] = L31; A1[1] = L32; 
  return 0;
}

static int potrf2(double* A,const int n)
{
  double* A1 = A+n+1;
  double L11 = *A;
  double L21 = A[1], L22 = *A1;

  if (L11 < POTRF_LIMIT) {
    return 1;
  }
  if (L11 < POTRF_NONZERO) {
    CHOLESKY_ADJUST(L11);
    L11 = POTRF_ASSIGN;
  }

  *A = L11 = sqrt(L11);
  L21 /= L11;
  L22 -= L21*L21;
  if (L22 < POTRF_LIMIT) {
    return 2;
  }
  if (L22 < POTRF_NONZERO) {
    CHOLESKY_ADJUST(L22);
    L22 = POTRF_ASSIGN;
  }
  *A = L11;
  A[1] = L21; *A1 = sqrt(L22);

  return 0;
}

int rATL_dpotrfL(int N, double *A,int lda)
{
  double *An, *Ar;
  int Nleft, Nright, ierr;

  if (N > 4) {
    Nleft = N >> 1;
#if 0
    int nb = ilaenv_fc(&IONE, "DPOTRF", "L", &N,
		     &IMONE,&IONE, &IMONE, strlen("DPOTRF"), strlen("L"));
    if (Nleft > nb<<1) Nleft = (Nleft/nb)*nb;
#endif
#if 0
    if (Nleft > 64) {
	Nleft = 64;
    }
#endif
    Nright = N - Nleft;
    ierr = rATL_dpotrfL(Nleft, A,lda);
    if (!ierr) {
      Ar = A + Nleft;
      An = Ar + lda * Nleft;
      dtrsm_fc ((char *)"R",(char *)"L",(char *)"T",(char *)"N",
		 &Nright,&Nleft,&DONE,A,&lda,
		 Ar, &lda, strlen("R"),strlen("L"),
		 strlen("T"),strlen("N"));
      dsyrk_fc ((char *)"L",(char *)"N",&Nright,&Nleft,&DMONE,
		 Ar, &lda, &DONE,An,&lda,strlen("L"),strlen("N"));
      ierr = rATL_dpotrfL(Nright, An,lda);
      if (ierr) return(ierr+Nleft);
    }
    else return(ierr);
  }
  else if (N==4) return(potrf4(A,lda));
  else if (N==3) return(potrf3(A,lda));
  else if (N==2) return(potrf2(A,lda));
  else if (N==1) {
    if (*A < POTRF_LIMIT) {
      return 1;
    }
    if (*A < POTRF_NONZERO) {
      CHOLESKY_ADJUST(*A);
      *A = POTRF_ASSIGN;
    }
    *A = sqrt(*A);
  }
  return(0);
}

void rdpotrfl_(int* N, double *A,int* lda,int* info)
{
  *info = rATL_dpotrfL(*N,A,*lda);
}

}; // end of extern "C"

} // end of namespace 'sdpa'
