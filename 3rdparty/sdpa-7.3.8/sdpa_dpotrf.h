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
  int rATL_dpotrfL(int N, double *A,int lda)
  $Id: rsdpa_dpotrf.h,v 1.2 2004/09/01 06:34:12 makoto Exp $
-----------------------------------------------*/

#ifndef __sdpa_dpotrf_h__
#define __sdpa_dpotrf_h__

namespace sdpa {

#ifdef __cplusplus
extern "C" int rATL_dpotrfL(int N, double *A,int lda);
#else
extern int rATL_dpotrfL(int N, double *A,int lda);
#endif

} // end of namespace 'sdpa'

#endif // __sdpa_dpotrf_h__


