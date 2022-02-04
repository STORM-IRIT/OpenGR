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
// Here is the start of ``example1.cpp''
  
#include <cstdio>
#include <cstdlib>
#include <sdpa_call.h>

/*
*Example 2:
*mDim = 5, nBLOCK = 3, {2,3,-2}
   5  =  mDIM
   3  =  nBLOCK
   2    3   -2   = bLOCKsTRUCT
{1.1, -10, 6.6 , 19 , 4.1}
{
{ { -1.4, -3.2 },
  { -3.2,-28   }   }
{ { 15,  -12,    2.1 },
  {-12,   16,   -3.8 },
  {  2.1, -3.8, 15   }   }
  {  1.8, -4.0 }
}
{
{ {  0.5,  5.2 },
  {  5.2, -5.3 }   }
{ {  7.8, -2.4,  6.0 },
  { -2.4,  4.2,  6.5 },
  {  6.0,  6.5,  2.1 }   }
  { -4.5, -3.5 }
}
{
{ { 1.7,  7.0 },
  { 7.0, -9.3 }   }
{ {-1.9, -0.9, -1.3 },
  {-0.9, -0.8, -2.1 },
  {-1.3, -2.1,  4.0 }   }
  {-0.2, -3.7 }
}
{
{ { 6.3, -7.5 },
  {-7.5, -3.3 }   }
{ { 0.2,  8.8,  5.4 },
  { 8.8,  3.4, -0.4 },
  { 5.4, -0.4,  7.5 }   }
  {-3.3, -4.0 }
}
{
{ { -2.4, -2.5 },
  { -2.5, -2.9 }   }
{ {  3.4, -3.2, -4.5 },
  { -3.2,  3.0, -4.8 },
  { -4.5, -4.8,  3.6 }   }
  {  4.8 , 9.7 }
}
{
{ { -6.5, -5.4 },
  { -5.4, -6.6 }   }
{ {  6.7, -7.2, -3.6 },
  { -7.2,  7.3, -3.0 },
  { -3.6, -3.0, -1.4 }   }
  {  6.1, -1.5 }
}
*/

void printVector(double* ele, int dim, char* printFormat,
		 FILE* fpout);
void printMatrix(double* ele, int dim, char* printFormat,
		 FILE* fpout);
void printDimacsError(double dimacs_error[7],char* printFormat,
		      FILE* fpout);

int main ()
{
  SDPA::printSDPAVersion(stdout);
  SDPA	Problem1;
  Problem1.setDisplay(stdout);
  
  // All parameteres are renewed
  Problem1.setParameterType(SDPA::PARAMETER_DEFAULT);
  // Problem1.printParameters(stdout);

  int mDIM   = 5;
  int nBlock = 3;
  Problem1.inputConstraintNumber(mDIM);
  Problem1.inputBlockNumber(nBlock);
  // bLOCKsTRUCT :: 2(SDP)  3(SDP) -2(LP)
  Problem1.inputBlockSize(1,2);
  Problem1.inputBlockSize(2,3);
  Problem1.inputBlockSize(3,-2);
  Problem1.inputBlockType(1,SDPA::SDP);
  Problem1.inputBlockType(2,SDPA::SDP);
  Problem1.inputBlockType(3,SDPA::LP);

  Problem1.initializeUpperTriangleSpace();
  
  //cVECT = {1.1, -10, 6.6 , 19 , 4.1}
  Problem1.inputCVec(1,1.1);
  Problem1.inputCVec(2,-10);
  Problem1.inputCVec(3,6.6);
  Problem1.inputCVec(4,19);
  Problem1.inputCVec(5,4.1);

  // ---------  Input F_0  --------------------

  // 1st block
  // { { -1.4, -3.2},
  //   { -3.2,  -28} }
  
  Problem1.inputElement(0, 1, 1, 1, -1.4);
  Problem1.inputElement(0, 1, 1, 2, -3.2);
  Problem1.inputElement(0, 1, 2, 2,  -28);
  
  // 2nd block
  // { { 15,  -12,    2.1 },
  //   {-12,   16,   -3.8 },
  //   {  2.1, -3.8, 15   }   }

  Problem1.inputElement(0, 2, 1, 1,  15);
  Problem1.inputElement(0, 2, 1, 2, -12);
  Problem1.inputElement(0, 2, 1, 3, 2.1);
  Problem1.inputElement(0, 2, 2, 2,  16);
  Problem1.inputElement(0, 2, 2, 3,-3.8);
  Problem1.inputElement(0, 2, 3, 3,  15);

  // 3rd block
  // {  1.8, -4.0 }
  Problem1.inputElement(0, 3, 1, 1, 1.8);
  Problem1.inputElement(0, 3, 2, 2,-4.0);

  // ---------  Input F_1  --------------------

  // 1st block
  // { {  0.5,  5.2},
  //   {  5.2, -5.3} }
  
  Problem1.inputElement(1, 1, 1, 1, 0.5);
  Problem1.inputElement(1, 1, 1, 2, 5.2);
  Problem1.inputElement(1, 1, 2, 2,-5.3);
  
  // 2nd block
  // { { 7.8, -2.4, 6.0 },
  //   {-2.4,  4.2, 6.5 },
  //   { 6.0,  6.5, 2.1 } }

  Problem1.inputElement(1, 2, 1, 1, 7.8);
  Problem1.inputElement(1, 2, 1, 2,-2.4);
  Problem1.inputElement(1, 2, 1, 3, 6.0);
  Problem1.inputElement(1, 2, 2, 2, 4.2);
  Problem1.inputElement(1, 2, 2, 3, 6.5);
  Problem1.inputElement(1, 2, 3, 3, 2.1);

  // 3rd block
  // { -4.5, -3.5 }
  Problem1.inputElement(1, 3, 1, 1, -4.5);
  Problem1.inputElement(1, 3, 2, 2, -3.5);

  // ---------  Input F_2  --------------------

  // 1st block
  // { { 1.7, 7.0},
  //   { 7.0,-9.3} }
  
  Problem1.inputElement(2, 1, 1, 1,  1.7);
  Problem1.inputElement(2, 1, 1, 2,  7.0);
  Problem1.inputElement(2, 1, 2, 2, -9.3);
  
  // 2nd block
  // { {-1.9, -0.9, -1.3 },
  //   {-0.9, -0.8, -2.1 },
  //   {-1.3, -2.1,  4.0 } }

  Problem1.inputElement(2, 2, 1, 1, -1.9);
  Problem1.inputElement(2, 2, 1, 2, -0.9);
  Problem1.inputElement(2, 2, 1, 3, -1.3);
  Problem1.inputElement(2, 2, 2, 2, -0.8);
  Problem1.inputElement(2, 2, 2, 3, -2.1);
  Problem1.inputElement(2, 2, 3, 3,  4.0);

  // 3rd block
  // { -0.2, -3.7 }
  Problem1.inputElement(2, 3, 1, 1, -0.2);
  Problem1.inputElement(2, 3, 2, 2, -3.7);

  // ---------  Input F_3  --------------------

  // 1st block
  // { {  6.3, -7.5},
  //   { -7.5, -3.3} }
  
  Problem1.inputElement(3, 1, 1, 1,  6.3);
  Problem1.inputElement(3, 1, 1, 2, -7.5);
  Problem1.inputElement(3, 1, 2, 2, -3.3);
  
  // 2nd block
  // { {  0.2,  8.8,  5.4 },
  //   {  8.8,  3.4, -0.4 },
  //   {  5.4, -0.4,  7.5 } }

  Problem1.inputElement(3, 2, 1, 1,  0.2);
  Problem1.inputElement(3, 2, 1, 2,  8.8);
  Problem1.inputElement(3, 2, 1, 3,  5.4);
  Problem1.inputElement(3, 2, 2, 2,  3.4);
  Problem1.inputElement(3, 2, 2, 3, -0.4);
  Problem1.inputElement(3, 2, 3, 3,  7.5);

  // 3rd block
  // { -3.3, -4.0 }
  Problem1.inputElement(3, 3, 1, 1, -3.3);
  Problem1.inputElement(3, 3, 2, 2, -4.0);

  // ---------  Input F_4  --------------------

  // 1st block
  // { { -2.4, -2.5},
  //   { -2.5, -2.9} }
  
  Problem1.inputElement(4, 1, 1, 1, -2.4);
  Problem1.inputElement(4, 1, 1, 2, -2.5);
  Problem1.inputElement(4, 1, 2, 2, -2.9);
  
  // 2nd block
  // { {  3.4, -3.2, -4.5 },
  //   { -3.2,  3.0, -4.8 },
  //   { -4.5, -4.8,  3.6 } }

  Problem1.inputElement(4, 2, 1, 1,  3.4);
  Problem1.inputElement(4, 2, 1, 2, -3.2);
  Problem1.inputElement(4, 2, 1, 3, -4.5);
  Problem1.inputElement(4, 2, 2, 2,  3.0);
  Problem1.inputElement(4, 2, 2, 3, -4.8);
  Problem1.inputElement(4, 2, 3, 3,  3.6);

  // 3rd block
  // { 4.8, 9.7 }
  Problem1.inputElement(4, 3, 1, 1, 4.8);
  Problem1.inputElement(4, 3, 2, 2, 9.7);

  // ---------  Input F_5  --------------------

  // 1st block
  // { { -6.5, -5.4},
  //   { -5.4, -6.6} }
  
  Problem1.inputElement(5, 1, 1, 1, -6.5);
  Problem1.inputElement(5, 1, 1, 2, -5.4);
  Problem1.inputElement(5, 1, 2, 2, -6.6);
  
  // 2nd block
  // { {  6.7, -7.2, -3.6 },
  //   { -7.2,  7.3, -3.0 },
  //   { -3.6, -3.0, -1.4 } }

  Problem1.inputElement(5, 2, 1, 1,  6.7);
  Problem1.inputElement(5, 2, 1, 2, -7.2);
  Problem1.inputElement(5, 2, 1, 3, -3.6);
  Problem1.inputElement(5, 2, 2, 2,  7.3);
  Problem1.inputElement(5, 2, 2, 3, -3.0);
  Problem1.inputElement(5, 2, 3, 3, -1.4);

  // 3rd block
  // {  6.1, -1.5 }
  Problem1.inputElement(5, 3, 1, 1, 6.1);
  Problem1.inputElement(5, 3, 2, 2,-1.5);

  Problem1.initializeUpperTriangle();
  Problem1.initializeSolve();
  
  // if necessary, dump input data and initial point
  // Problem1.writeInputSparse((char*)"tmp.dat-s",(char*)"%+8.3e");
  // Problem1.writeInitSparse((char*)"tmp.ini-s",(char*)"%+8.3e");
  
  Problem1.solve();
	
  fprintf(stdout, "\nStop iteration = %d\n",
	  Problem1.getIteration());
  char phase_string[30];
  Problem1.getPhaseString(phase_string);
  fprintf(stdout, "Phase          = %s\n", phase_string);
  fprintf(stdout, "objValPrimal   = %+10.6e\n",
	  Problem1.getPrimalObj());
  fprintf(stdout, "objValDual     = %+10.6e\n",
	  Problem1.getDualObj());
  fprintf(stdout, "p. feas. error = %+10.6e\n",
	  Problem1.getPrimalError());
  fprintf(stdout, "d. feas. error = %+10.6e\n\n",
	  Problem1.getDualError());

  
  fprintf(stdout, "xVec = \n");
  // Problem1.printResultXVec();
  printVector(Problem1.getResultXVec(),
	      Problem1.getConstraintNumber(), (char*)"%+8.3e",
	      stdout);
  
  fprintf(stdout, "xMat = \n");
  // Problem1.printResultXMat();
  for (int l=0; l<Problem1.getBlockNumber(); ++l) {
    if (Problem1.getBlockType(l+1) == SDPA::SDP) {
      printMatrix(Problem1.getResultXMat(l+1),
		  Problem1.getBlockSize(l+1), (char*)"%+8.3e",
		  stdout);
    }
    else if (Problem1.getBlockType(l+1) == SDPA::SOCP) {
      printf("current version does not support SOCP\n");
    }
    if (Problem1.getBlockType(l+1) == SDPA::LP) {
      printVector(Problem1.getResultXMat(l+1),
		  Problem1.getBlockSize(l+1), (char*)"%+8.3e",
		  stdout);
    }
  }
		  
  fprintf(stdout, "yMat = \n");
  // Problem1.printResultYMat();
  for (int l=0; l<Problem1.getBlockNumber(); ++l) {
    if (Problem1.getBlockType(l+1) == SDPA::SDP) {
      printMatrix(Problem1.getResultYMat(l+1),
		  Problem1.getBlockSize(l+1), (char*)"%+8.3e",
		  stdout);
    }
    else if (Problem1.getBlockType(l+1) == SDPA::SOCP) {
      printf("current version does not support SOCP\n");
    }
    if (Problem1.getBlockType(l+1) == SDPA::LP) {
      printVector(Problem1.getResultYMat(l+1),
		  Problem1.getBlockSize(l+1), (char*)"%+8.3e",
		  stdout);
    }
  }

  double dimacs_error[7];
  Problem1.getDimacsError(dimacs_error);
  printDimacsError(dimacs_error,(char*)"%+8.3e",stdout);

  // Problem1.printComputationTime(stdout);

  Problem1.terminate();
  exit(0);
};	

void printVector(double* ele, int dim, char* printFormat, FILE* fpout)
{
  fprintf(fpout,"[ ");
  for (int k=0; k<dim-1; ++k) {
    fprintf(fpout,printFormat,ele[k]);
    fprintf(fpout," ");
  }
  fprintf(fpout,printFormat,ele[dim-1]);
  fprintf(fpout,"]; \n");
}

void printMatrix(double* ele, int dim, char* printFormat, FILE* fpout)
{
  fprintf(fpout,"[\n");
  for (int i=0; i<dim; ++i) {
    fprintf(fpout,"[ ");
    for (int j=0; j<dim-1; ++j) {
      fprintf(fpout,printFormat,ele[i+dim*j]);
      fprintf(fpout," ");
    }
    fprintf(fpout,printFormat,ele[i+dim*(dim-1)]);
    fprintf(fpout,"]; \n");
  }
  fprintf(fpout,"]; \n");
}

void printDimacsError(double dimacs_error[7],char* printFormat,
		      FILE* fpout)
{
  fprintf(fpout,  "\n");
  fprintf(fpout,  "* DIMACS_ERRORS * \n");
  fprintf(fpout,  "err1 = ");
  fprintf(fpout,  printFormat, dimacs_error[1]);
  fprintf(fpout, "  [||Ax-b|| / (1+||b||_1)]\n");
  fprintf(fpout,  "err2 = ");
  fprintf(fpout,  printFormat, dimacs_error[2]);
  fprintf(fpout, "  [max(0, -lambda(x)/(1+||b||_1))]\n");
  fprintf(fpout,  "err3 = ");
  fprintf(fpout,  printFormat, dimacs_error[3]);
  fprintf(fpout, "  [||A^Ty + z - c || / (1+||c||_1)]\n");
  fprintf(fpout,  "err4 = ");
  fprintf(fpout,  printFormat, dimacs_error[4]);
  fprintf(fpout, "  [max(0, -lambda(z)/(1+||c||_1))]\n");
  fprintf(fpout,  "err5 = ");
  fprintf(fpout,  printFormat, dimacs_error[5]);
  fprintf(fpout, "  [(<c,x> - <b,y>) / (1 + |<c,x>| + |<b,y>|)]\n");
  fprintf(fpout,  "err6 = ");
  fprintf(fpout,  printFormat, dimacs_error[6]);
  fprintf(fpout, "  [<x,z> / (1 + |<c,x>| + |<b,y>|)]\n");
  fprintf(fpout,  "\n");
}

