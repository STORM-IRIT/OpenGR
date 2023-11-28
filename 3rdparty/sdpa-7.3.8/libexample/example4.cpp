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
// Here is the start of ``example4.cpp''
  
#include <cstdio>
#include <cstdlib>
#include <sdpa_call.h>

/*
example1.dat:

"Example 1: mDim = 3, nBLOCK = 1, {2}"
   3  =  mDIM
   1  =  nBLOCK
   2  = bLOCKsTRUCT
{48, -8, 20}
{ {-11,  0}, { 0, 23} }
{ { 10,  4}, { 4,  0} }
{ {  0,  0}, { 0, -8} }
{ {  0, -8}, {-8, -2} }
*/

void printVector(double* ele, int dim, char* printFormat,
		 FILE* fpout);
void printMatrix(double* ele, int dim, char* printFormat,
		 FILE* fpout);
void printDimacsError(double dimacs_error[7],char* printFormat,
		      FILE* fpout);
void printResult(SDPA& Problem1);

int main ()
{
  SDPA::printSDPAVersion(stdout);
  SDPA	Problem1;
  Problem1.setDisplay(stdout);

  // All parameteres are renewed
  Problem1.setParameterType(SDPA::PARAMETER_DEFAULT);

  // If necessary, each parameter can be set independently
  
  Problem1.setParameterEpsilonStar(1.0e-7);
  Problem1.printParameters(stdout);

  int mDIM   = 3;
  int nBlock = 1;
  Problem1.inputConstraintNumber(mDIM);
  Problem1.inputBlockNumber(nBlock);
  Problem1.inputBlockSize(1,2);
  Problem1.inputBlockType(1,SDPA::SDP);

  Problem1.initializeUpperTriangleSpace();

  Problem1.inputCVec(1,48);
  Problem1.inputCVec(2,-8);
  Problem1.inputCVec(3,20);

  Problem1.inputElement(0, 1, 1, 1, -11);
  Problem1.inputElement(0, 1, 2, 2,  23);
  
  Problem1.inputElement(1, 1, 1, 1,  10);
  Problem1.inputElement(1, 1, 1, 2,   4);
  
  Problem1.inputElement(2, 1, 2, 2,  -8);
  
  Problem1.inputElement(3, 1, 1, 2,  -8);
  Problem1.inputElement(3, 1, 2, 2,  -2);

  Problem1.initializeUpperTriangle();
  Problem1.initializeSolve();

  // if necessary, dump input data and initial point
  // Problem1.writeInputSparse((char*)"tmp.dat-s",(char*)"%+8.3e");
  // Problem1.writeInitSparse((char*)"tmp.ini-s",(char*)"%+8.3e");

  Problem1.solve();
  printResult(Problem1);

  // Solve the problem again with smaller epsilonStar
  Problem1.setParameterEpsilonStar(1.0e-8);

  Problem1.setInitPoint(true);
  Problem1.copyCurrentToInit();

  Problem1.solve();
  printResult(Problem1);

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

void printResult(SDPA& Problem1)
{
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

  #if 0
  
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
  #endif

  double dimacs_error[7];
  Problem1.getDimacsError(dimacs_error);
  printDimacsError(dimacs_error,(char*)"%+8.3e",stdout);

}
