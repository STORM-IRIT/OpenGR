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

#ifndef __sdpa_io_h__
#define __sdpa_io_h__

#define lengthOfString 256

#include "sdpa_block.h"
#include "sdpa_parts.h"
  
namespace sdpa {

class IO
{
public:
  static void read(FILE* fpData, FILE* fpout, int& m, char* str);
  static void read(FILE* fpData, int& nBlock);
  static void read(FILE* fpData, BlockStruct& bs);
  static void read(FILE* fpData, Vector& b);
  static void read(FILE* fpData, DenseLinearSpace& xMat,
		   Vector& yVec, DenseLinearSpace& zMat,
		   BlockStruct& bs, bool inputSparse);
  static void read(FILE* fpData, int m,
		   BlockStruct& bs,
		   InputData& inputData, bool isDataSparse);


  // 2008/02/27 kazuhide nakata   
  // without LP_ANonZeroCount
  static void setBlockStruct(FILE* fpData, InputData& inputData,
			     int m, BlockStruct& bs,
                             long position, bool isDataSparse);
  
  // 2008/02/27 kazuhide nakata   
  // without LP_ANonZeroCount
  static void setElement(FILE* fpData, InputData& inputData, int m,
			 BlockStruct& bs,
                         long position, bool isDataSparse);

  static void printHeader(FILE* fpout, FILE* Display);

  static void printOneIteration(int pIteration,
				AverageComplementarity& mu,
				RatioInitResCurrentRes& theta,
				SolveInfo& solveInfo,
				StepLength& alpha,
				DirectionParameter& beta,
				FILE* fpout,
				FILE* Display);

  static void printLastInfo(int pIteration,
			    AverageComplementarity& mu,
			    RatioInitResCurrentRes& theta,
			    SolveInfo& solveInfo,
			    StepLength& alpha,
			    DirectionParameter& beta,
			    Residuals& currentRes,
			    Phase & phase,
			    Solutions& currentPt,
			    InputData& inputData,
                            WorkVariables& work,
			    double cputime,
			    ComputeTime& com,
			    Parameter& param,
			    FILE* fpout,
			    FILE* Display,
			    bool printTime = true);

  static void computeDimacs(double* dimacs_error,
			    SolveInfo& solveInfo,
			    Residuals& currentRes,
			    Solutions& currentPt,
			    InputData& inputData,
			    WorkVariables& work);
  
  static void printDimacs(double* dimacs_error,
			  char* printFormat,
			  FILE* fpout);

  static void printSolution(BlockStruct& bs, Solutions& currentPt,
			    Parameter& param, FILE* fpout);

};

} // end of namespace 'sdpa'

#endif // __sdpa_io_h__
