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

#include <cstdio>
#include <cstdlib>
#include <sdpa_call.h>
using namespace sdpa;

#define USER_PARAMETER_FILE ((char *)"./param.sdpa")
#define DEFAULT_PARAMETER_FILE ((char *)"/usr/share/sdpa/param.sdpa")
// PARAMETER_FILE is decided by the following priority
// 1: The file assigned by '-p' option of 'option type 2'.
//    For 'option type1', this is skipped.
// 2: USER_PARAMETER_FILE
//    For 'option type2', this is skipped.
// 3: DEFAULT_PARAMETER_FILE
// 4: Default parameter
static void message(char* argv0)
{
  cout << endl;
  cout << "*** Please assign data file and output file.***" << endl;
  cout << endl;
  cout << "---- option type 1 ------------" << endl;
  cout << argv0 <<" DataFile OutputFile [InitialPtFile]"
    " [-pt parameters] [-dimacs] [-numThreads numThreads]"<< endl;
  cout << "parameters = 0 default, 1 fast (unstable),"
    " 2 slow (stable)" << endl;
  cout << "  -dimacs : printout dimacs information incurring additional computation cost " << endl;
  cout << "  -numThreads: Number of pthreads for internal computation" << endl;
  cout << "example1-1: " << argv0
       << " example1.dat example1.result" << endl;
  cout << "example1-2: " << argv0
       << " example1.dat-s example1.result" << endl;
  cout << "example1-3: " << argv0
       << " example1.dat example1.result example1.ini" << endl;
  cout << "example1-4: " << argv0
       << " example1.dat example1.result -pt 2" << endl;
  cout << "example1-5: " << argv0
       << " example1.dat example1.result -dimacs" << endl;
  cout << "example1-6: " << argv0
       << " example1.dat example1.result -numThreads 4" << endl;

  cout << endl;
  cout << "---- option type 2 ------------" << endl;
  cout << argv0 << " [option filename]+ " << endl;
  cout << "  -dd : data dense :: -ds : data sparse     " << endl;
  cout << "  -id : init dense :: -is : init sparse     " << endl;
  cout << "  -o  : output     :: -p  : parameter       " << endl;
  cout << "  -pt : parameters , 0 default, 1 fast (unstable)" << endl;
  cout << "                     2 slow (stable)         " << endl;
  cout << "  -dimacs : printout dimacs information incurring additional computation cost " << endl;
  cout << "  -numThreads: Number of pthreads for internal computation" << endl;
  cout << "example2-1: " << argv0
       << " -o example1.result -dd example1.dat" << endl;
  cout << "example2-2: " << argv0
       << " -ds example1.dat-s -o example1.result "
       << "-p param.sdpa" << endl;
  cout << "example2-3: " << argv0
       << " -ds example1.dat-s -o example1.result "
       << "-pt 2" << endl;
  cout << "example2-4: " << argv0
       << " -ds example1.dat-s -o example1.result "
       << "-dimacs" << endl;
  cout << "example2-5: " << argv0
       << " -ds example1.dat-s -o example1.result "
       << "-numThreads 4" << endl;

  cout << endl;
  cout << "---- option type 3 ------------" << endl;
  cout << argv0 << " --version " << endl;
  cout << "   to print out version and exit." << endl;

  cout << endl << endl;
  cout << "PARAMETER_FILE is decided by the following priority" << endl;
  cout << "   1: The file assigned by '-p' option of 'option type 2'." << endl;
  cout << "       For 'option type1', this is skipped." << endl;
  cout << "   2: " << USER_PARAMETER_FILE  << endl;
  cout << "       For 'option type2', this is skipped." << endl;
  cout << "   3: " << DEFAULT_PARAMETER_FILE << endl;
  cout << "   4: Default parameter" << endl;
  exit(1);
}

static void argumentAnalysis(SDPA& Problem1,
			     int argc, char** argv,
			     char*& inputFileName,
			     char*& resultFileName,
			     char*& initFileName,
			     char*& paramFileName,
			     SDPA::SparseType& isInputSparse,
			     SDPA::SparseType& isInitSparse,
			     SDPA::ParameterType& parameterType,
			     bool& isDimacs, int& numThreads)
{
  if (argc == 1) {
    message(argv[0]);
  }
  if (strcmp(argv[1],"--version") == 0) {
    fprintf(stdout,"====\n");
    fprintf(stdout,"SDPA (SemiDefinite Programming Algorithm) %s\n",sdpa_version);
    fprintf(stdout,"     %s\n",sdpa_right);
    fprintf(stdout,"====\n");
    exit(0);
  }
  if (argv[1][0] == '-') {
    // rsdpa argument
    
    for (int index = 0; index < argc; ++index) {
      char* target = argv[index];
      if (strcmp(target,"-dd")==0 && index+1 < argc) {
	inputFileName = argv[index+1];
	isInputSparse = SDPA::DENSE;
	index++;
	continue;
      }
      if (strcmp(target,"-ds")==0 && index+1 < argc) {
	inputFileName = argv[index+1];
	isInputSparse = SDPA::SPARSE;
	continue;
      }
      if (strcmp(target,"-id")==0 && index+1 < argc) {
	initFileName = argv[index+1];
	isInitSparse = SDPA::DENSE;
	continue;
      }
      if (strcmp(target,"-is")==0 && index+1 < argc) {
	initFileName = argv[index+1];
	isInitSparse = SDPA::SPARSE;
	index++;
	continue;
      }
      if (strcmp(target,"-o")==0 && index+1 < argc) {
	resultFileName = argv[index+1];
	index++;
	continue;
      }
      if (strcmp(target,"-p")==0 && index+1 < argc) {
	paramFileName = argv[index+1];
	index++;
	continue;
      }
      if (strcmp(target,"-k")==0 && index+1 < argc) {
	double KAPPA = atof(argv[index+1]);
	Problem1.setKappa(KAPPA);
	rMessage("Kappa = " << KAPPA);
	index++;
	continue;
      }
      if (strcmp(target,"-dimacs")==0) {
	isDimacs = true;
	continue;
      }
      if (strcmp(target,"-pt")==0 && index+1 < argc) {
	int tmp = atoi(argv[index+1]);
	switch (tmp) {
	case 0:
	  parameterType = SDPA::PARAMETER_DEFAULT;
	  break;
	case 1:
	  parameterType = SDPA::PARAMETER_UNSTABLE_BUT_FAST;
	  break;
	case 2:
	  parameterType = SDPA::PARAMETER_STABLE_BUT_SLOW;
	  break;
	default:
	  parameterType = SDPA::PARAMETER_DEFAULT;
	}
	index++;
	paramFileName = NULL;
	continue;
      }
      if (strcmp(target,"-numThreads")==0 && index+1 < argc) {
	numThreads = atoi(argv[index+1]);
	index++;
	continue;
      }
    }
  }
  else { // SDPA argument
    inputFileName = argv[1];
    int len = strlen(inputFileName);
    if (inputFileName[len-1] == 's'
	&& inputFileName[len-2] == '-') {
      isInputSparse = SDPA::SPARSE;
    }
	
    resultFileName = argv[2];

    paramFileName = USER_PARAMETER_FILE;

    for (int index=3; index<argc; ++index) {
      if (strcmp(argv[index],"-dimacs")==0) {
	isDimacs = true;
      }
      else if (strcmp(argv[index],"-numThreads")==0 && index+1 < argc) {
	numThreads = atoi(argv[index+1]);
	++index;
      }
      else if (strcmp(argv[index],"-pt")==0 && index+1 < argc) {
	int tmp = atoi(argv[index+1]);
	switch (tmp) {
	case 0:
	  parameterType = SDPA::PARAMETER_DEFAULT;
	  break;
	case 1:
	  parameterType = SDPA::PARAMETER_UNSTABLE_BUT_FAST;
	  break;
	case 2:
	  parameterType = SDPA::PARAMETER_STABLE_BUT_SLOW;
	  break;
	default:
	  parameterType = SDPA::PARAMETER_DEFAULT;
	}
	index++;
	paramFileName = NULL;
      } // end of "-pt"
      else {
	initFileName = argv[index];
	int len = strlen(initFileName);
	if (initFileName[len-1] == 's'
	    && initFileName[len-2] == '-') {
	  isInitSparse = SDPA::SPARSE;
	}
      }
    } // end of 'for'
  }

  if (paramFileName != NULL) {
    // check the availability for paramFileName,
    // usually USER_PARAMETER_FILE
    FILE* fptmp = NULL;
    if ((fptmp=fopen(paramFileName,"r"))==NULL) {
      // we try DEFAULT_PARAMETER_FILE
      if ((fptmp=fopen(DEFAULT_PARAMETER_FILE,"r"))==NULL) {
	rMessage("Cannot Open user parameter File " << paramFileName
	       << " and default parameter file " << DEFAULT_PARAMETER_FILE);
	rMessage("Default parameter will be used.");
	paramFileName = NULL;
      }
      else {
	paramFileName = DEFAULT_PARAMETER_FILE;
      }
    }
  } // end of 'if (paramFileName != NULL)'

  
  if (inputFileName == NULL || resultFileName == NULL) {
    message(argv[0]);
  }
}

 
int main(int argc, char** argv)
{
  TimeStart(ALL_START1);
  SDPA Problem1;
  time_t ltime;
  time(&ltime);
  char string_time[1024];
  strcpy(string_time,ctime(&ltime));
  string_time[strlen(string_time)-1]='\0';
  fprintf(stdout,"SDPA (Version %s) start at [%s]\n", sdpa_version, string_time);
  // cout << "let me see your ..." << endl;
  if (argc == 1) {
    message(argv[0]);
  }
  char* inputFileName  = NULL;
  char* resultFileName = NULL;
  char* initFileName   = NULL;
  char* paramFileName  = NULL;

  SDPA::SparseType isInputSparse = SDPA::DENSE;
  SDPA::SparseType isInitSparse  = SDPA::DENSE;
  SDPA::ParameterType parameterType = SDPA::PARAMETER_DEFAULT;

  bool isDimacs = false;
  int  numThreads = 0; // 0 means automatic computation

  argumentAnalysis(Problem1, argc, argv,
		   inputFileName, resultFileName, initFileName,
		   paramFileName, isInputSparse, isInitSparse,
		   parameterType, isDimacs, numThreads);
  Problem1.setDisplay(stdout);
  FILE* fpresult;
  if ((fpresult = fopen(resultFileName,"w")) == NULL) {
    rError("Cannot Open Result File : " << resultFileName);
  }
  fprintf(fpresult,"SDPA start at [%s]\n",string_time);

  Problem1.setResultFile(fpresult);
  if (paramFileName == NULL) {
    if (parameterType == SDPA::PARAMETER_DEFAULT) {
      fprintf(stdout  ,"set   is DEFAULT\n");
      fprintf(fpresult,"set   is DEFAULT\n");
    }
    else if (parameterType == SDPA::PARAMETER_UNSTABLE_BUT_FAST) {
      fprintf(stdout  ,"set   is UNSTABLE_BUT_FAST\n");
      fprintf(fpresult,"set   is UNSTABLE_BUT_FAST\n");
    }
    else if (parameterType == SDPA::PARAMETER_STABLE_BUT_SLOW) {
      fprintf(stdout  ,"set   is STABLE_BUT_SLOW\n");
      fprintf(fpresult,"set   is STABLE_BUT_SLOW\n");
    }
    Problem1.setParameterType(parameterType);
  }

  if (paramFileName) {
    fprintf(stdout  ,"param is %s\n", paramFileName);
    Problem1.readParameter(paramFileName,fpresult);
  }
  
  fprintf(stdout  ,"data  is %s", inputFileName);
  if (isInputSparse == SDPA::SPARSE) {
    fprintf(stdout  ," : sparse\n");
  }
  else {
    fprintf(stdout  ," : dense\n");
  }
  Problem1.readInput(inputFileName, fpresult, isInputSparse);
  
  if (initFileName) {
    Problem1.setInitPoint(true);
    fprintf(stdout  ,"init  is %s", initFileName);
    if (isInitSparse == SDPA::SPARSE) {
      fprintf(stdout  ," : sparse\n");
    }
    else {
      fprintf(stdout  ," : dense\n");
    }
    Problem1.readInit(initFileName, fpresult, isInitSparse);
  }
  fprintf(stdout  ,"out   is %s\n", resultFileName);
  fprintf(fpresult,"out    is %s\n", resultFileName);
  if (isDimacs) {
    fprintf(stdout  ,"Dimacs information will be computed after the iteration.\n");
    fprintf(fpresult,"Dimacs information will be computed after the iteration.\n");
  }

  Problem1.setNumThreads(numThreads);
  Problem1.initializeSolve();
  Problem1.solve();

  if (isDimacs) {
    double dimacs_error[7];
    Problem1.getDimacsError(dimacs_error);
    Problem1.printDimacsError(dimacs_error,
			      Problem1.getParameterPrintInformation(),
			      stdout);
    Problem1.printDimacsError(dimacs_error,
			      Problem1.getParameterPrintInformation(),
			      fpresult);
  }
  
  Problem1.terminate();
  
  time(&ltime);
  strcpy(string_time,ctime(&ltime));
  string_time[strlen(string_time)-1]='\0';
  fprintf(stdout  ,"SDPA end at [%s]\n",string_time);
  fprintf(fpresult,"SDPA end at [%s]\n",string_time);
  TimeEnd(ALL_END1);
  double all_time = TimeCal(ALL_START1,ALL_END1);
  fprintf(stdout  ,"ALL TIME = %.6lf\n", all_time);
  fprintf(fpresult,"ALL TIME = %.6lf\n", all_time);
  fclose(fpresult);
  
  return 0;
}

