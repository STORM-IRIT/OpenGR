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

#include "sdpa_parts.h"
#include "sdpa_jordan.h"
#include "sdpa_linear.h"
#include "sdpa_newton.h"  

namespace sdpa {

ComputeTime::ComputeTime()
{
  Predictor = 0.0;
  Corrector = 0.0;
  
  StepPredictor = 0.0;
  StepCorrector = 0.0;
  
  xMatTime = 0.0;
  zMatTime = 0.0;
  xMatzMatTime = 0.0;
  invzMatTime = 0.0;
  
  EigxMatTime = 0.0;
  EigzMatTime = 0.0;
  EigxMatzMatTime = 0.0;

  makebMat = 0.0;
  B_DIAG   = 0.0;
  B_F1     = 0.0;
  B_F2     = 0.0;
  B_F3     = 0.0;
  B_PRE    = 0.0;

  makegVecMul = 0.0;
  makegVec = 0.0;
  makerMat = 0.0;
  choleskybMat = 0.0;
  
  solve    = 0.0;
  sumDz    = 0.0;
  makedX    = 0.0;
  symmetriseDx = 0.0;
  makedXdZ = 0.0;
  updateRes = 0.0;

  MainLoop  = 0.0;
  FileRead  = 0.0;
  FileCheck = 0.0;
  FileChange= 0.0;
  TotalTime = 0.0;
}

ComputeTime::~ComputeTime()
{
  // Nothing needs.
}


void ComputeTime::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  fprintf(fpout,"\n");
  #if 0
  if (TotalTime <= 0.0) {
    return;
  }
  #endif
  fprintf(fpout, "                         Time(sec) ");
  fprintf(fpout," Ratio(%% : MainLoop) \n");
  fprintf(fpout, " Predictor time  =       %f,  %f\n",
	  Predictor, Predictor/MainLoop*100.0);
  fprintf(fpout, " Corrector time  =       %f,  %f\n",
	  Corrector, Corrector/MainLoop*100.0);
  fprintf(fpout, " Make bMat time  =       %f,  %f\n",
	  makebMat, makebMat/MainLoop*100.0);
  fprintf(fpout, " Make bDia time  =       %f,  %f\n",
	  B_DIAG,B_DIAG/MainLoop*100.0);
  fprintf(fpout, " Make bF1  time  =       %f,  %f\n",
	  B_F1,B_F1/MainLoop*100.0);
  fprintf(fpout, " Make bF2  time  =       %f,  %f\n",
	  B_F2,B_F2/MainLoop*100.0);
  fprintf(fpout, " Make bF3  time  =       %f,  %f\n",
	  B_F3,B_F3/MainLoop*100.0);
  fprintf(fpout, " Make bPRE time  =       %f,  %f\n",
	  B_PRE,B_PRE/MainLoop*100.0);
  fprintf(fpout, " Make rMat time  =       %f,  %f\n",
	  makerMat, makerMat/MainLoop*100.0);
  fprintf(fpout, " Make gVec Mul   =       %f,  %f\n",
	  makegVecMul, makegVecMul/MainLoop*100.0);
  fprintf(fpout, " Make gVec time  =       %f,  %f\n",
	  makegVec, makegVec/MainLoop*100.0);
  fprintf(fpout, " Cholesky bMat   =       %f,  %f\n",
	  choleskybMat, choleskybMat/MainLoop*100.0);
  fprintf(fpout, " Ste Pre time    =       %f,  %f\n",
	  StepPredictor, StepPredictor/MainLoop*100.0);
  fprintf(fpout, " Ste Cor time    =       %f,  %f\n",
	  StepCorrector, StepCorrector/MainLoop*100.0);
  fprintf(fpout, " solve           =       %f,  %f\n",
	  solve, solve/MainLoop*100.0);
  fprintf(fpout, " sumDz           =       %f,  %f\n",
  	  sumDz, sumDz/MainLoop*100.0);
  fprintf(fpout, " makedX          =       %f,  %f\n",
  	  makedX, makedX/MainLoop*100.0);
  fprintf(fpout, " symmetriseDx    =       %f,  %f\n",
  	  symmetriseDx, symmetriseDx/MainLoop*100.0);
  fprintf(fpout, " makedXdZ        =       %f,  %f\n",
	  makedXdZ, makedXdZ/MainLoop*100.0);
  fprintf(fpout, " xMatTime        =       %f,  %f\n",
	  xMatTime, xMatTime/MainLoop*100.0);
  fprintf(fpout, " zMatTime        =       %f,  %f\n",
	  zMatTime, zMatTime/MainLoop*100.0);
  fprintf(fpout, " invzMatTime     =       %f,  %f\n",
  	  invzMatTime, invzMatTime/MainLoop*100.0);
  fprintf(fpout, " xMatzMatTime    =       %f,  %f\n",
	  xMatzMatTime, xMatzMatTime/MainLoop*100.0);
  fprintf(fpout, " EigxMatTime     =       %f,  %f\n",
	  EigxMatTime, EigxMatTime/MainLoop*100.0);
  fprintf(fpout, " EigzMatTime     =       %f,  %f\n",
	  EigzMatTime, EigzMatTime/MainLoop*100.0);
  fprintf(fpout, " EigxMatzMatTime =       %f,  %f\n",
	  EigxMatzMatTime, EigxMatzMatTime/MainLoop*100.0);
  fprintf(fpout, " updateRes       =       %f,  %f\n",
	  updateRes, updateRes/MainLoop*100.0);
  double total_eigen = EigxMatTime + EigzMatTime + EigxMatzMatTime;
  fprintf(fpout, " EigTime         =       %f,  %f\n",
	  total_eigen, total_eigen/MainLoop*100.0);
  double sub_total_bMat = MainLoop - makebMat;
  fprintf(fpout, " sub_total_bMat  =       %f,  %f\n",
	  sub_total_bMat, sub_total_bMat/MainLoop*100.0);
  fprintf(fpout, " Main Loop       =       %f,  %f\n",
	  MainLoop, MainLoop/MainLoop*100.0);
  fprintf(fpout, " File Check      =       %f,  %f\n",
	  FileCheck, FileCheck/MainLoop*100.0);
  fprintf(fpout, " File Change     =       %f,  %f\n",
	  FileChange, FileChange/MainLoop*100.0);
  fprintf(fpout, " File Read       =       %f,  %f\n",
	  FileRead, FileRead/MainLoop*100.0);
  fprintf(fpout, " Total           =       %f,  %f\n",
	  TotalTime, TotalTime/MainLoop*100.0);
  fprintf(fpout, "\n");

  return;
}

//-------------------------------------------------------------

Parameter::Parameter()
{
  // setDefaultParameter();
}
Parameter::Parameter(FILE* parameterFile)
{
  readFile(parameterFile);
}

Parameter::~Parameter()
{
  // Nothings needs.
}

void Parameter::setDefaultParameter(Parameter::parameterType type)
{
  if (type == PARAMETER_STABLE_BUT_SLOW) {
    maxIteration =  1000;
    epsilonStar  =  1.0e-7;
    lambdaStar   =  1.0e+4;
    omegaStar    =  2.0;
    lowerBound   = -1.0e+5;
    upperBound   =  1.0e+5;
    betaStar     =  0.1;
    betaBar      =  0.3;
    gammaStar    =  0.8;
    epsilonDash  =  1.0e-7;
  }
  else if (type == PARAMETER_UNSTABLE_BUT_FAST) {
    maxIteration =  100;
    epsilonStar  =  1.0e-7;
    lambdaStar   =  1.0e+2;
    omegaStar    =  2.0;
    lowerBound   = -1.0e+5;
    upperBound   =  1.0e+5;
    betaStar     =  0.01;
    betaBar      =  0.02;
    gammaStar    =  0.95;
    epsilonDash  =  1.0e-7;
  }
  else {
    maxIteration =  100;
    epsilonStar  =  1.0e-7;
    lambdaStar   =  1.0e+2;
    omegaStar    =  2.0;
    lowerBound   = -1.0e+5;
    upperBound   =  1.0e+5;
    betaStar     =  0.1;
    betaBar      =  0.2;
    gammaStar    =  0.9;
    epsilonDash  =  1.0e-7;
  }    
  strcpy(xPrint,xPRINT_DEFAULT);
  strcpy(XPrint,XPRINT_DEFAULT);
  strcpy(YPrint,YPRINT_DEFAULT);
  strcpy(infPrint,infPRINT_DEFAULT);
}

char Parameter::xPRINT_DEFAULT[PRINT_DEFAULT_LENGTH] = "%+8.3e";
char Parameter::XPRINT_DEFAULT[PRINT_DEFAULT_LENGTH] = "%+8.3e";
char Parameter::YPRINT_DEFAULT[PRINT_DEFAULT_LENGTH] = "%+8.3e";
char Parameter::infPRINT_DEFAULT[PRINT_DEFAULT_LENGTH] = "%+10.16e";

void Parameter::readFile(FILE* parameterFile)
{
  fscanf(parameterFile,"%d%*[^\n]",&maxIteration);
  fscanf(parameterFile,"%lf%*[^\n]",&epsilonStar);
  fscanf(parameterFile,"%lf%*[^\n]",&lambdaStar);
  fscanf(parameterFile,"%lf%*[^\n]",&omegaStar);
  fscanf(parameterFile,"%lf%*[^\n]",&lowerBound);
  fscanf(parameterFile,"%lf%*[^\n]",&upperBound);
  fscanf(parameterFile,"%lf%*[^\n]",&betaStar);
  fscanf(parameterFile,"%lf%*[^\n]",&betaBar);
  fscanf(parameterFile,"%lf%*[^\n]",&gammaStar);
  fscanf(parameterFile,"%lf%*[^\n]",&epsilonDash);
  fscanf(parameterFile,"%s %*[^\n]",xPrint);
  fscanf(parameterFile,"%s %*[^\n]",XPrint);
  fscanf(parameterFile,"%s %*[^\n]",YPrint);
  fscanf(parameterFile,"%s %*[^\n]",infPrint);
  if (strcmp(xPrint,NO_P_FORMAT)!=0 && xPrint[0]!='%') {
    rMessage("Strange xPrint[" << xPrint << "]"
	     " migh cause trouble when printing x");
  }
  if (strcmp(XPrint,NO_P_FORMAT)!=0 && XPrint[0]!='%') {
    rMessage("Strange XPrint[" << XPrint << "]"
	     " migh cause trouble when printing X.");
  }
  if (strcmp(YPrint,NO_P_FORMAT)!=0 && YPrint[0]!='%') {
    rMessage("Strange YPrint[" << YPrint << "]"
	     " migh cause trouble when printing Y.");
  }
  if (strcmp(infPrint,NO_P_FORMAT)!=0 && infPrint[0]!='%') {
    rMessage("Strange infPrint[" << infPrint << "]"
	     " migh cause trouble when printing information.");
  }
}


void Parameter::display(FILE* fpout, char* printFormat)
{
  if (fpout == NULL) {
    return;
  }
  if (strcmp(printFormat,NO_P_FORMAT) == 0) {
    fprintf(fpout,"%s\n",NO_P_FORMAT);
    return;
  }
  fprintf(fpout, "** Parameters **\n");
  fprintf(fpout, "maxIteration = %d\n",maxIteration);
  fprintf(fpout, "epsilonStar  = ");
  fprintf(fpout, printFormat, epsilonStar );
  fprintf(fpout, "\n");
  fprintf(fpout, "lambdaStar   = ");
  fprintf(fpout, printFormat, lambdaStar );
  fprintf(fpout, "\n");
  fprintf(fpout, "omegaStar    = ");
  fprintf(fpout, printFormat, omegaStar  );
  fprintf(fpout, "\n");
  fprintf(fpout, "lowerBound   = ");
  fprintf(fpout, printFormat, lowerBound);
  fprintf(fpout, "\n");
  fprintf(fpout, "upperBound   = ");
  fprintf(fpout, printFormat, upperBound);
  fprintf(fpout, "\n");
  fprintf(fpout, "betaStar     = ");
  fprintf(fpout, printFormat, betaStar   );
  fprintf(fpout, "\n");
  fprintf(fpout, "betaBar      = ");
  fprintf(fpout, printFormat, betaBar    );
  fprintf(fpout, "\n");
  fprintf(fpout, "gammaStar    = ");
  fprintf(fpout, printFormat, gammaStar  );
  fprintf(fpout, "\n");
  fprintf(fpout, "epsilonDash  = ");
  fprintf(fpout, printFormat, epsilonDash );
  fprintf(fpout, "\n");
  #if 1
  fprintf(fpout, "xPrint       = %s \n", xPrint );
  fprintf(fpout, "XPrint       = %s \n", XPrint );
  fprintf(fpout, "YPrint       = %s \n", YPrint );
  fprintf(fpout, "infPrint     = %s \n", infPrint );
  #endif
  return;
}

//----------------------------------------------------------

StepLength::StepLength()
{
  primal = 0.0;
  dual   = 0.0;

}

StepLength::~StepLength()
{
  terminate();
}


void StepLength::initialize(double alphaP, double alphaD)
{
  primal = alphaP;
  dual   = alphaD;
}

void StepLength::terminate()
{
  // Nothing needs.
}

double StepLength::minBlockVector(BlockVector& aVec)
{
  int nBlock = aVec.nBlock;
  double ret = aVec.ele[0].ele[0];
  double tmp;
  int size = aVec.ele[0].nDim;
  for (int j=1; j<size; ++j) {
    tmp = aVec.ele[0].ele[j];
    if (tmp < ret) {
      ret = tmp;
    }
  }
  for (int k=1; k<nBlock; ++k) {
    size = aVec.ele[k].nDim;
    for (int j=0; j<size; ++j) {
      tmp = aVec.ele[k].ele[j];
      if (tmp < ret) {
	ret = tmp;
      }
    }
  }
  return ret;
}

void StepLength::computeStepLength(Solutions& currentPt,
				   Newton& newton,
				   WorkVariables& work,
				   ComputeTime& com)
{
  double alphaBD = 100.0;

  // calculate  eigenvalues of X^{-1} dX
  TimeStart(START1);
  // rMessage("invCholeskyX=");
  // currentPt.invCholeskyX.display();
  // rMessage("Dx=");
  // newton.DxMat.display();
  double minxInvDxEigenValue;
  #define ALL_EIGEN 0
  #if ALL_EIGEN
  Lal::let(work.DLS2,'=',newton.DxMat,'T',currentPt.invCholeskyX);
  Lal::let(work.DLS1,'=',currentPt.invCholeskyX,'*',work.DLS2);
  // rMessage("LDxLT=");
  // workMat1.display();
  // rMessage("eigen test");
  Lal::getMinEigenValue(work.DLS1,xInvDxEigenValues,workVec);
  minxInvDxEigenValue = minBlockVector(xInvDxEigenValues);
  // rMessage("minxInvDxEigenValue = " << minxInvDxEigenValue);
  #else
  minxInvDxEigenValue
    = Jal::getMinEigen(currentPt.invCholeskyX,newton.DxMat,work);
  // rMessage("minxInvDxEigenValue = " << minxInvDxEigenValue);
  #endif
  if (-minxInvDxEigenValue > 1.0 /alphaBD) {
    primal = - 1.0/minxInvDxEigenValue;
    // the limint of primal steplength
  } else {
    primal = alphaBD;
  }
  TimeEnd(END1);
  com.EigxMatTime += TimeCal(START1,END1);

  
  // calculate  eigenvalues of Z^{-1} dZ
  TimeStart(START2);
  // rMessage("invCholeskyZ=");
  // currentPt.invCholeskyZ.display();
  // rMessage("Dz=");
  // newton.DzMat.display();

  double minzInvDzEigenValue;
  #if ALL_EIGEN
  Lal::let(work.DLS2,'=',newton.DzMat,'T',currentPt.invCholeskyZ);
  Lal::let(work.DLS1,'=',currentPt.invCholeskyZ,'*',work.DLS2);
  // rMessage("LDzLT=");
  // workMat1.display();
  // rMessage("eigen test");
  Lal::getMinEigenValue(work.DLS1,zInvDzEigenValues,workVec);
  minzInvDzEigenValue = minBlockVector(zInvDzEigenValues);
  // rMessage("minzInvDzEigenValue = " << minzInvDzEigenValue);
  #else
  minzInvDzEigenValue
    = Jal::getMinEigen(currentPt.invCholeskyZ,newton.DzMat,work);
  // rMessage("minzInvDzEigenValue = " << minzInvDzEigenValue);
  #endif
  if (-minzInvDzEigenValue > 1.0 /alphaBD) {
    dual = - 1.0/minzInvDzEigenValue;
    // the limint of dual steplength
  } else {
    dual = alphaBD;
  }
  TimeEnd(END2);
  com.EigzMatTime += TimeCal(START2,END2);
}


void StepLength::MehrotraPredictor(InputData& inputData,
				   Solutions& currentPt,
				   Phase& phase,
				   Newton& newton,
				   WorkVariables& work,
				   ComputeTime& com)
{
  #if 1
  primal = dual = 0.9;
  #else
  computeStepLength(currentPt, newton, work, com);
  #endif

  Vector& b = inputData.b;
  SparseLinearSpace& C = inputData.C;

  if (phase.value==SolveInfo::noINFO
      || phase.value==SolveInfo::dFEAS) {
    // primal is infeasible
    if (primal>1.0) {
      primal = 1.0;
    }
  } else {
    // when primal is feasible,
    // check stepP1 is effective or not.
    double incPrimalObj;
    Lal::let(incPrimalObj,'=',C,'.',newton.DxMat);
    if (incPrimalObj>0.0) {
      if (primal>dual) {
	primal = dual;
      }
      if (primal>1.0) {
	primal = 1.0;
      }
    }
  }
  if (phase.value==SolveInfo::noINFO
      || phase.value==SolveInfo::pFEAS) {
    // dual is infeasible
    if (dual>1.0) {
      dual = 1.0;
    }
  } else {
    // when dual is feasible
    // check stepD1 is effective or not.
    double incDualObj;
    Lal::let(incDualObj,'=',b,'.',newton.DyVec);
    if(incDualObj<0.0) {
      if (dual>primal) {
	dual = primal;
      }
      if (dual>1.0) {
	dual = 1.0;
      }
    }
  }
}

void StepLength::MehrotraCorrector(InputData& inputData,
				   Solutions& currentPt,
				   Phase& phase,
				   Switch& reduction,
				   Newton& newton,
				   AverageComplementarity& mu,
				   RatioInitResCurrentRes& theta,
				   WorkVariables& work,
				   Parameter& param,
				   ComputeTime& com)
{
  double xi      = 3.0;
  
  Vector& b = inputData.b;
  SparseLinearSpace& C = inputData.C;
  int nDim = currentPt.nDim;

  computeStepLength(currentPt, newton, work, com);

  // adjust steplength with param.gammaStar
  // param.gammaStar = 0.5;
  primal = param.gammaStar * primal;
  dual   = param.gammaStar * dual;
  // rMessage("primal = " << primal);
  // rMessage("dual = " << dual);
  // rMessage("phase = ");
  // phase.display();
  
  if (phase.value==SolveInfo::noINFO
      || phase.value==SolveInfo::dFEAS) {
    // primal is infeasible.
    if (primal>1.0) {
      primal = 1.0;
    }
  } else {
    double incPrimalObj;
    Lal::let(incPrimalObj,'=',C,'.',newton.DxMat);
    if(incPrimalObj>0.0) {
      // when primal is feasible
      // check stepD1 is effective or not.
      if (primal>dual) {
	primal = dual;
      }
      if (primal>1.0) {
	primal = 1.0;
      }
    }
  }
  if (phase.value==SolveInfo::noINFO
      || phase.value==SolveInfo::pFEAS) {
    // dual is infeasible
    if (dual>1.0) {
      dual = 1.0;
    }
  } else {
    // when dual is feasible
    // check stepD1 is effective or not.
    double incDualObj;
    Lal::let(incDualObj,'=',b,'.',newton.DyVec);
    if(incDualObj<0.0) {
      if (dual>primal) {
	// change because noneffective
	dual = primal;
      }
      if (dual>1.0) {
	dual = 1.0;
      }
    }
  }
#if 1
  // attain feasibility before mu reduction
  if (reduction.switchType==Switch::CENTERING
      && (phase.value == SolveInfo::noINFO
	  || phase.value == SolveInfo::pFEAS
	  || phase.value == SolveInfo::dFEAS) ) {
    double xMatvMat;
    Lal::let(xMatvMat,'=',currentPt.xMat,'.',newton.DzMat);
    double uMatzMat;
    Lal::let(uMatzMat,'=',newton.DxMat,'.',currentPt.zMat);
    double uMatvMat;
    Lal::let(uMatvMat,'=',newton.DxMat,'.',newton.DzMat);

    double thetaMax = max((1.0-primal)*theta.primal,
			  (1.0-dual  )*theta.dual);
    double muNew = mu.current
      + (primal*uMatzMat + dual*xMatvMat
	 + primal*dual*uMatvMat) / nDim;
    double alphaMax;
    //   memo by kazuhide nakata
    //   thetaMax*mu.initial -> thetamax*thetaMax*mu.initial ???
	//    while (thetaMax*mu.initial > xi*muNew) {
	while (thetaMax*thetaMax*mu.initial > xi*muNew) {
      alphaMax = 0.95 * max(primal,dual);
      primal = min(primal,alphaMax);
      dual   = min(dual  ,alphaMax);
      thetaMax = max((1.0-primal)*theta.primal,
		     (1.0-dual  )*theta.dual);
      muNew = mu.current + (primal*uMatzMat + dual*xMatvMat
			    + primal*dual*uMatvMat) / nDim;
      // if "too short step", then break down the algorithm.
      if (primal < 1.0e-6 && dual < 1.0e-6) {
	break;
      }
    }
  }
#endif

  // 2007/09/18 kazuhide nakata
  // keep primal objective value >= dual objective value
  if (phase.value == SolveInfo::pdFEAS){
	// if (mu.current < 1.0){
	
	double objValDual,objValPrimal,incDualObj,incPrimalObj,maxRatio;

	Lal::let(objValDual,'=',inputData.b,'.',currentPt.yVec);
	Lal::let(objValPrimal,'=',inputData.C,'.',currentPt.xMat);
	Lal::let(incDualObj,'=',b,'.',newton.DyVec);
	incDualObj *= dual;
	Lal::let(incPrimalObj,'=',C,'.',newton.DxMat);
	incPrimalObj *= primal;
	maxRatio = (objValDual - objValPrimal)
	  / (incPrimalObj - incDualObj);

	if ((maxRatio > 0.0) && (maxRatio < 1.0)){
	  primal *= maxRatio;
	  dual *= maxRatio;
#if 0
	  printf("max stepsise ratio: %9.1e\n",maxRatio);
	  printf("new stepsize  primal:%9.1e, dual:%9.1Fe\n",primal,dual);
#endif
	}
  }
}

void StepLength::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }

  fprintf(fpout,"alpha.primal = %8.3e\n",primal);
  fprintf(fpout,"alpha.dual   = %8.3e\n",dual);
}

//-------------------------------------------------
DirectionParameter::DirectionParameter(double betaStar)
{
  initialize(betaStar);
}

DirectionParameter::~DirectionParameter()
{
  // Nothing needs.
}

void DirectionParameter::initialize(double betaStar)
{
  value = betaStar;
}

void DirectionParameter::MehrotraPredictor(Phase& phase,
					   Switch& reduction,
					   Parameter& param)
{
  const double nu = 2.0;
  if (phase.value == SolveInfo::pdFEAS) {
    value = 0.0;
  } else {
    value = param.betaBar;
    if (reduction.switchType==Switch::AFFINE) {
      value = nu;
    }
  }
}

void DirectionParameter::
MehrotraCorrector(Phase& phase,StepLength& alpha,
		  Solutions& currentPt,Newton& newton,
		  AverageComplementarity& mu,Parameter& param)
{
  int nDim = currentPt.nDim;

  double xMatvMat;
  Lal::let(xMatvMat,'=',currentPt.xMat,'.',newton.DzMat);
  double uMatzMat;
  Lal::let(uMatzMat,'=',newton.DxMat,'.',currentPt.zMat);
  double uMatvMat;
  Lal::let(uMatvMat,'=',newton.DxMat,'.',newton.DzMat);

  double muTarget = mu.current
    + (alpha.primal*uMatzMat + alpha.dual*xMatvMat
       + alpha.primal*alpha.dual*uMatvMat) / nDim;
  // rMessage("muTarget  : " << muTarget);
  // rMessage("muCurrent : " << mu.current);
  value = muTarget/mu.current;
  // rMessage("muValue : " << value);
  if (value < 1.0) {
    value = value*value;
  }
  if (phase.value==SolveInfo::pdFEAS) {
    // rMessage("MehrotraCorrector : pdFEAS" << value);
    if (value < param.betaStar) {
      value = param.betaStar;
    }
    if (value > 1.0) {
      value = 1.0;
    }
  } else {
    if (value < param.betaBar) {
      value = param.betaBar;
    }
  }
  // rMessage("MehrotraCorrector : " << value);
}

void DirectionParameter::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  fprintf(fpout,"beta.value = %8.3e\n",value);
}

//---------------------------------------------------

Switch::Switch(SwitchType switchType)
{
  initialize(switchType);
}

Switch::~Switch()
{
  // Nothing needs.
}

void Switch::initialize(SwitchType switchType)
{
  this->switchType = switchType;
}

void Switch::MehrotraPredictor(Phase& phase)
{
  if (phase.value==SolveInfo::noINFO
      || phase.value==SolveInfo::pFEAS
      || phase.value==SolveInfo::dFEAS) {
    // At least one of primal or dual is infeasible.
    switchType = CENTERING;
  } else {
    switchType = AFFINE;
  }
}

void Switch::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }

  if (switchType == CENTERING) {
    fprintf(fpout,"reduction.switchType == CENTERING\n");
  } else {
    fprintf(fpout,"reduction.switchType == AFFINE\n");
  }
}

// ----------------------------------------

AverageComplementarity::AverageComplementarity(double lambdaStar)
{
  initialize(lambdaStar);
}

AverageComplementarity::~AverageComplementarity()
{
  // Nothing needs.
}

void AverageComplementarity::initialize(double lambdaStar)
{
  initial = lambdaStar*lambdaStar;
  current = initial;
  // rMessage("initial average = " << initial);
}

void AverageComplementarity::initialize(Solutions& initPt)
{
  int nDim = initPt.nDim;
  Lal::let(initial,'=',initPt.xMat,'.',initPt.zMat);
  initial /= nDim;
  current = initial;
}

void AverageComplementarity::update(Solutions& currentPt)
{
  int nDim = currentPt.nDim;
  Lal::let(current,'=',currentPt.xMat,'.',currentPt.zMat);
  current /= nDim;
}

void AverageComplementarity::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }

  fprintf(fpout,"mu0 = %8.3e\n",initial);
  fprintf(fpout,"mu  = %8.3e\n",current);
}

//--------------------------------------------------


RatioInitResCurrentRes::RatioInitResCurrentRes()
{
  primal = 0.0;
  dual   = 0.0;
}

RatioInitResCurrentRes::~RatioInitResCurrentRes()
{
  // Nothing needs.
}

RatioInitResCurrentRes::RatioInitResCurrentRes(Parameter& param,
					       Residuals& initRes)
{
  initialize(param,initRes);
}

void RatioInitResCurrentRes::initialize(Parameter& param,
					Residuals& initRes)
{
  double accuracy = param.epsilonDash;
  if (initRes.normPrimalVec < accuracy) {
    primal = 0.0;
  } else {
    primal = 1.0;
  }
  if (initRes.normDualMat < accuracy) {
    dual = 0.0;
  } else {
    dual = 1.0;
  }
}

void RatioInitResCurrentRes::update(Switch& reduction,
				    StepLength& alpha)
{
  if (reduction.switchType==Switch::CENTERING) {
    // At least one of primal or dual is infeasible
    primal = fabs((1.0-alpha.primal)*primal);
    dual   = fabs((1.0-alpha.dual  )*dual  );
  }
}

void RatioInitResCurrentRes::update_exact(Residuals& initRes,
					  Residuals& currentRes,
					  Parameter& param)
{
  if (initRes.normPrimalVec
      > param.epsilonDash * 1.0e-2) {
    primal = currentRes.normPrimalVec / initRes.normPrimalVec;
  }
  else {
    primal = 0.0;
  }
  if (initRes.normDualMat
      > param.epsilonDash * 1.0e-2) {
    dual = currentRes.normDualMat / initRes.normDualMat;
  }
  else {
    dual = 0.0;
  }
}

void RatioInitResCurrentRes::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }

  fprintf(fpout,"theta.primal = %8.3e\n",primal);
  fprintf(fpout,"theta.dual   = %8.3e\n",dual);
}

//---------------------------------------------------

SolveInfo::SolveInfo()
{
  rho = 0.0;
  etaPrimal = 0.0;
  etaDual   = 0.0;
  objValPrimal = 0.0;
  objValDual   = 0.0;
}

SolveInfo::SolveInfo(InputData& inputData, Solutions& currentPt, 
		     double mu0, double omegaStar)
{
  initialize(inputData,currentPt,mu0,omegaStar);
}

SolveInfo::~SolveInfo()
{
  // Nothing needs.
}

void SolveInfo::initialize(InputData& inputData, Solutions& currentPt, 
			   double mu0, double omegaStar)
{
  int nDim = currentPt.nDim;
  Vector& b = inputData.b;
  SparseLinearSpace& C = inputData.C;

  rho = 1.0;
  etaPrimal = omegaStar * nDim * mu0;
  etaDual   = omegaStar * nDim * mu0;
  Lal::let(objValPrimal,'=',C,'.',currentPt.xMat);
  Lal::let(objValDual  ,'=',b,'.',currentPt.yVec);
}

void SolveInfo::update(InputData& inputData,
		       DenseLinearSpace& initPt_xMat, 
		       DenseLinearSpace& initPt_zMat, 
		       Solutions& currentPt,
		       Residuals& currentRes,
		       AverageComplementarity& mu,
		       RatioInitResCurrentRes& theta,
		       Parameter& param)
{
  int nDim = currentPt.nDim;
  Vector& b = inputData.b;
  SparseLinearSpace& C = inputData.C;

  Lal::let(objValPrimal,'=',C,'.',currentPt.xMat);
  Lal::let(objValDual  ,'=',b,'.',currentPt.yVec);
  double primal = theta.primal;
  double dual   = theta.dual;
  double omega  = param.omegaStar;
  rho = 0.0;
  double x0z0     = nDim*mu.initial;
  double xMatzMat = nDim*mu.current;
  double x0zMat   = 0.0;
  double xMatz0   = 0.0;
  Lal::let(x0zMat,'=',initPt_xMat,'.',currentPt.zMat);
  Lal::let(xMatz0,'=',currentPt.xMat,'.',initPt_zMat);

  double accuracy = param.epsilonDash;

  if (currentRes.normPrimalVec <= accuracy) {
    // rMessage("primal accuracy");
    if (xMatz0 < etaPrimal) {
      etaPrimal = xMatz0;
    }
  }
  if (currentRes.normDualMat <= accuracy) {
    // rMessage("dual accuracy");
    if (x0zMat < etaDual) {
      etaDual = x0zMat;
    }
  }

  // primal is infeasible and dual is feasible
  if (currentRes.normPrimalVec > accuracy
      && currentRes.normDualMat <= accuracy) {
    rho = primal*x0zMat
      / ((primal+(1.0-primal)*omega)*etaDual + xMatzMat);
  }

  // primal is feasible and dual is infeasible
  if (currentRes.normPrimalVec <= accuracy
      && currentRes.normDualMat > accuracy) {
    rho = dual*xMatz0
      / ((dual+(1.0-dual)*omega)* etaPrimal + xMatzMat);
  }
  
  // primal and dual are infeasible
  if (currentRes.normPrimalVec > accuracy
      && currentRes.normDualMat > accuracy) {
    rho = (dual*xMatz0+primal*x0zMat)
      / ((primal*dual
	  + omega *(primal*(1.0-dual) + (1.0-primal)*dual))* x0z0
	 + xMatzMat);
  }
  // rMessage("eta Primal = " << etaPrimal);
  // rMessage("eta Dual = " << etaDual);
}


void SolveInfo::update(double& lambda, 
		       InputData& inputData,
		       Solutions& currentPt,
		       Residuals& currentRes,
		       AverageComplementarity& mu,
		       RatioInitResCurrentRes& theta,
		       Parameter& param)
{
  int nDim = currentPt.nDim;
  Vector& b = inputData.b;
  SparseLinearSpace& C = inputData.C;

  Lal::let(objValPrimal,'=',C,'.',currentPt.xMat);
  Lal::let(objValDual  ,'=',b,'.',currentPt.yVec);
  double primal = theta.primal;
  double dual   = theta.dual;
  double omega  = param.omegaStar;
  rho = 0.0;
  double x0z0     = nDim*mu.initial;
  double xMatzMat = nDim*mu.current;
  double x0zMat   = 0.0;
  double xMatz0   = 0.0;

  for (int l=0; l<currentPt.xMat.SDP_nBlock; l++){
    int dim = currentPt.xMat.SDP_block[l].nRow; 
    for (int k=0; k<dim; k++){
      x0zMat += lambda * currentPt.zMat.SDP_block[l].de_ele[k*dim+k];
      xMatz0 += lambda * currentPt.xMat.SDP_block[l].de_ele[k*dim+k];
    }
  }
  for (int l=0; l<currentPt.xMat.SOCP_nBlock; l++){
    rError("no support for SOCP");
  }
  for (int l=0; l<currentPt.xMat.LP_nBlock; l++){
    x0zMat += lambda * currentPt.zMat.LP_block[l];
    xMatz0 += lambda * currentPt.xMat.LP_block[l];
  }


  double accuracy = param.epsilonDash;

  if (currentRes.normPrimalVec <= accuracy) {
    // rMessage("primal accuracy");
    if (xMatz0 < etaPrimal) {
      etaPrimal = xMatz0;
    }
  }
  if (currentRes.normDualMat <= accuracy) {
    // rMessage("dual accuracy");
    if (x0zMat < etaDual) {
      etaDual = x0zMat;
    }
  }

  // primal is infeasible and dual is feasible
  if (currentRes.normPrimalVec > accuracy
      && currentRes.normDualMat <= accuracy) {
    rho = primal*x0zMat
      / ((primal+(1.0-primal)*omega)*etaDual + xMatzMat);
  }

  // primal is feasible and dual is infeasible
  if (currentRes.normPrimalVec <= accuracy
      && currentRes.normDualMat > accuracy) {
    rho = dual*xMatz0
      / ((dual+(1.0-dual)*omega)* etaPrimal + xMatzMat);
  }
  
  // primal and dual are infeasible
  if (currentRes.normPrimalVec > accuracy
      && currentRes.normDualMat > accuracy) {
    rho = (dual*xMatz0+primal*x0zMat)
      / ((primal*dual
	  + omega *(primal*(1.0-dual) + (1.0-primal)*dual))* x0z0
	 + xMatzMat);
  }
  // rMessage("eta Primal = " << etaPrimal);
  // rMessage("eta Dual = " << etaDual);
}

// 2007/09/13 kazuhide nakata  
// print information of ObjVal, residual, gap, complementarity
//   b^T y       + R \bullet X =  value,    norm(r),   norm(Z)
//   C \bullet X + r^T y       =  value,    norm(R),   norm(X)
//     gap                         gap,     mu * nDim    
void SolveInfo::check(InputData& inputData,
					  Solutions& currentPt,
					  Residuals& currentRes,
					  AverageComplementarity& mu,
					  RatioInitResCurrentRes& theta,
					  Parameter& param)
{
  double tmp,tmp1p,tmp1d,tmp2p,tmp2d,tmp3p,tmp3d,tmp4,tmp5p,tmp5d;

  Lal::let(tmp,'=',inputData.b,'.',currentPt.yVec);
  tmp1p = - tmp;
  printf("Primal: %9.1e",tmp1p);
  Lal::let(tmp,'=',currentRes.dualMat,'.',currentPt.xMat);
  tmp2p = -tmp;
  printf(" + %9.1e",tmp2p);
  tmp3p = tmp1p + tmp2p;
  printf(" = %9.1e",tmp3p);
  printf(",   residual:%-9.1e",currentRes.normDualMat);
  tmp5p = currentRes.computeMaxNorm(currentPt.zMat);
  printf(" norm:%-9.1e\n",tmp5p);

  Lal::let(tmp,'=',inputData.C,'.',currentPt.xMat);
  tmp1d = - tmp;
  printf("Dual:   %9.1e",tmp1d);
  Lal::let(tmp,'=',currentRes.primalVec,'.',currentPt.yVec);
  tmp2d = -tmp;
  printf(" + %9.1e",tmp2d);
  tmp3d = tmp1d + tmp2d;
  printf(" = %9.1e",tmp3d);
  printf(",   residual:%-9.1e", currentRes.normPrimalVec);
  tmp5d = currentRes.computeMaxNorm(currentPt.xMat);
  printf(" norm:%-9.1e\n",tmp5d);

  tmp4 = tmp1p - tmp1d;
  printf("P-D:    %9.1e",tmp4);
  tmp4 = tmp3p - tmp3d;
  printf("               %9.1e",tmp4);
  tmp4 = mu.current * currentPt.nDim;
  printf(",    mu * n:%-9.1e\n",tmp4);

}


void SolveInfo::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }

  fprintf(fpout,"rSolveInfo.rho          = %8.3e\n",rho);
  fprintf(fpout,"rSolveInfo.etaPrimal    = %8.3e\n",etaPrimal);
  fprintf(fpout,"rSolveInfo.etaDual      = %8.3e\n",etaDual);
  fprintf(fpout,"rSolveInfo.objValPrimal = %8.3e\n",objValPrimal);
  fprintf(fpout,"rSolveInfo.objValDual   = %8.3e\n",objValDual);
}

// ----------------------------------------------------

Phase::Phase()
{
  nDim = 0;
  value = SolveInfo::noINFO;
}

Phase::Phase(Residuals& initRes, SolveInfo& solveInfo,
	       Parameter& param, int nDim)
{
  initialize(initRes,solveInfo,param,nDim);
}

Phase::~Phase()
{
  // Nothing needs.
}

bool Phase::initialize(Residuals& initRes,
			SolveInfo& solveInfo,
			Parameter& param, int nDim)
{
  this->nDim = nDim;
  return updateCheck(initRes,solveInfo,param);
}

bool Phase::updateCheck(Residuals& currentRes,
			  SolveInfo& solveInfo,
			  Parameter& param)
{
  const double NONZERO = 1.0e-6;
  double accuracy = param.epsilonDash;
  value = SolveInfo::noINFO;

  if (currentRes.normPrimalVec <= accuracy) {
    if (currentRes.normDualMat <= accuracy) {
      value = SolveInfo::pdFEAS;
    } else {
      value = SolveInfo::pFEAS;
    }
  }
  if (value==SolveInfo::noINFO
      && currentRes.normDualMat <= accuracy) {
    value = SolveInfo::dFEAS;
  }
  if (value==SolveInfo::pdFEAS) {
    double mean = (fabs(solveInfo.objValPrimal)+
		   fabs(solveInfo.objValDual)) / 2.0;
    double PDgap = fabs(solveInfo.objValPrimal - solveInfo.objValDual);

    double dominator;
    if (mean < 1.0) {
      dominator = 1.0;
    } else {
      dominator = mean;
    }
    #if 0
    rMessage("PDgap = " << PDgap);
    rMessage("dominator = " << dominator);
    rMessage("PDgap/dominator = " << PDgap/dominator);
    #endif
    if (PDgap/dominator <= param.epsilonStar) {
      value = SolveInfo::pdOPT;
      return false;
    }
  }
  if (value == SolveInfo::noINFO
      && solveInfo.rho > 1.0+NONZERO) {
    rMessage("pdINF criteria");
    value = SolveInfo::pdINF;
    return false;
  }
  if (value == SolveInfo::pFEAS) {
    #if REVERSE_PRIMAL_DUAL
    if (solveInfo.objValPrimal<=-param.upperBound) {
      rMessage("pUNBD criteria");
      value = SolveInfo::pUNBD;
      return false;
    }
    #else
    if (solveInfo.objValPrimal<=param.lowerBound) {
      rMessage("pdINF criteria");
      value = SolveInfo::pUNBD;
      return false;
    }
    #endif
    if (solveInfo.rho > 1.0+NONZERO) {
      rMessage("pFEAS_dINF criteria");
      value = SolveInfo::pFEAS_dINF;
      return false;
    }
  }

  if (value == SolveInfo::dFEAS) {
    #if REVERSE_PRIMAL_DUAL
    if (solveInfo.objValDual>=-param.lowerBound) {
      rMessage("dUNBD criteria");
      value = SolveInfo::dUNBD;
      return false;
    }
    #else
    if (solveInfo.objValDual>=param.upperBound) {
      rMessage("dUNBD criteria");
      value = SolveInfo::dUNBD;
      return false;
    }
    #endif
    if (solveInfo.rho > 1.0+NONZERO) {
      rMessage("pINF_dFEAD criteria");
      value = SolveInfo::pINF_dFEAS;
      return false;
    }
  }
  #if 0
  rMessage("phase =");
  display();
  #endif
  return true;
}

void Phase::reverse()
{
  #if REVERSE_PRIMAL_DUAL
  switch (value) {
  case SolveInfo::noINFO    :                              ; break;
  case SolveInfo::pFEAS     : value = SolveInfo::dFEAS     ; break;
  case SolveInfo::dFEAS     : value = SolveInfo::pFEAS     ; break;
  case SolveInfo::pdFEAS    :                              ; break;
  case SolveInfo::pdINF     :                              ; break;
  case SolveInfo::pFEAS_dINF: value = SolveInfo::pINF_dFEAS; break;
  case SolveInfo::pINF_dFEAS: value = SolveInfo::pFEAS_dINF; break;
  case SolveInfo::pdOPT     :                              ; break;
  case SolveInfo::pUNBD     : value = SolveInfo::dUNBD     ; break;
  case SolveInfo::dUNBD     : value = SolveInfo::pUNBD     ; break;
  default: break;
  }
  #else
  // do nothing
  #endif
}


void Phase::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  char* str;
  switch (value) {
  case SolveInfo::noINFO    : str = (char *)"noINFO    "; break;
  case SolveInfo::pFEAS     : str = (char *)"pFEAS     "; break;
  case SolveInfo::dFEAS     : str = (char *)"dFEAS     "; break;
  case SolveInfo::pdFEAS    : str = (char *)"pdFEAS    "; break;
  case SolveInfo::pdINF     : str = (char *)"pdINF     "; break;
  case SolveInfo::pFEAS_dINF: str = (char *)"pFEAS_dINF"; break;
  case SolveInfo::pINF_dFEAS: str = (char *)"pINF_dFEAS"; break;
  case SolveInfo::pdOPT     : str = (char *)"pdOPT     "; break;
  case SolveInfo::pUNBD     : str = (char *)"pUNBD     "; break;
  case SolveInfo::dUNBD     : str = (char *)"dUNBD     "; break;
  default:
    str = (char *)"phase error";
    rMessage("rPhase:: phase error");
    break;
  }
  fprintf(fpout,"phase.value  = %s\n",str);
}

} // end of namespace 'sdpa'

