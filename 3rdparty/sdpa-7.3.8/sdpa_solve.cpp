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
/*--------------------------------------------------
  sdpa_solve.cpp
--------------------------------------------------*/

#include "sdpa_call.h"
#include "sdpa_linear.h"
#include "sdpa_io.h"
using namespace sdpa;

void SDPA::initializeSolve()
{
  TimeStart(FILE_CHANGE_START1);
  // if possible , change C and A to Dense
  inputData.C.changeToDense();
  for (int k=0; k<m; ++k) {
    inputData.A[k].changeToDense();
  }
  TimeEnd(FILE_CHANGE_END1);
  com.FileChange += TimeCal(FILE_CHANGE_START1,
			    FILE_CHANGE_END1);
  com.TotalTime += TimeCal(FILE_CHANGE_START1,
			    FILE_CHANGE_END1);
  inputData.initialize_index();
  newton.initialize(m,bs);
  int nBlock2 = bs.SDP_nBlock + bs.SOCP_nBlock + bs.LP_nBlock;
  chordal.initialize(&newton.sparse_bMat);
  chordal.ordering_bMat(m, nBlock2, inputData, Display, fpout);
  newton.initialize_bMat(m, chordal, inputData, Display, fpout);
  newton.computeFormula_SDP(inputData,0.0,KAPPA);

  work.initialize(m,bs);
  if (isInitPoint == false) {
    mu.initialize(param.lambdaStar);
    initRes.initialize(m, bs, inputData, currentPt);
    currentRes.copyFrom(initRes);
    beta.initialize(param.betaStar);
    theta.initialize(param, initRes);
    solveInfo.initialize(inputData, currentPt, mu.initial,
			 param.omegaStar);
    phase.initialize(initRes, solveInfo, param, currentPt.nDim);
  }
  // writeInputSparse((char*)"tmp.dat-s",(char*)"%+8.3e");
}

void SDPA::solve()
{
  if (isInitPoint == true) {
    mu.initialize(currentPt);
    currentPt.computeInverse(work,com);
    initPt_xMat.copyFrom(currentPt.xMat);
    initPt_zMat.copyFrom(currentPt.zMat);
    initRes.initialize(m, bs, inputData, currentPt);
    currentRes.copyFrom(initRes);
    theta.initialize(param, initRes);
    solveInfo.initialize(inputData, currentPt, mu.initial,
			 param.omegaStar);
    phase.initialize(initRes, solveInfo, param, currentPt.nDim);
  }
  
  pIteration = 0;
  TimeStart(MAIN_LOOP_START1);
  IO::printHeader(fpout,Display);
  while (phase.updateCheck(currentRes, solveInfo, param)
	 && pIteration < param.maxIteration) {
    // rMessage(" turn hajimari " << pIteration );
    // Mehrotra's Predictor
    TimeStart(MEHROTRA_PREDICTOR_START1);
    // set variable of Mehrotra
    reduction.MehrotraPredictor(phase);
    beta.MehrotraPredictor(phase, reduction, param);

    // rMessage("reduction = "); reduction.display();
    // rMessage("phase = "); phase.display();
    // rMessage("beta.predictor.value = " << beta.value);
    // rMessage(" mu = " << mu.current);
    // rMessage("currentPt = "); currentPt.display();
    // rMessage("initRes = "); initRes.display();
    // rMessage("currentRes = "); currentRes.display();
    // inputData.display();

    bool isSuccessCholesky;
    isSuccessCholesky = newton.Mehrotra(Newton::PREDICTOR,
					m, inputData, chordal,
					currentPt, currentRes,
					mu, beta, reduction,
					phase, work, com,
					Display, fpout);
    if (isSuccessCholesky == false) {
      break;
    }
    // rMessage("newton predictor = "); newton.display();

    TimeEnd(MEHROTRA_PREDICTOR_END1);
    com.Predictor += TimeCal(MEHROTRA_PREDICTOR_START1,
			     MEHROTRA_PREDICTOR_END1);

    TimeStart(STEP_PRE_START1);
    alpha.MehrotraPredictor(inputData, currentPt, phase, newton,
			    work, com);
    // rMessage("alpha predictor = "); alpha.display();

    TimeStart(STEP_PRE_END1);
    com.StepPredictor += TimeCal(STEP_PRE_START1,STEP_PRE_END1);

    // rMessage("alphaStar = " << param.alphaStar);
    // Mehrotra's Corrector
    // rMessage(" Corrector ");

    TimeStart(CORRECTOR_START1);
    beta.MehrotraCorrector(phase,alpha,currentPt,
			   newton,mu,param);
    // beta.display();

    // rMessage("beta corrector = " << beta.value);

#if 1 // 2007/08/29 kazuhide nakata
    // add stopping criteria: objValPrimal < ObjValDual
    //	if ((pIteration > 10) &&
    if (phase.value == SolveInfo::pdFEAS
	&& ( beta.value> 5.0
	     || solveInfo.objValPrimal < solveInfo.objValDual + (1.0e-6))
	&& (fabs(solveInfo.objValPrimal) > 1.0e-4
	    && fabs(solveInfo.objValDual) > 1.0e-4) ) {
      rMessage("Strange behavior : primal < dual");
      break;
    }
#endif

    newton.Mehrotra(Newton::CORRECTOR, m, 
		    inputData, chordal, currentPt, currentRes,
		    mu, beta, reduction, phase,work,com,
		    Display, fpout);

    // rMessage("currentPt = "); currentPt.display();
    // rMessage("newton corrector = "); newton.display();

    TimeEnd(CORRECTOR_END1);
    com.Corrector += TimeCal(CORRECTOR_START1,
			     CORRECTOR_END1);
    TimeStart(CORRECTOR_STEP_START1);
    alpha.MehrotraCorrector(inputData, currentPt, phase,
			    reduction, newton, mu, theta,
			    work, param, com);
    // rMessage("alpha corrector = "); alpha.display();
    TimeEnd(CORRECTOR_STEP_END1);
    com.StepCorrector += TimeCal(CORRECTOR_STEP_START1,
				 CORRECTOR_STEP_END1);
    // the end of Corrector
    
    IO::printOneIteration(pIteration, mu, theta, solveInfo,
			  alpha, beta, fpout, Display);

    if (currentPt.update(alpha,newton,work,com)==false) {
      // if step length is too short,
      // we finish algorithm
      rMessage("cannot move: step length is too short");
      //   memo by kazuhide nakata
      //   StepLength::MehrotraCorrector
      //   thetaMax*mu.initial -> thetamax*thetaMax*mu.initial
      break;
    }

    // rMessage("currentPt = "); currentPt.display();
    // rMessage("updated");

    theta.update(reduction,alpha);
    mu.update(currentPt);
    currentRes.update(m,inputData, currentPt, com);
    theta.update_exact(initRes,currentRes, param);

    if (isInitPoint) {
      solveInfo.update(inputData, initPt_xMat, initPt_zMat, currentPt,
		       currentRes, mu, theta, param);
    } else {
      solveInfo.update(param.lambdaStar,inputData, currentPt,
		       currentRes, mu, theta, param);
    }
    // 2007/09/18 kazuhide nakata
    // print information of ObjVal, residual, gap, complementarity
    // solveInfo.check(inputData, currentPt,
    //                 currentRes, mu, theta, param);
    pIteration++;
  } // end of MAIN_LOOP

  if (pIteration == param.maxIteration) {
    rMessage("maxIteration is reached");
  }
  TimeEnd(MAIN_LOOP_END1);

  com.MainLoop = TimeCal(MAIN_LOOP_START1,
			 MAIN_LOOP_END1);
  com.TotalTime += com.MainLoop;
  currentRes.compute(m,inputData,currentPt);
#if REVERSE_PRIMAL_DUAL
  Lal::let(currentPt.yVec,'=',currentPt.yVec,'*',&DMONE);
  phase.reverse();
#endif
  IO::printLastInfo(pIteration, mu, theta, solveInfo, alpha, beta,
		    currentRes, phase, currentPt, 
		    inputData, work, com.TotalTime, com,
		    param, fpout, Display);
  IO::printSolution(bs, currentPt, param, fpout);
  // com.display(fpout);

  if (Display) {
    fprintf(Display,   "  main loop time = %.6f\n",com.MainLoop);
    fprintf(Display,   "      total time = %.6f\n",com.TotalTime);
    fprintf(Display,   "file  check time = %.6f\n",com.FileCheck);
    fprintf(Display,   "file change time = %.6f\n",com.FileChange);
    fprintf(Display,   "file   read time = %.6f\n",com.FileRead);
  }
  if (fpout) {
    fprintf(fpout,   "    main loop time = %.6f\n",com.MainLoop);
    fprintf(fpout,   "        total time = %.6f\n",com.TotalTime);
    fprintf(fpout,   "  file  check time = %.6f\n",com.FileCheck);
    fprintf(fpout,   "  file change time = %.6f\n",com.FileChange);
    fprintf(fpout,   "  file   read time = %.6f\n",com.FileRead);
  }
}
