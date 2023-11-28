#ifndef _OPENGR_ACCELERATORS_SDPA_WRAPPER_H
#define _OPENGR_ACCELERATORS_SDPA_WRAPPER_H

#include <Eigen/Dense>
#include <sdpa_call.h>
#include <cstdio>

#include <gr/utils/logger.h>

namespace gr{

    template <typename Scalar_>
    class SDPA_WRAPPER{
    
    public:
        using MatrixX = typename Eigen::Matrix<Scalar_, Eigen::Dynamic, Eigen::Dynamic>;

        void Solve(Eigen::Ref<const MatrixX> C, Eigen::Ref<MatrixX> G, const Utils::Logger& logger, const int d, const int m);

    private:
        void LogResults(SDPA& Problem);
    };

    /// Computes the SDP (P2) as described in [this paper](https://arxiv.org/abs/1306.5226).
    /// P2: min(Tr(CG)) subject to G >= 0, G_ii = I_d (1<=i<=m).
    /// @param [in] C "patch-stress" matrix.
    /// @param [in] d Dimension.
    /// @param [in] m Number of patches.
    /// @param [out] G Solution of the SDP (P2)
    template <typename Scalar_>
    void SDPA_WRAPPER<Scalar_>::Solve(Eigen::Ref<const MatrixX> C, Eigen::Ref<MatrixX> G, const Utils::Logger& logger, const int d, const int m){
        SDPA	Problem;

        // All parameteres are renewed
        Problem.setParameterType(SDPA::PARAMETER_STABLE_BUT_SLOW);
        Problem.printParameters(stdout);

        int mDIM   = d*(d+1)/2*m;
        int nBlock = 1;
        Problem.inputConstraintNumber(mDIM);
        Problem.inputBlockNumber(nBlock);
        Problem.inputBlockSize(1,d*m);
        Problem.inputBlockType(1,SDPA::SDP);

        Problem.initializeUpperTriangleSpace();

        // c vec
        int cnt = 1;
        for(int i = 0; i < m; i++){
            for(int j = 0; j < d; j++){
            Problem.inputCVec(cnt++,1);
            for(int k = j+1; k < d; k++)
                Problem.inputCVec(cnt++,0);
            }
        } 

        // F0 = -C
        for(int i = 0; i < m*d; i++)
            for(int j = i; j < m*d; j++)
            Problem.inputElement(0, 1, i+1, j+1, -C(i,j));

        // Fi
        cnt = 1;
        for(int k = 0; k < m; k++){
            for(int i = 0; i < d; i++){
            for(int j = i; j < d; j++){
                Problem.inputElement(cnt++, 1, k*d+i+1, k*d+j+1, 1);
            }
            }
        }

        Problem.initializeUpperTriangle();
        Problem.initializeSolve();
        Problem.solve();

        if(logger.logLevel() == Utils::LogLevel::Verbose)
            LogResults(Problem);


        double* yMat = Problem.getResultYMat(1);
        Eigen::Map<Eigen::MatrixXd> Gmap(yMat, d*m, d*m);
        G = Gmap;
    }


    template <typename Scalar_>
    void SDPA_WRAPPER<Scalar_>::LogResults(SDPA& Problem){

        fprintf(stdout, "\nStop iteration = %d\n",
            Problem.getIteration());
        char phase_string[30];
        Problem.getPhaseString(phase_string);
        fprintf(stdout, "Phase          = %s\n", phase_string);
        fprintf(stdout, "objValPrimal   = %+10.6e\n",
            Problem.getPrimalObj());
        fprintf(stdout, "objValDual     = %+10.6e\n",
            Problem.getDualObj());
        fprintf(stdout, "p. feas. error = %+10.6e\n",
            Problem.getPrimalError());
        fprintf(stdout, "d. feas. error = %+10.6e\n\n",
            Problem.getDualError());
    }

}
#endif // SDPA_WRAPPER_H