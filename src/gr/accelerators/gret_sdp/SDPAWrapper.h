#ifndef _OPENGR_ACCELERATORS_SDPA_WRAPPER_H
#define _OPENGR_ACCELERATORS_SDPA_WRAPPER_H

#include <Eigen/Dense>
#include <sdpa_call.h>

namespace gr{

    template <typename Scalar_>
    class SDPA_WRAPPER{
    
    public:
        using MatrixX = typename Eigen::Matrix<Scalar_, Eigen::Dynamic, Eigen::Dynamic>;

        void Solve(Eigen::Ref<const MatrixX> C, Eigen::Ref<MatrixX> G, const int d, const int m);
    };
    
    template <typename Scalar_>
    void SDPA_WRAPPER<Scalar_>::Solve(Eigen::Ref<const MatrixX> C, Eigen::Ref<MatrixX> G, const int d, const int m){
        SDPA	Problem;

        // All parameteres are renewed
        Problem.setParameterType(SDPA::PARAMETER_STABLE_BUT_SLOW);

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


        double* yMat = Problem.getResultYMat(1);
        Eigen::Map<Eigen::MatrixXd> Gmap(yMat, d*m, d*m);
        G = Gmap;
    }

}
#endif // SDPA_WRAPPER_H