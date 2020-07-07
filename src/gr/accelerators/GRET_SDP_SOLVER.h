#ifndef _OPENGR_ACCELERATORS_GRET_SDP_SOLVER_H
#define _OPENGR_ACCELERATORS_GRET_SDP_SOLVER_H

namespace gr{

    template <typename Scalar_>
    class GRET_SDP_SOLVER  {

    public:
        using Scalar = Scalar_;
        using MatrixX = typename Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

        // computes registered points and corresponding transformations
        virtual void Solve(Eigen::Ref<const MatrixX> C, Eigen::Ref<MatrixX> G, const int d, const int m) = 0;

    };
}


#endif // GRET_SDP_SOLVER_H