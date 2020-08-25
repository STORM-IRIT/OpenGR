#ifndef _OPENGR_ACCELERATORS_MOSEK_WRAPPER_H
#define _OPENGR_ACCELERATORS_MOSEK_WRAPPER_H

#include <Eigen/Dense>
#include <fusion.h>

#include <gr/utils/logger.h>

namespace gr{

    template <typename Scalar_>
    class MOSEK_WRAPPER{
    
    public:
        using MatrixX = typename Eigen::Matrix<Scalar_, Eigen::Dynamic, Eigen::Dynamic>;

        void Solve(Eigen::Ref<const MatrixX> C, Eigen::Ref<MatrixX> G, const Utils::Logger& logger, const int d, const int m);
    };
    
    /// Computes the SDP (P2) as described in [this paper](https://arxiv.org/abs/1306.5226).
    /// P2: min(Tr(CG)) subject to G >= 0, G_ii = I_d (1<=i<=m).
    /// @param [in] C "patch-stress" matrix.
    /// @param [in] d Dimension.
    /// @param [in] m Number of patches.
    /// @param [out] G Solution of the SDP (P2)
    template <typename Scalar>
    void MOSEK_WRAPPER<Scalar>::Solve(Eigen::Ref<const MatrixX> C, Eigen::Ref<MatrixX> G, const Utils::Logger& logger, int d, const int m){
        using namespace mosek::fusion;
        using namespace monty;

        Model::t M = new Model("gret-sdp"); auto _M = finally([&]() { M->dispose(); });

        auto c_ptr = new_array_ptr<Scalar, 2>(shape_t<2>(m*d,m*d));
        std::copy(C.data(), C.data()+C.size(), c_ptr->begin());
        Matrix::t C_ = Matrix::dense(c_ptr);

        if(logger.logLevel() == Utils::LogLevel::Verbose)
            M->setLogHandler([ = ](const std::string & msg) { std::cout << msg << std::flush; } );

        auto G_ = M->variable(Domain::inPSDCone(m*d));

        for (int i = 0; i < m; i++){
            for (int j = 0; j < d; j++){
                M->constraint(G_->index(i*d+j, i*d+j),Domain::equalsTo(1.0));
            for (int k = j+1; k < d; k++)
                M->constraint(G_->index(i*d+j, i*d+k), Domain::equalsTo(0.0));
            }
        }  
        
        // Set the objective function to (Tr(CG)=sum(C⋅G))
        M->objective("obj", ObjectiveSense::Minimize, Expr::dot(C_, G_));

        // Solve the problem
        M->solve();
        
        // Get the solution values
        auto sol = G_->level();
        Eigen::Map<Eigen::MatrixXd> Gmap(sol->begin(), d*m, d*m);
        G = Gmap;
    }

}
#endif // MOSEK_WRAPPER_H