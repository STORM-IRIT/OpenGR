
#include <algorithm> // std::transform
#include <utility>   // std::pair
#include <Eigen/Dense>


#include "gr/algorithms/GRET_SDP.h"
#include "gr/utils/logger.h"


namespace gr {


template <typename PointType, typename TransformVisitor, template < class, class > class ... OptExts>
template <typename InputRange, template<typename> class Sampler>
void
GRET_SDP<PointType, TransformVisitor, OptExts ... >
::ComputeTransformations(const std::vector<InputRange>& P,
                         int N_,
                         int M_,
                         int d_,
                         Eigen::MatrixXd& L,
                         std::vector<MatrixType>& transformations_,
                         std::vector<VectorType>& translations_,
                         const Sampler<PointType>& sampler,
                         TransformVisitor& v){
    d = d_;
    N = N_;
    M = M_;

    // TODO: compute matrices B, D
    Eigen::Matrix<Scalar, M*d, N+M> B;
    Eigen::Matrix<Scalar, M*d, M*d> D;
    Eigen::Matrix<Scalar, N+M, N+M> Linv;
    Eigen::Matrix<Scalar, M*d, M*d> C;


    // construct B and D
    VectorType x = 0;
    for(int i = 0; i < M; i++){
        int num = 0;
        auto x_it = P[i].begin();
        for(int k = 0; k < N; k++){
            if(L(k,N+i) != 0){
                x = *x_it.pos(); x_it++;
                B.block<d, 1>(M*i, k) = x;
                B.block<d, 1>(M*i, N+i) += x;
                D.block<d, d>(M*i, M*i) += x * x.transpose();

            }
        }
    }

    // compute Linv, the Moore-Penrose pseudoinverse of L
    Linv = L.completeOrthogonalDecomposition().pseudoInverse();

    // compute C
    C = D - B * Linv * B.transpose();

    // TODO: solve the SDP (P2) using C
    Eigen::Matrix<Scalar, M*d, M*d> G;


    // TODO: compute W
    // compute top d eigenvalues and eigenvectors
    Eigen::EigenSolver<Eigen::Matrix<Scalar, M*d, M*d>> s(G);
    Eigen::VectorXcd eigvals = s.eigenvalues();
    Eigen::MatrixXcd eigvecs = s.eigenvectors();
    std::vector<double> re_eigvals;
    std::vector<Eigen::VectorXd> re_eigvecs;

    for(int i = 0; i < eigvals.rows(); i++) {
        if(vals(i).imag() == 0){
            re_vals.push_back(vals(i).real());
            re_vecs.push_back(vec.col(i).real());
        }
    }

    std::vector<std::pair<double, Eigen::VectorXd>> eig_pair;
    eig_pair.reserve(eigenvals.cols());
    std::transform(re_vals.begin(), re_vals.end(), re_vecs.begin(), std::back_inserter(eig_pair),
                [](double a, VectorXd b) { return std::make_pair(a, b); });

    sort(eig_pair.begin(), eig_pair.end(),
        [&](std::pair<double, VectorXd>& a, std::pair<double, VectorXd>& b) {
            return (a.first > b.first);
        }
    );

    // construct W
    Eigen::Matrix<Scalar, d, M*d> W;
    for(int di = 0; di < d; d++){
        W.row(di) = Eigen::sqrt(eig_pair[di].first * eig_pair[di].second.transpose());
    }


    // TODO: compute transformations O
    Eigen::Matrix<Scalar, d, M*d> O;
    for(int i = 0; i < M; i++){
        Eigen::Matrix<Scalar, d, d>& w = W.block<d, d>(0, i*d);
        Eigen::JacobiSVD<Eigen::Matrix<Scalar, d, d>> svd(w, Eigen::ComputeFullU | Eigen::ComputeFullV);
        O.block<d, d>(0, i*d) = svd.matrixU * svd.matrixV.transpose();
    }

    // TODO: compute Z
    Eigen::Matrix<Scalar, d, N+M> Z;
    Z = O*B*Linv;

    for(int i = 0; i < M; i++){
        transformations_.push_back(O.block<d, d>(0, i*d));
        translations_.push_back(Z.col(N+i));
    }
}


} // namespace gr
