
#include <algorithm> // std::transform
#include <utility>   // std::pair
#include <Eigen/Dense>


#include "gr/algorithms/GRET_SDP.h"
#include "gr/utils/logger.h"
namespace gr {


template <typename PointType, typename TransformVisitor, template < class, class > typename ... OptExts>
template<typename Solver, typename PatchRange>
void GRET_SDP<PointType, TransformVisitor, OptExts ... >::RegisterPatches(const PatchRange& patches, const int num_global_coordinates, TransformVisitor& v){
    n = num_global_coordinates;
    d = PointType::dim();
    m = patches.size();

    MatrixX L = MatrixX::Zero(n+m, n+m);
    MatrixX D = MatrixX::Zero(m*d, m*d);
    MatrixX C = MatrixX::Zero(m*d, m*d);
    MatrixX B = MatrixX::Zero(m*d, n+m);
    Linv = MatrixX::Zero(n+m, n+m);

    // Compute matrices B, D and C
    for(int i = 0; i < patches.size(); i++){
        for(const std::pair<PointType,int> &point_with_index : patches[i]){
            int k = point_with_index.second;
            VectorType pos = point_with_index.first.pos();
            
            L(k,k)++;
            L(k,n+i)--;
            L(n+i,k)--;
            L(n+i,n+i)++;

            B.block(i*d,k,d,1) += pos;
            B.block(i*d,n+i,d,1) -= pos;

            D.block(i*d, i*d, d, d) += pos*pos.transpose();
        }
    }

    // compute Linv, the Moore-Penrose pseudoinverse of L
    Linv = L.completeOrthogonalDecomposition().pseudoInverse();

    // compute C
    C = D - B * Linv * B.transpose();

    // solve the SDP (P2) using C
    MatrixX G(m*d, m*d);

    
    Solver solver;
    solver.Solve(C,G,MatchBaseType::logger_,d,m);

    // TODO: compute W
    // compute top d eigenvalues and eigenvectors
    Eigen::EigenSolver<MatrixX> s(G);
    Eigen::VectorXcd eigvals = s.eigenvalues();
    Eigen::MatrixXcd eigvecs = s.eigenvectors();
    std::vector<double> re_eigvals;
    std::vector<Eigen::VectorXd> re_eigvecs;

    for(int i = 0; i < eigvals.rows(); i++) {
        if(eigvals(i).imag() == 0){
            re_eigvals.push_back(eigvals(i).real());
            re_eigvecs.push_back(eigvecs.col(i).real());
        }
    }   

  std::vector<std::pair<double, Eigen::VectorXd>> eig_pairs;
  eig_pairs.reserve(d);
  std::transform(re_eigvals.begin(), re_eigvals.end(), re_eigvecs.begin(), std::back_inserter(eig_pairs),
               [](double a, const Eigen::VectorXd& b) { return std::make_pair(a, b); });

  sort(eig_pairs.begin(), eig_pairs.end(),
    [&](std::pair<double, Eigen::VectorXd>& a, std::pair<double, Eigen::VectorXd>& b) {
        return (a.first > b.first);
    }
  );

  // construct W
  MatrixX W(d, m*d);
  for(int i = 0; i < d; i++){
      W.row(i) = std::sqrt(eig_pairs[i].first) * eig_pairs[i].second.transpose();
  }

  // compute transformations O
  O = MatrixX(d, m*d);
  for(int i = 0; i < m; i++){
      Eigen::Ref<MatrixX> w(W.block(0, i*d, d, d));
      Eigen::JacobiSVD<MatrixX> svd(w, Eigen::ComputeFullU | Eigen::ComputeFullV);
      O.block(0, i*d, d, d) = svd.matrixU() * svd.matrixV().transpose();
  }

  OB = O*B;
}

// returns registered points
template <typename PointType, typename TransformVisitor, template < class, class > typename ... OptExts>
template<typename PointRange>
void GRET_SDP<PointType, TransformVisitor, OptExts ... >::getRegisteredPatches(PointRange& registered_points){
  MatrixX R(d, n);
  R = OB*Linv.block(0,0,n+m,n);
  registered_points.reserve(n);
  for(int k = 0; k < n; k++){
    registered_points.push_back(PointType(VectorType(R.block(0,k,d,1))));
  }
}

// returns transformations
template <typename PointType, typename TransformVisitor, template < class, class > typename ... OptExts>
template<typename TrRange>
void GRET_SDP<PointType, TransformVisitor, OptExts ... >::getTransformations(TrRange& transformations){
  MatrixX T(d, m);
  T = OB*Linv.block(0, n, n+m, m);

  // store transformations
  MatrixType trafo;
  transformations.reserve(m);
  for(int i = 0; i < m; i++){
    trafo = MatrixType::Zero();
    trafo.block(0, 0, d, d) = O.block(0, i*d, d, d);
    trafo.block(0, d, d, 1) = T.block(0, i, d, 1);
    trafo(d,d) = 1;
    transformations.emplace_back(trafo);
  }
}

} // namespace gr
