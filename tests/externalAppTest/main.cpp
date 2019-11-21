#include "gr/algorithms/match4pcsBase.h"
#include "gr/algorithms/FunctorSuper4pcs.h"
#include "gr/utils/geometry.h"
#include <gr/algorithms/PointPairFilter.h>

#include <Eigen/Dense>


int main(int argc, char **argv) {
  using namespace gr;
  using namespace std;

  using TrVisitor = gr::DummyTransformVisitor;

  using MatcherType = gr::Match4pcsBase<gr::FunctorSuper4PCS, gr::Point3D<float>, TrVisitor, gr::AdaptivePointFilter, gr::AdaptivePointFilter::Options>;
  using OptionType  = typename MatcherType::OptionsType;
  using SamplerType = gr::UniformDistSampler<gr::Point3D<float> >;

  vector<Point3D<float> > set1, set2;
  vector<Eigen::Matrix2f> tex_coords1, tex_coords2;
  vector<typename Point3D<float>::VectorType> normals1, normals2;
  vector<std::string> mtls1, mtls2;

  // dummy calls, to test symbols accessibility
  // check availability of the Utils functions
  Utils::CleanInvalidNormals(set1, normals1);

  // Our matcher.
  OptionType options;

  // Set parameters.
  typename MatcherType::MatrixType mat;
  double overlap (1);
  options.configureOverlap(overlap);

  typename Point3D<float>::Scalar score = 0;

  constexpr Utils::LogLevel loglvl = Utils::Verbose;
  Utils::Logger logger(loglvl);
  SamplerType sampler;
  TrVisitor visitor;

  MatcherType matcher(options, logger);
  score = matcher.ComputeTransformation(set1, set2, mat, sampler, visitor);

  logger.Log<Utils::Verbose>( "Score: ", score );

  return 0;
}

