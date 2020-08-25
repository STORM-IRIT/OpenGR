// Copyright 2014 Nicolas Mellado
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// -------------------------------------------------------------------------- //
//
// Authors: Nicolas Mellado
//
// This test check the validity of the pair extraction subroutine on random data
// of different dimensions (2,3 and 4).
//
//
// This test is part of the implementation of the Super 4-points Congruent Sets
// (Super 4PCS) algorithm presented in:

#include "gr/accelerators/gret_sdp/wrappers.h"
#include "gr/algorithms/GRET_SDP.h"
#include "gr/utils/timer.h"
#include "gr/accelerators/kdtree.h"
#include "gr/io/io.h"


#include <Eigen/Dense>

#include <fstream>
#include <iostream>
#include <string>

#include <stdlib.h>
#include <utility> // pair

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/filesystem.hpp>

#include "testing.h"

using namespace std;
using namespace gr;
namespace pt = boost::property_tree;
namespace fs = boost::filesystem;

using PointType = Point3D<double>;
using Scalar = PointType::Scalar;
enum {Dim = PointType::dim()};
using MatrixType = Eigen::Matrix<Scalar, Dim+1, Dim+1>;
using VectorType = Eigen::Matrix<Scalar, Dim, 1>;
using MatrixX = typename Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
using Pwi = pair<PointType, int>;
using PwiRange = vector<Pwi>;
using PointRange = vector<PointType>;
typedef Eigen::Transform<Scalar, Dim, Eigen::Affine> Transform;

constexpr Utils::LogLevel loglvl = Utils::Verbose;
Utils::Logger logger(loglvl);

constexpr int max_num_point_clouds = 5;

struct SamplingOptions {
    double delta = 0.01;
    SamplingOptions() {};
    SamplingOptions(double delta_) : delta(delta_) {} ;
};

// extracts point clouds and transformations from a stanford config file
void extractPCAndTrFromStandfordConfFile(
        const std::string &confFilePath,
        std::vector<Transform>& transforms,
        std::vector<PointRange>& point_clouds
        ){
    using namespace boost;
    using namespace std;


    std::vector<string> files;
    
    if(!(filesystem::exists(confFilePath) && filesystem::is_regular_file(confFilePath)))
      throw std::runtime_error("Config file does not exist or is no regular file.");

    // extract the working directory for the configuration path
    const std::string workingDir = filesystem::path(confFilePath).parent_path().native();    
    if(!filesystem::exists(confFilePath))
      throw std::runtime_error("Directory \"" + workingDir + "\" does not exist.");

    // read the configuration file and call the matching process
    std::string line;
    std::ifstream confFile;
    confFile.open(confFilePath);
    if(!confFile.is_open())
      throw std::runtime_error("Could not open config file.");

    int read_point_clouds = 0;
    while ( getline (confFile,line) && read_point_clouds < max_num_point_clouds)
    {
        std::istringstream iss (line);
        std::vector<string> tokens{istream_iterator<string>{iss},
                              istream_iterator<string>{}};

        // here we know that the tokens are:
        // [0]: keyword, must be bmesh
        // [1]: 3D object filename
        // [2-4]: target translation with previous object
        // [5-8]: target quaternion with previous object
        
        if (tokens.size() == 9){
            if (tokens[0].compare("bmesh") == 0){
                std::string inputfile = filesystem::path(confFilePath).parent_path().string()+string("/")+tokens[1];
                if(!(filesystem::exists(inputfile) && filesystem::is_regular_file(inputfile)))
                  throw std::runtime_error("File \"" + inputfile + "\" does not exist or is no regular file.");

                // build the Eigen rotation matrix from the rotation and translation stored in the files
                Eigen::Matrix<double, 3, 1> tr (
                            std::atof(tokens[2].c_str()),
                            std::atof(tokens[3].c_str()),
                            std::atof(tokens[4].c_str()));

                Eigen::Quaternion<double> quat(
                            std::atof(tokens[8].c_str()), // eigen starts by w
                            std::atof(tokens[5].c_str()),
                            std::atof(tokens[6].c_str()),
                            std::atof(tokens[7].c_str()));

                quat.normalize();

                Transform transform (Transform::Identity());
                transform.rotate(quat);
                transform.translate(-tr);

                transforms.push_back(transform);
                files.push_back(inputfile);
                
                read_point_clouds++;
            }
        }
    }
    confFile.close();

    std::ifstream pc_file;
    int num_point_clouds = files.size();
    point_clouds.resize(num_point_clouds);

    vector<Eigen::Matrix2f> tex_coords;
    vector<VectorType> normals;
    vector<tripple> tris;
    vector<string> mtls;

    IOManager iomanager;
    for(int i = 0; i < num_point_clouds; i++){
        const string& file = files[i];
         if(!iomanager.ReadObject((char *)file.c_str(), point_clouds[i], tex_coords, normals, tris, mtls) )
            throw std::runtime_error("Couldn't read ply file \"" + file + "\"");
    }
}

template <typename PointRange>
gr::KdTree<Scalar> constructKdTree(const PointRange& Q){
  size_t number_of_points = Q.size();
  // Build the kdtree.
  gr::KdTree<Scalar> kd_tree(number_of_points);

  for (size_t i = 0; i < number_of_points; ++i) {
      kd_tree.add(Q[i].pos());
  }
  kd_tree.finalize();
  return kd_tree;
}

// computes patches from point_clouds and returns number of global coordinates
int computePatches(const std::vector<PointRange>& point_clouds, const std::vector<Transform>& transforms, vector<PwiRange>& patches){
    using PointAdapter = gr::Point3D<Scalar>;
    UniformDistSampler<PointAdapter> sampler;

    int num_point_clouds = point_clouds.size();

    std::vector<PointRange> transformed_point_clouds(num_point_clouds);
    PointRange merged_point_cloud;
    for(size_t i = 0; i < num_point_clouds; i++){
        for(size_t j = 0; j < point_clouds[i].size(); j++){
            PointType transformed_point(transforms[i].inverse() * point_clouds[i][j].pos());
            transformed_point_clouds[i].push_back(transformed_point);
            merged_point_cloud.push_back(transformed_point);
        }
    }

    PointRange sampled_point_cloud;
    sampler(merged_point_cloud, SamplingOptions(0.01), sampled_point_cloud);

    int num_global_coordinates = sampled_point_cloud.size();
    using RangeQuery = typename gr::KdTree<Scalar>::template RangeQuery<>;
    const Scalar max_dist = 0.00001;
    const Scalar epsilon = 0.001;
    const Scalar sq_eps = epsilon*epsilon;
    for(size_t i = 0; i < num_point_clouds; i++){
        gr::KdTree<Scalar> tree = constructKdTree(transformed_point_clouds[i]);
        for(int j = 0; j < num_global_coordinates; j++){
            RangeQuery query;
            query.queryPoint = sampled_point_cloud[j].pos();
            query.sqdist     = sq_eps;

            auto result = tree.doQueryRestrictedClosestIndex( query );

            if ( result.first != gr::KdTree<Scalar>::invalidIndex() ) {
                if( result.second < max_dist )
                    patches[i].emplace_back(point_clouds[i][result.first], j);
            }
        }
    }

    return num_global_coordinates;
}


template <typename PointRange>
Scalar compute_lcp( const gr::KdTree<Scalar>& P, const PointRange& Q){
  using RangeQuery = typename gr::KdTree<Scalar>::template RangeQuery<>;
  const Scalar epsilon = 0.01;
  std::atomic_uint good_points(0);
  const size_t number_of_points = Q.size();
  Scalar best_LCP_ = 0;
  const size_t terminate_value = best_LCP_ * number_of_points;
  const Scalar sq_eps = epsilon*epsilon;

  for (size_t i = 0; i < number_of_points; ++i) {
    RangeQuery query;
    query.queryPoint = Q[i].pos();
    query.sqdist     = sq_eps;

    auto result = P.doQueryRestrictedClosestIndex( query );

    if ( result.first != gr::KdTree<Scalar>::invalidIndex() ) {
        good_points++;
    }

    // We can terminate if there is no longer chance to get better than the
    // current best LCP.
    if (number_of_points - i + good_points < terminate_value) {
        break;
    }
  }
  return Scalar(good_points) / Scalar(number_of_points);
}


int main(int argc, const char **argv) {
    if(argc < 2){
        std::cout << "execute program using: " << "./gret-sdp" << " <config/file/path>" << std::endl;
        return EXIT_FAILURE;
    }

    std::string config_file(argv[1]);

    using MatcherType = GRET_SDP<PointType, DummyTransformVisitor, GRET_SDP_Options>;
    using OptionType  = typename MatcherType::OptionsType;
    OptionType options;

    // read points and transformations from file
    vector<Transform> ground_truth_transformations;
    vector<PointRange> point_clouds;
    extractPCAndTrFromStandfordConfFile(config_file, ground_truth_transformations, point_clouds);

    // compute patches
    int num_point_clouds = point_clouds.size();
    vector<PwiRange> patches(num_point_clouds);
    int num_global_coordinates = computePatches(point_clouds, ground_truth_transformations, patches);

    MatcherType matcher(options, logger);

    DummyTransformVisitor tr_visitor;

    using WrapperType =
#ifdef OpenGR_USE_MOSEK
        MOSEK_WRAPPER<Scalar>;
#elif OpenGR_USE_SDPA
        SDPA_WRAPPER<Scalar>;
#else
        void;
#error Could not find a wrapper. Either use SDPA or MOSEK
#endif

    // register patches
    matcher.RegisterPatches< WrapperType >(patches, num_global_coordinates, tr_visitor);

    // get transformations and registered patches
    std::vector<MatrixType> transformations;
    matcher.getTransformations(transformations);
    std::vector<PointType> registered_patches;
    matcher.getRegisteredPatches(registered_patches);

    // verify results
    int num_total_points = std::accumulate(point_clouds.begin(), point_clouds.end(), 0, [](int sum, const auto& pc){ return sum + pc.size(); });
    PointRange registered_points(num_total_points);
    PointRange merged_points(num_total_points);
    // transform to common frame
    for (size_t i = 0; i < point_clouds.size(); i++){
        for(const auto& point : point_clouds[i]){
            MatrixType trafo =  ground_truth_transformations[0].inverse().matrix() * transformations[0].inverse() * transformations[i];
            registered_points.emplace_back(VectorType((trafo * point.pos().homogeneous()).template head<3>()));
            merged_points.emplace_back(ground_truth_transformations[i].inverse() * point.pos());
        }
    }


    gr::KdTree<Scalar> tree(constructKdTree(merged_points));

    // compute lcp
     std::cout << "compute lcp between registered points and ground truth registration" << std::endl;
    Scalar lcp = compute_lcp(tree, registered_points);
    std::cout << "lcp = " << lcp << std::endl;

#define WRITE_OUTPUT_FILES
    vector<Eigen::Matrix2f> tex_coords;
    vector<VectorType> normals;
    vector<tripple> tris;
    vector<string> mtls;
    stringstream iss;
    iss << "registered_point_clouds.ply";
    std::cout << "Exporting file " << iss.str().c_str() << "\n";
    IOManager iomanager;
    iomanager.WriteObject(iss.str().c_str(),
                           registered_points,
                           tex_coords,
                           normals,
                           tris,
                           mtls);
#define WRITE_OUTPUT_FILES

    return EXIT_SUCCESS;
}
