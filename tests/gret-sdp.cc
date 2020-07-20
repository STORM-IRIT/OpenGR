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
//
// Super 4PCS: Fast Global Pointcloud Registration via Smart Indexing
// Nicolas Mellado, Dror Aiger, Niloy J. Mitra
// Symposium on Geometry Processing 2014.
//
// Data acquisition in large-scale scenes regularly involves accumulating
// information across multiple scans. A common approach is to locally align scan
// pairs using Iterative Closest Point (ICP) algorithm (or its variants), but
// requires static scenes and small motion between scan pairs. This prevents
// accumulating data across multiple scan sessions and/or different acquisition
// modalities (e.g., stereo, depth scans). Alternatively, one can use a global
// registration algorithm allowing scans to be in arbitrary initial poses. The
// state-of-the-art global registration algorithm, 4PCS, however has a quadratic
// time complexity in the number of data points. This vastly limits its
// applicability to acquisition of large environments. We present Super 4PCS for
// global pointcloud registration that is optimal, i.e., runs in linear time (in
// the number of data points) and is also output sensitive in the complexity of
// the alignment problem based on the (unknown) overlap across scan pairs.
// Technically, we map the algorithm as an ‘instance problem’ and solve it
// efficiently using a smart indexing data organization. The algorithm is
// simple, memory-efficient, and fast. We demonstrate that Super 4PCS results in
// significant speedup over alternative approaches and allows unstructured
// efficient acquisition of scenes at scales previously not possible. Complete
// source code and datasets are available for research use at
// http://geometry.cs.ucl.ac.uk/projects/2014/super4PCS/.

#include "gr/accelerators/gret_sdp/wrappers.h"
#include "gr/algorithms/GRET_SDP.h"
#include "gr/utils/timer.h"
#include "gr/accelerators/kdtree.h"


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
using PatchType = vector<pair<PointType, int>>;

constexpr Utils::LogLevel loglvl = Utils::Verbose;
Utils::Logger logger(loglvl);

struct RegistrationProblem {
    int n;
    int m;
    int d;
    vector<PatchType> patches;
};

template <typename PatchRange>
void readPatch(const string& file, PatchRange& patch){
    int num_points;
    fs::ifstream file_stream(file);
    file_stream >> num_points;
    patch.reserve(num_points);

    VectorType vec;
    int index;
    while(file_stream >> index){
        for(int i = 0; i < Dim; i++)
            file_stream >> vec(i);
        patch.emplace_back(PointType(vec), index);
  }
}

void readTransformation(const string& filePath, MatrixType& trafo){
    int rows, cols;
    fs::ifstream file(filePath);
    if(file.is_open()){
        file >> rows >> cols;
        if(cols != Dim+1 || rows != Dim+1)
            throw std::runtime_error("matrices have to be of size " + to_string(Dim+1) + "x" + to_string(Dim+1));
        for(int i = 0; i < cols; i++)
            for (int j = 0; j < rows; j++)
                file >> trafo(i, j);
    }
}

template <typename TrRange>
void extractPatchesAndTrFromConfigFile(const string& configFilePath,  RegistrationProblem& problem, TrRange& transformations){
    const string workingDir = fs::path(configFilePath).parent_path().native();

    pt::ptree root;
    pt::read_json(configFilePath, root);

    int n = root.get<int>("n");
    int m = root.get<int>("m");
    int d = root.get<int>("d");

    problem.n = n;
    problem.m = m;
    problem.d = d;

    vector< string  > patchFiles;

    for (pt::ptree::value_type &item : root.get_child("patches"))
    {
        patchFiles.push_back(item.second.data());
    }

    if(patchFiles.size() != m)
        throw runtime_error("Number of patches m and number of given patch files is not the same.");

    if(d != Dim)
        throw runtime_error("Dimension of point type has to be " + to_string(Dim));

    // read patch files
    problem.patches.resize(m);
    ifstream patch_file;
    for(int i = 0; i < m; i++){
        readPatch(workingDir + "/" + patchFiles[i], problem.patches[i]);
    }
    
    vector< string  > transformationFiles;
    for (pt::ptree::value_type &item : root.get_child("gt_trafos"))
        transformationFiles.push_back(item.second.data());

    if(transformationFiles.size() != m)
        throw runtime_error("Number of transformations and number of given transformation files is not the same.");

    transformations.reserve(m);
    MatrixType trafo;
    for(int i = 0; i < m; i++){
        readTransformation(workingDir + "/" + transformationFiles[i], trafo);
        transformations.emplace_back(trafo);
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
    if(argc != 2){
        std::cout << "execute program using: " << "./gret-sdp" << " <config/file/path>" << std::endl;
        return EXIT_FAILURE;
    }

    std::string config_file(argv[1]);


    using MatcherType = GRET_SDP<PointType, DummyTransformVisitor, GRET_SDP_Options>;
    using OptionType  = typename MatcherType::OptionsType;
    OptionType options;

    RegistrationProblem problem;

    vector<MatrixType> gt_transformations;
    extractPatchesAndTrFromConfigFile(config_file, problem, gt_transformations);

    const int d = problem.d;
    const int n = problem.n;
    const int m = problem.m;
    const vector<PatchType>& patches = problem.patches;

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

    matcher.RegisterPatches< WrapperType >(patches, n, tr_visitor);

    std::vector<MatrixType> transformations;
    matcher.getTransformations(transformations);
    std::vector<PointType> registered_patches;
    matcher.getRegisteredPatches(registered_patches);

    // verify results
    vector<PointType> reg_transformed_patches;
    vector<PointType> ori_transformed_patches;
    int accum_patch_size = 0;
    for(int i = 0; i < m; i++) accum_patch_size += patches[i].size();

    reg_transformed_patches.reserve(accum_patch_size);
    ori_transformed_patches.reserve(accum_patch_size);

    VectorType tmp;
    Eigen::Matrix3d reg_O_0(transformations[0].block(0,0,d,d).transpose());
    Eigen::Vector3d reg_t_0(transformations[0].block(0,d,d,1));
    Eigen::Matrix3d ori_O_0(gt_transformations[0].block(0,0,d,d).transpose());
    Eigen::Vector3d ori_t_0(gt_transformations[0].block(0,d,d,1));

    // computing the transformed patches for reference frame of patch one
    PointType point;
    int index;
    for(int i = 0; i < m; i++){
        for (const pair<PointType, int>& point_with_index : patches[i])
        {
            point = point_with_index.first;
            index = point_with_index.second;
            // using computed transformations
            tmp = (transformations[i]*point.pos().homogeneous()).template head<3>();
            tmp = reg_O_0 * (tmp - reg_t_0);
            reg_transformed_patches.emplace_back(PointType(tmp));

            // using ground truth transformations
            tmp = (gt_transformations[i]*point.pos().homogeneous()).template head<3>();
            tmp = ori_O_0 * (tmp - ori_t_0);
            ori_transformed_patches.emplace_back(PointType(tmp));
        }
        
    }

    // construct kd_tree
    gr::KdTree<Scalar> kd_tree(constructKdTree(ori_transformed_patches));
    // compute lcp
    Scalar lcp = compute_lcp(kd_tree, reg_transformed_patches);
    std::cout << "lcp = " << lcp << std::endl;

    return EXIT_SUCCESS;
}
