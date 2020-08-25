// Copyright 2017 Nicolas Mellado
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
// Authors: Dror Aiger, Yoni Weill, Nicolas Mellado
//
// This file is part of the implementation of the 4-points Congruent Sets (4PCS)
// algorithm presented in:
//
// 4-points Congruent Sets for Robust Surface Registration
// Dror Aiger, Niloy J. Mitra, Daniel Cohen-Or
// ACM SIGGRAPH 2008 and ACM Transaction of Graphics.
//
// Given two sets of points in 3-space, P and Q, the algorithm applies RANSAC
// in roughly O(n^2) time instead of O(n^3) for standard RANSAC, using an
// efficient method based on invariants, to find the set of all 4-points in Q
// that can be matched by rigid transformation to a given set of 4-points in P
// called a base. This avoids the need to examine all sets of 3-points in Q
// against any base of 3-points in P as in standard RANSAC.
// The algorithm can use colors and normals to speed-up the matching
// and to improve the quality. It can be easily extended to affine/similarity
// transformation but then the speed-up is smaller because of the large number
// of congruent sets. The algorithm can also limit the range of transformations
// when the application knows something on the initial pose but this is not
// necessary in general (though can speed the runtime significantly).

// Home page of the 4PCS project (containing the paper, presentations and a
// demo): http://graphics.stanford.edu/~niloy/research/fpcs/fpcs_sig_08.html
// Use google search on "4-points congruent sets" to see many related papers
// and applications.

#ifndef _OPENGR_ALGO_ABSTRACT_MATCH_BASE_
#define _OPENGR_ALGO_ABSTRACT_MATCH_BASE_

#include <vector>

#ifdef OpenGR_USE_OPENMP
#include <omp.h>
#endif

#include "gr/utils/shared.h"
#include "gr/utils/sampling.h"
#include "gr/utils/logger.h"
#include "gr/utils/crtp.h"

#ifdef TEST_GLOBAL_TIMINGS
#   include "gr/utils/timer.h"
#endif

namespace gr{

struct DummyTransformVisitor {
    template <typename Derived>
    inline void operator() (float, float, const Eigen::MatrixBase<Derived>&) const {}
    constexpr bool needsGlobalTransformation() const { return false; }
};

/// \brief Abstract class for registration algorithms
template <typename PointType, typename _TransformVisitor = DummyTransformVisitor,
          template < class, class > class ... OptExts>
class AbstractMatchBase {

public:
    using Scalar = typename PointType::Scalar;
    using VectorType = typename PointType::VectorType;
    using MatrixType = Eigen::Matrix<Scalar, 4, 4>;
    using LogLevel = Utils::LogLevel;
    using TransformVisitor = _TransformVisitor;

    template < class Derived, class TBase>
    class Options : public TBase
    {
    public:
        using Scalar = typename PointType::Scalar;

        /// The number of points in the sample. We sample this number of points
        /// uniformly from P and Q.
        size_t sample_size = 200;
        /// Maximum time we allow the computation to take. This makes the algorithm
        /// an ANY TIME algorithm that can be stopped at any time, producing the best
        /// solution so far.
        /// \warning Max. computation time must be handled in child classes
        int max_time_seconds = 60;
        /// use a constant default seed by default
        unsigned int randomSeed = std::mt19937::default_seed;
    };

    using OptionsType = gr::Utils::CRTP < OptExts ... , Options >;

    /// A convenience class used to wrap (any) PointType to allow mutation of position
    /// of point samples for internal computations.
    struct PosMutablePoint : public PointType
    {
        using VectorType = typename PointType::VectorType;

        private:
            VectorType posCopy;

        public:
            template<typename ExternalType>
            PosMutablePoint(const ExternalType& i)
                : PointType(i), posCopy(PointType(i).pos()) { }

            inline VectorType & pos() { return posCopy; }

            inline VectorType pos() const { return posCopy; }
    };

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline AbstractMatchBase(const OptionsType& options, const Utils::Logger &logger)
        : randomGenerator_(options.randomSeed), logger_(logger), options_(options) {}

    virtual inline ~AbstractMatchBase() {}


protected:
    std::mt19937 randomGenerator_;
    const Utils::Logger &logger_;

    OptionsType options_;

    /// \todo Rationnalize use and name of this variable
    static constexpr int kNumberOfDiameterTrials = 1000;

protected :
    template <Utils::LogLevel level, typename...Args>
    inline void Log(Args...args) const { logger_.Log<level>(args...); }


    /// Initializes the data structures and needed values before the match
    /// computation.
    /// This method is called once the internal state of the Base class as been
    /// set.
    virtual void Initialize() { }


}; /// class AbstractMatchBase
} /// namespace gr

#endif // _OPENGR_ALGO_ABSTRACT_MATCH_BASE_
