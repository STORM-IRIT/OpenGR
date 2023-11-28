#ifndef _OPENGR_ALGO_GRET_SDP_
#define _OPENGR_ALGO_GRET_SDP_

#include <vector>
#include <Eigen/Dense>

#include "gr/utils/shared.h"
#include "gr/algorithms/NAryMatchBase.h"

namespace gr{


template < class Derived, class TBase>
class GRET_SDP_Options : public  TBase
{
public:
    using Scalar = typename TBase::Scalar;

};

/// Class for the computation of the GRET-SDP algorithm as described in [this paper]((https://arxiv.org/abs/1306.5226)).
template <typename _PointType,
          typename _TransformVisitor,
          template < class, class > class ... OptExts >
class GRET_SDP : public NAryMatchBase<_PointType, _TransformVisitor, OptExts ..., GRET_SDP_Options> {

public:
    using TransformVisitor = _TransformVisitor;

    using MatchBaseType = NAryMatchBase<_PointType, _TransformVisitor, OptExts ..., GRET_SDP_Options>;
    using PosMutablePoint = typename MatchBaseType::PosMutablePoint;
    using OptionsType = typename MatchBaseType::OptionsType;

    using Scalar = typename MatchBaseType::Scalar;
    using VectorType = typename MatchBaseType::VectorType;
    using MatrixType = typename MatchBaseType::MatrixType;
    
    using MatrixX = typename Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

    using LogLevel = typename MatchBaseType::LogLevel;

    
    inline GRET_SDP(const OptionsType& options , const Utils::Logger &logger)
        : MatchBaseType(options, logger) {}

    inline virtual ~GRET_SDP() {}

    /// Computes the global coordinates of patches (subsets of the global coordinates)
    /// and its corresponding transformations.
    /// @param [in] patches is a range of point clouds with indexes assigned to each point. The 
    /// index of a point refers to its corresponding global coordinate.
    /// @param [in] n number of global coordinates.
    template<typename Solver, typename PatchRange>
    void RegisterPatches(const PatchRange& patches, const int n, TransformVisitor& v);

    /// returns the registered patches
    /// @param [out] registered_points Range of the registered patches. 
    template<typename PointRange>
    void getRegisteredPatches(PointRange& registered_points);

    /// returns the transformations
    /// @param [out] transformations Range of transformations that register the patches. 
    /// the i-th transformation corresponds to the i-th patch in order to register it.
    template<typename TrRange>
    void getTransformations(TrRange& transformations);

    private:
    int d; /// Dimension of points. 
    int m; /// Number of patches.
    int n; /// Number of global coordinates.

    /// 
    MatrixX O;
    MatrixX OB;
    MatrixX Linv;

}; /// class GRET_SDP
} /// namespace gr
#include "GRET_SDP.hpp"

#endif // _OPENGR_ALGO_GRET_SDP_
