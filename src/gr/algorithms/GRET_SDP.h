#ifndef _OPENGR_ALGO_GRET_SDP_
#define _OPENGR_ALGO_GRET_SDP_

#include <vector>

#include "gr/utils/shared.h"
#include "gr/algorithms/NAryMatchBase.h"

namespace gr{


template < class Derived, class TBase>
class GRETSDPOptions : public  TBase
{
public:
    using Scalar = typename TBase::Scalar;

};


/// \brief Base class for Congruent Sec Exploration algorithms
/// \tparam _Traits Defines properties of the Base used to build the congruent set.
template <typename _PointType,
          typename _TransformVisitor,
          template < class, class > class ... OptExts >
class GRET_SDP : public NAryMatchBase<_PointType, _TransformVisitor, OptExts ..., GRETSDPOptions> {

public:
    using TransformVisitor = _TransformVisitor;

    using MatchBaseType = NAryMatchBase<_PointType, _TransformVisitor, OptExts ..., GRETSDPOptions>;
    using PosMutablePoint = typename MatchBaseType::PosMutablePoint;
    using OptionsType = typename MatchBaseType::OptionsType;

    using Scalar = typename MatchBaseType::Scalar;
    using VectorType = typename MatchBaseType::VectorType;
    using MatrixType = typename MatchBaseType::MatrixType;

    using LogLevel = typename MatchBaseType::LogLevel;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW


    inline GRET_SDP(const OptionsType& options , const Utils::Logger &logger)
        : MatchBaseType(options, logger) {}

    inline virtual ~GRET_SDP() {}

    template <typename InputRange, template<typename> class Sampler>
    inline
    void ComputeTransformations(const std::vector<InputRange>& P,
                                int N_,
                                int M_,
                                int d_,
                                Eigen::MatrixXd& L,
                                std::vector<MatrixType>& transformations_,
                                std::vector<VectorType>& translations_,
                                const Sampler<_PointType>& sampler,
                                TransformVisitor& v);

    private:
    /// dimension, currently set to 3
    int d;
    /// number of points
    int N;
    /// number of patches
    int M;

    template <typename InputRange>
    void ComputeMatricesBDC(const std::vector<InputRange>& P, 
                            Eigen::MatrixXd& L);

}; /// class GRET_SDP
} /// namespace gr
#include "GRET_SDP.hpp"

#endif // _OPENGR_ALGO_GRET_SDP_
