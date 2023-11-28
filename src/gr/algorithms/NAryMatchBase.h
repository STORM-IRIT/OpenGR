
#ifndef _OPENGR_ALGO_N_ARY_MATCH_BASE_
#define _OPENGR_ALGO_N_ARY_MATCH_BASE_

#include <vector>
#include "gr/algorithms/AbstractMatchBase.h"
#include "gr/utils/shared.h"
#include "gr/utils/logger.h"

namespace gr{

template < class Derived, class TBase>
    class NAryMatchOptions : public TBase
    {
    public:
        using Scalar = typename TBase::Scalar;

        // n-ary matching options here
    };


/// \brief Abstract class for n-ary registration algorithms (multiple point clouds)
template   <typename PointType,
            typename _TransformVisitor,
            template < class, class > class ... OptExts>
class NAryMatchBase : public AbstractMatchBase<PointType, _TransformVisitor, OptExts ..., NAryMatchOptions> {

public:
    using MatchBaseType = AbstractMatchBase<PointType, _TransformVisitor, OptExts ..., NAryMatchOptions>;

    using Scalar = typename MatchBaseType::Scalar;
    using VectorType = typename MatchBaseType::VectorType;
    using MatrixType = typename MatchBaseType::MatrixType;
    using LogLevel = typename MatchBaseType::LogLevel;
    using OptionsType = typename MatchBaseType::OptionsType;
    using PosMutablePoint = typename MatchBaseType::PosMutablePoint;
    using TransformVisitor = _TransformVisitor;


    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline NAryMatchBase(const OptionsType& options, const Utils::Logger &logger)
        : MatchBaseType(options, logger) {}

    virtual inline ~NAryMatchBase() {}

    protected:

}; /// class NAryMatchBase
} /// namespace gr

#endif // _OPENGR_ALGO_N_ARY_MATCH_BASE_
