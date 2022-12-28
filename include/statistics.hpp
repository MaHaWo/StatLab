#ifndef STATLAB_STATISTICS_HPP
#define STATLAB_STATISTICS_HPP

#include <cmath>

namespace StatLab
{
namespace Statistics
{
template < typename T, std::size_t order, class BinaryFunc >
struct MomentOperator
{
};

template < typename T, std::size_t order, class BinaryFunc >
struct CentralMomentOperator
{
};

template < typename T, std::size_t order, class BinaryFunc >
struct CentralMomentOnlineOperator
{
};

template < typename T, class BinaryFunc >
struct HarmonicMeanOperator
{
};

template < typename T, class BinaryFunc >
struct GeometricMeanOperator
{
};

template < typename T, class BinaryFunc >
struct SampleVarianceOperator
{
};

template < typename T, class BinaryFunc >
struct SampleVarianceOnlineOperator
{
};

template < typename T, class BinaryFunc >
struct SampleStandardDeviationOperator
{};

template < typename T, class BinaryFunc >
struct SampleStandardDeviationOnline
{};

}
}
#endif