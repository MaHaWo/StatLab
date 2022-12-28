#ifndef STATLAB_HPP
#define STATLAB_HPP

#include "arithmetics.hpp"
#include "statistics.hpp"
#include "utils.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <type_traits>
#include <unordered_map>

#include <boost/hana.hpp>

namespace StatLab
{
/**
 * @brief TODO
 *
 * @tparam T
 * @tparam order
 * @tparam SumPairwise
 * @tparam NaNPolicy::propagate
 */
template < typename T,
           std::size_t      order,
           Utils::NaNPolicy nan_p = Utils::NaNPolicy::propagate,
           typename BinaryFunc    = std::plus< T >,
           typename RangeAlgo     = Arithmetics::naive_iteration >
using Moment = Arithmetics::ArithmeticOperation<
    T,
    nan_p,
    Statistics::MomentOperator< T, order, BinaryFunc >,
    RangeAlgo >;

/**
 * @brief
 *
 * @tparam T
 * @tparam order
 * @tparam nan_p
 * @tparam BinaryFunc
 * @tparam RangeAlgo
 */
template < typename T,
           std::size_t      order,
           Utils::NaNPolicy nan_p = Utils::NaNPolicy::propagate,
           typename BinaryFunc    = std::plus< T >,
           typename RangeAlgo     = Arithmetics::naive_iteration >
using CentralMoment = Arithmetics::ArithmeticOperation<
    T,
    nan_p,
    Statistics::CentralMomentOperator< T, order, BinaryFunc >,
    RangeAlgo >;

/**
 * @brief
 *
 * @tparam T
 * @tparam order
 * @tparam nan_p
 * @tparam BinaryFunc
 * @tparam RangeAlgo
 */
template < typename T,
           std::size_t      order,
           Utils::NaNPolicy nan_p = Utils::NaNPolicy::propagate,
           typename BinaryFunc    = std::plus< T >,
           typename RangeAlgo     = Arithmetics::naive_iteration >
using CentralMomentOnline = Arithmetics::ArithmeticOperation<
    T,
    nan_p,
    Statistics::CentralMomentOnlineOperator< T, order, BinaryFunc >,
    RangeAlgo >;

/**
 * @brief
 *
 * @tparam T
 * @tparam nan_p
 * @tparam BinaryFunc
 * @tparam RangeAlgo
 */
template < typename T,
           Utils::NaNPolicy nan_p = Utils::NaNPolicy::propagate,
           typename BinaryFunc    = std::plus< T >,
           typename RangeAlgo     = Arithmetics::naive_iteration >
using ArithmeticMean = Moment< T, 1, nan_p, BinaryFunc, RangeAlgo >;

/**
 * @brief
 *
 * @tparam T
 * @tparam nan_p
 * @tparam BinaryFunc
 * @tparam RangeAlgo
 */
template < typename T,
           Utils::NaNPolicy nan_p = Utils::NaNPolicy::propagate,
           typename BinaryFunc    = std::plus< T >,
           typename RangeAlgo     = Arithmetics::naive_iteration >
using HarmonicMean = Arithmetics::ArithmeticOperation<
    T,
    nan_p,
    Statistics::HarmonicMeanOperator< T, BinaryFunc >,
    RangeAlgo >;

/**
 * @brief
 *
 * @tparam T
 * @tparam nan_p
 * @tparam BinaryFunc
 * @tparam RangeAlgo
 */
template < typename T,
           Utils::NaNPolicy nan_p = Utils::NaNPolicy::propagate,
           typename BinaryFunc    = std::plus< T >,
           typename RangeAlgo     = Arithmetics::naive_iteration >
using GeometricMean = Arithmetics::ArithmeticOperation<
    T,
    nan_p,
    Statistics::GeometricMeanOperator< T, BinaryFunc >,
    RangeAlgo >;

/**
 * @brief
 *
 * @tparam T
 * @tparam nan_p
 * @tparam BinaryFunc
 * @tparam RangeAlgo
 */
template < typename T,
           Utils::NaNPolicy nan_p = Utils::NaNPolicy::propagate,
           typename BinaryFunc    = std::plus< T >,
           typename RangeAlgo     = Arithmetics::naive_iteration >
using SampleVariance = Arithmetics::ArithmeticOperation<
    T,
    nan_p,
    Statistics::SampleVarianceOperator< T, BinaryFunc >,
    RangeAlgo >;

/**
 * @brief
 *
 * @tparam T
 * @tparam nan_p
 * @tparam BinaryFunc
 * @tparam RangeAlgo
 */
template < typename T,
           Utils::NaNPolicy nan_p = Utils::NaNPolicy::propagate,
           typename BinaryFunc    = std::plus< T >,
           typename RangeAlgo     = Arithmetics::naive_iteration >
using SampleVarianceOnline = Arithmetics::ArithmeticOperation<
    T,
    nan_p,
    Statistics::SampleVarianceOnlineOperator< T, BinaryFunc >,
    RangeAlgo >;

/**
 * @brief
 *
 * @tparam T
 * @tparam nan_p
 * @tparam BinaryFunc
 * @tparam RangeAlgo
 */
template < typename T,
           Utils::NaNPolicy nan_p = Utils::NaNPolicy::propagate,
           typename BinaryFunc    = std::plus< T >,
           typename RangeAlgo     = Arithmetics::naive_iteration >
using SampleStandardDeviation = Arithmetics::ArithmeticOperation<
    T,
    nan_p,
    Statistics::SampleStandardDeviationOperator< T, BinaryFunc >,
    RangeAlgo >;

/**
 * @brief
 *
 * @tparam T
 * @tparam nan_p
 * @tparam BinaryFunc
 * @tparam RangeAlgo
 */
template < typename T,
           Utils::NaNPolicy nan_p = Utils::NaNPolicy::propagate,
           typename BinaryFunc    = std::plus< T >,
           typename RangeAlgo     = Arithmetics::naive_iteration >
using SampleStandardDeviationOnline = Arithmetics::ArithmeticOperation<
    T,
    nan_p,
    Statistics::SampleStandardDeviationOnline< T, BinaryFunc >,
    RangeAlgo >;

/**
 * @brief TODO
 *
 * @tparam T
 * @tparam NaNPolicy::propagate
 * @tparam std::hash
 * @tparam std::equal_to
 */
template < typename T,
           template < typename, NaNPolicy, template < typename... > class >
           class SumAlgo                          = SumPairwise,
           template < typename... > class Hasher  = std::hash,
           template < typename... > class Compare = std::equal_to >
class Mode
{
    std::unordered_map< T, std::size_t, Hasher< T >, Compare< T > > _counts;

  public:
    using result_type = T;

    /**
     * @brief TODO
     *
     */
    inline auto
    reset() -> void
    {
        _counts.clear();
    }

    /**
     * @brief TODO
     *
     * @return T
     */
    inline auto
    result() -> T
    {
        return std::max_element(_counts.begin(),
                                _counts.end(),
                                [](auto&& lhs, auto&& rhs)
                                { return lhs.second < rhs.second; })
            ->second;
    }

    /**
     * @brief TODO
     *
     * @param v
     * @return T
     */
    inline auto
    operator()(const T& v) -> T
    {
        ++_counts[v];

        return result();
    }

    /**
     * @brief TODO
     *
     * @tparam Iter
     * @tparam Getter
     * @param begin
     * @param end
     * @param getter
     * @return T
     */
    template < typename Iter, typename Getter >
    inline auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter) -> T
    {
        for (; begin != end; ++begin)
        {
            ++_counts[getter(*begin)];
        }

        return result();
    }

    /**
     * @brief TODO
     *
     * @tparam Iter
     * @param begin
     * @param end
     * @return T
     */
    template < typename Iter >
    inline auto
    operator()(Iter&& begin, Iter&& end) -> T
    {
        return operator()(std::forward< Iter >(begin),
                          std::forward< Iter >(end),
                          [](auto&& v) { return v; });
    }
};

/**
 * @brief TODO
 *
 * @tparam T
 * @tparam percent
 * @tparam NaNPolicy::propagate
 */
template < typename T,
           std::size_t percent,
           NaNPolicy   nan_p = NaNPolicy::propagate >
class Quantile

{
  public:
    using result_type = T;

    /**
     * @brief TODO
     *
     * @param v
     * @return T
     */
    inline auto
    operator()(const T& v) -> T
    {
        // static_assert(false, "Not implemented yet");
    }

    /**
     * @brief TODO
     *
     * @tparam Iter
     * @tparam Getter
     * @param begin
     * @param end
     * @param getter
     * @return T
     */
    template < typename Iter, typename Getter >
    inline auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter) -> T
    {
        // static_assert(false, "Not implemented yet");
    }

    /**
     * @brief TODO
     *
     * @tparam Iter
     * @param begin
     * @param end
     * @return T
     */
    template < typename Iter >
    inline auto
    operator()(Iter&& begin, Iter&& end) -> T
    {
        operator()(std::forward< Iter >(begin),
                   std::forward< Iter >(end),
                   [](auto&& v) { return v; });
    }
};

/**
 * @brief TODO
 *
 * @tparam T
 * @tparam percent
 * @tparam NaNPolicy::propagate
 */
template < typename T,
           std::size_t percent,
           NaNPolicy   nan_p = NaNPolicy::propagate >
using Median = Quantile< T, 50, nan_p >;

template < typename T,
           NaNPolicy nan_p                           = NaNPolicy::propagate,
           template < typename... > class BinaryFunc = std::plus >
class Describe
{
  public:
    struct Result
    {
        T count;
        T mean;
        T std;
        T skewness;
        T kurtosis;
        T min;
        T q25;
        T median;
        T q75;
        T max;
    };

    using result_type = Result;

    inline auto
    operator()(const T& v) -> T
    {
        // // static_assert(false, "Not implemented yet");
    }

    template < typename Iter, typename Getter >
    inline auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter) -> result_type
    {
        // static_assert(false, "Not implemented yet");
    }

    template < typename Iter >
    inline auto
    operator()(Iter&& begin, Iter&& end) -> result_type
    {
        operator()(std::forward< Iter >(begin),
                   std::forward< Iter >(end),
                   [](auto&& v) { return v; });
    }
};

template < typename T, NaNPolicy nan_p = NaNPolicy::propagate >
class TukeyFiveNumbers
{

  public:
    struct Result
    {
        T count;
        T mean;
        T std;
        T skewness;
        T kurtosis;
        T min;
        T q25;
        T median;
        T q75;
        T max;
    };

    using result_type = Result;

    inline auto
    operator()(const T& v)
    {
        // static_assert(false, "Not implemented yet");
    }

    template < typename Iter, typename Getter >
    inline auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter) -> result_type
    {
        // static_assert(false, "Not implemented yet");
    }

    template < typename Iter >
    inline auto
    operator()(Iter&& begin, Iter&& end) -> result_type
    {
        operator()(std::forward< Iter >(begin),
                   std::forward< Iter >(end),
                   [](auto&& v) { return v; });
    }
};

template < typename T,
           NaNPolicy nan_p,
           template < typename... >
           class BinaryFunc,
           typename... Operators >
class describe
{

    boost::hana::tuple< Operators... > _operators;

    boost::hana::tuple< typename Operators::result_type... > _results;

  public:
    using result_type =
        boost::hana::tuple< typename Operators::result_type... >;

    inline auto
    operator()(const T& v) -> T
    {
        // // static_assert(false, "Not implemented yet");
    }

    template < typename Iter, typename Getter >
    inline auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter) -> result_type
    {
        // static_assert(false, "Not implemented yet");
    }

    template < typename Iter >
    inline auto
    operator()(Iter&& begin, Iter&& end) -> result_type
    {
        operator()(std::forward< Iter >(begin),
                   std::forward< Iter >(end),
                   [](auto&& v) { return v; });
    }
};

}

#endif