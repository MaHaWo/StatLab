/** @file statistics.hpp
 *  @brief Implements functions for computing statistics.
 *
 *  @bug Quantiles return only exact quantiles, that is only
 *       values which are contained in the container.
 *  @todo Implement interpolation for quantiles, since they are defined such
 *
 */
#ifndef UTOPIA_MODELS_AMEEMULTI_STATISTICS_HH
#define UTOPIA_MODELS_AMEEMULTI_STATISTICS_HH

#include "utils.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <type_traits>
#include <unistd.h>
#include <unordered_map>
#include <boost/hana.hpp>

using namespace Statistics::Utils;

namespace Statistics
{

/**
 * @brief Indicate how NaN or Inf shall be treated:
 *  - propagate: default, this propagates nan or inf to the result
 *  - skip: skip them
 *  - error: throw invalid_agument exception when nan or inf is encountered
 */
enum class NaNPolicy
{
    propagate,
    skip,
    error
};

/**
 * @brief TODO
 *
 * @tparam T
 * @tparam NaNPolicy::propagate
 */
template < typename T,
           NaNPolicy nan_p                           = NaNPolicy::propagate,
           template < typename... > class BinaryFunc = std::plus >
class SumPairwise final
{
  private:
    T               _s;
    BinaryFunc< T > _plus;

  public:
    using result_type = T;
    /**
     * @brief Construct a new Sum object
     *
     * @tparam U
     * @param v
     */
    template < class U >
    SumPairwise(U&& v) noexcept : _s(std::forward< U >(v))
    {
    }

    /**
     * @brief Construct a new Sum object
     *
     */
    SumPairwise() noexcept = default;

    /**
     * @brief Construct a new Sum object
     *
     * @param other
     */
    SumPairwise(const SumPairwise& other) noexcept = default;

    /**
     * @brief Construct a new Sum object
     *
     * @param other
     */
    SumPairwise(SumPairwise&& other) noexcept = default;

    /**
     * @brief TODO
     *
     * @param other
     * @return Sum&
     */
    SumPairwise&
    operator=(const SumPairwise& other) noexcept = default;

    /**
     * @brief TODO
     *
     * @param other
     * @return Sum&
     */
    SumPairwise&
    operator=(SumPairwise&& other) noexcept = default;

    /**
     * @brief Destroy the Sum object
     *
     */
    ~SumPairwise() noexcept = default;

    inline constexpr void
    reset() noexcept
    {
        _s = T{};
    }

    /**
     * @brief TODO
     *
     * @return constexpr T
     */
    inline constexpr T
    result() noexcept
    {
        return _s;
    }

    /**
     * @brief TODO
     *
     * @tparam U
     * @param value
     * @return constexpr T
     */
    inline constexpr T
    operator()(const T& value)
    {
        static_assert(std::is_convertible< U, T >::value,
                      "Error, type T in struct 'sum' not comptabile with type "
                      "'U' used by call operator");

        if constexpr (nan_p == NaNPolicy::propagate)
        {
            _s = _plus(_s, value);
        }
        else if constexpr (nan_p == NaNPolicy::skip)
        {
            if (not std::isnan(value) and not std::isinf(value))
            {
                _s = _plus(_s, value);
            }
        }
        else
        {
            if (not std::isnan(value) and not std::isinf(value))
            {
                _s = _plus(_s, value);
            }
            else
            {
                throw std::invalid_argument(
                    "Error: invalid value found in sum: " +
                    std::to_string(value));
            }
        }

        return _s;
    }

    /**
     * @brief TODO
     *
     * @tparam InputIterator
     * @tparam Getter
     * @param begin
     * @param end
     * @param getter
     * @return constexpr T
     */
    template < typename InputIterator, typename Getter >
    inline constexpr T
    operator()(InputIterator begin, InputIterator end, Getter&& getter)
    {
        if (begin != end)
        {
            if (std::distance(begin, end) <= 1000)
            {
                for (; begin != end; ++begin)
                {
                    operator()(getter(*begin));
                }
            }
            else
            {
                std::size_t t =
                    std::floor(double(std::distance(begin, end)) / 2.);

                this->operator()(
                    begin, std::next(begin, t), std::forward< Getter >(getter));

                this->operator()(
                    std::next(begin, t), end, std::forward< Getter >(getter));
            }
        }
        return _s;
    }

    /**
     * @brief TODO
     *
     * @tparam InputIterator
     * @param begin
     * @param end
     * @return constexpr T
     */
    template < typename InputIterator >
    inline constexpr T
    operator()(InputIterator begin, InputIterator end)
    {
        if (begin != end)
        {
            this->operator()(begin, end, [](auto&& value) { return value; });
        }
        return _s;
    }
};

/**
 * @brief TODO
 *
 * @tparam T
 * @tparam NaNPolicy::propagate
 */
template < typename T,
           NaNPolicy nan_p                           = NaNPolicy::propagate,
           template < typename... > class BinaryFunc = std::plus >
class SumKahan final
{
  private:
    T _y;
    T _t;
    T _s;
    T _comp;

    BinaryFunc< T > _plus;

  public:
    using result_type = T;

    /**
     * @brief Construct a new Sum Kahan object
     *
     * @tparam U
     * @param v
     */
    template < class U >
    SumKahan(U&& v) noexcept : _s(std::forward< U >(v))
    {
    }

    /**
     * @brief Construct a new Sum Kahan object
     *
     */
    SumKahan() noexcept = default;

    /**
     * @brief Construct a new Sum Kahan object
     *
     * @param other
     */
    SumKahan(const SumKahan& other) noexcept = default;

    /**
     * @brief Construct a new Sum Kahan object
     *
     * @param other
     */
    SumKahan(SumKahan&& other) noexcept = default;

    /**
     * @brief TODO
     *
     * @param other
     * @return SumKahan&
     */
    SumKahan&
    operator=(const SumKahan& other) noexcept = default;

    /**
     * @brief TODO
     *
     * @param other
     * @return SumKahan&
     */
    SumKahan&
    operator=(SumKahan&& other) noexcept = default;

    /**
     * @brief Destroy the Sum Kahan object
     *
     */
    ~SumKahan() noexcept = default;

    /**
     * @brief Reset internal state of the object to T{}
     *
     */
    inline constexpr void
    reset() noexcept
    {
        _s    = T{};
        _y    = T{};
        _t    = T{};
        _comp = T{};
    }

    /**
     * @brief Get the result of the operations carried out until now.
     *
     * @return constexpr T
     */
    inline constexpr T
    result() noexcept
    {
        return _s;
    }

    /**
     * @brief TODO
     *
     * @tparam InputIterator
     * @tparam Getter
     * @param begin
     * @param end
     * @param getter
     * @return constexpr T
     */
    template < typename InputIterator, typename Getter >
    inline constexpr T
    operator()(InputIterator&& begin, InputIterator&& end, Getter&& getter)
    {
        if (begin != end)
        {
            _comp = 0.;
            for (; begin != end; ++begin)
            {
                operator()(getter(*begin));
            }
        }
        return _s;
    }

    /**
     * @brief TODO
     *
     * @tparam InputIterator
     * @param begin
     * @param end
     * @return constexpr T
     */
    template < typename InputIterator >
    inline constexpr T
    operator()(InputIterator&& begin, InputIterator&& end)
    {
        static_assert(std::is_convertible_v<
                      T,
                      std::iterator_traits< InputIterator >::value_type >);

        return this->operator()(std::forward< InputIterator >(begin),
                                std::forward< InputIterator >(end),
                                [](auto&& value) { return value; });
    }

    /**
     * @brief TODO
     *
     * @tparam U
     * @param value
     * @return constexpr T
     */
    template < typename U >
    inline constexpr T
    operator()(U&& value)
    {
        static_assert(std::is_convertible_v< T, std::decay_t< U > >);

        if constexpr (nan_p == NaNPolicy::propagate)
        {
            _y = _plus(value, T(-1.) * _comp);
            _t = _plus(_s, _y);
            _s = _t;
        }
        else if constexpr (nan_p == NaNPolicy::skip)
        {
            if (not std::isnan(value) and not std::isinf(value))
            {
                _y = _plus(value, T(-1.) * _comp);
                _t = _plus(_s, _y);
                _s = _t;
            }
        }
        else
        {
            if (not std::isnan(value) and not std::isinf(value))
            {
                _y = _plus(value, T(-1.) * _comp);
                _t = _plus(_s, _y);
                _s = _t;
            }
            else
            {
                throw std::invalid_argument(
                    "Error: invalid value found in sum: " +
                    std::to_string(value));
            }
        }

        return _s;
    }
};

/**
 * @brief TODO
 *
 * @tparam T
 * @tparam NaNPolicy::propagate
 */
template < typename T,
           NaNPolicy nan_p                           = NaNPolicy::propagate,
           template < typename... > class BinaryFunc = std::plus >
class Product final
{
  private:
    T _result;

  public:
    using result_type = T;
};

/**
 * @brief TODO
 *
 * @tparam T
 * @tparam order
 * @tparam SumPairwise
 * @tparam NaNPolicy::propagate
 */
template < typename T,
           std::size_t order,
           template < typename... > class SumAlgo    = SumPairwise,
           NaNPolicy nan_p                           = NaNPolicy::propagate,
           template < typename... > class BinaryFunc = std::plus >
class Moment final
{
  private:
    SumAlgo< T, nan_p, BinaryFunc > _sum;
    T                               _n;

  public:
    using result_type = T;

    Moment() noexcept : _sum(), _n(T{})
    {
    }
    Moment(const Moment&) = default;
    Moment(Moment&&)      = default;
    Moment&
    operator=(const Moment&) = default;
    Moment&
    operator=(Moment&&) = default;
    ~Moment()           = default;

    inline auto
    operator()(const T& v) -> T
    {
        ++_n;
        _sum(std::pow(v, order));

        return _sum.result() / _n;
    }

    template < typename Iter, typename Getter >
    inline auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter) -> T
    {
        for (; begin != end; ++begin)
        {
            ++_n;
            _sum(std::pow(getter(*begin), order));
        }

        return _sum.result() / _n;
    }

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
 * @tparam order
 * @tparam SumPairwise
 * @tparam NaNPolicy::propagate
 */
template < typename T,
           std::size_t order,
           template < typename... > class SumAlgo    = SumPairwise,
           NaNPolicy nan_p                           = NaNPolicy::propagate,
           template < typename... > class BinaryFunc = std::plus >
class CentralMoment final
{
  private:
    // todo

  public:
    using result_type = T;

    inline auto
    operator()(const T& v) -> T
    {
        static_assert(false, "Not implemented yet");
    }

    template < typename Iter, typename Getter >
    inline auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter) -> T
    {
        static_assert(false, "Not implemented yet");
    }

    template < typename Iter >
    inline auto
    operator()(Iter&& begin, Iter&& end) -> T
    {
        return operator()(std::forward< Iter >(begin),
                          std::forward< Iter >(end),
                          [](auto&& v) { return v; });
    }
};

template < typename T,
           template < typename... > class SumAlgo    = SumPairwise,
           NaNPolicy nan_p                           = NaNPolicy::propagate,
           template < typename... > class BinaryFunc = std::plus >
using ArithmeticMean = Moment< T, 1, SumAlgo, nan_p, BinaryFunc >;

template < typename T,
           template < typename... > class SumAlgo    = SumPairwise,
           NaNPolicy nan_p                           = NaNPolicy::propagate,
           template < typename... > class BinaryFunc = std::plus >
class HarmonicMean final
{
  public:
    inline auto
    operator()(const T& v) -> T
    {
        static_assert(false, "Not implemented yet");
    }

    template < typename Iter, typename Getter >
    inline auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter) -> T
    {
        static_assert(false, "Not implemented yet");
    }

    template < typename Iter >
    inline auto
    operator()(Iter&& begin, Iter&& end) -> T
    {
        return operator()(std::forward< Iter >(begin),
                          std::forward< Iter >(end),
                          [](auto&& v) { return v; });
    }
};

template < typename T,
           template < typename... > class SumAlgo    = SumPairwise,
           NaNPolicy nan_p                           = NaNPolicy::propagate,
           template < typename... > class BinaryFunc = std::plus >
class GeometricMean final

{
  public:
    using result_type = T;

    inline auto
    operator()(const T& v) -> T
    {
        static_assert(false, "Not implemented yet");
    }

    template < typename Iter, typename Getter >
    inline auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter) -> T
    {
        static_assert(false, "Not implemented yet");
    }

    template < typename Iter >
    inline auto
    operator()(Iter&& begin, Iter&& end) -> T
    {
        return operator()(std::forward< Iter >(begin),
                          std::forward< Iter >(end),
                          [](auto&& v) { return v; });
    }
};

template < typename T,
           template < typename... > class SumAlgo    = SumPairwise,
           NaNPolicy nan_p                           = NaNPolicy::propagate,
           template < typename... > class BinaryFunc = std::plus >
class Variance final

{
  public:
    using result_type = T;

    inline auto
    operator()(const T& v) -> T
    {
        static_assert(false, "Not implemented yet");
    }

    template < typename Iter, typename Getter >
    inline auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter) -> T
    {
        static_assert(false, "Not implemented yet");
    }

    template < typename Iter >
    inline auto
    operator()(Iter&& begin, Iter&& end) -> T
    {
        return operator()(std::forward< Iter >(begin),
                          std::forward< Iter >(end),
                          [](auto&& v) { return v; });
    }
};

template < typename T,
           template < typename... > class SumAlgo    = SumPairwise,
           NaNPolicy nan_p                           = NaNPolicy::propagate,
           template < typename... > class BinaryFunc = std::plus >
class Std final

{
    Variance< T, SumAlgo, nan_p > _var;

  public:
    using result_type = T;

    inline auto
    operator()(const T& v) -> T
    {
        _var(v);
        return std::sqrt(_var.result());
    }

    template < typename Iter, typename Getter >
    inline auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter) -> T
    {
        _var(std::forward< Iter >(begin),
             std::forward< Iter >(end),
             std::forward< Getter >(getter));
        return std::sqrt(_var.result());
    }

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
 * @tparam SumPairwise 
 * @tparam NaNPolicy::propagate 
 * @tparam std::plus 
 */
template < typename T,
           template < typename... > class SumAlgo    = SumPairwise,
           NaNPolicy nan_p                           = NaNPolicy::propagate,
           template < typename... > class BinaryFunc = std::plus >
class Skewness final

{
  public:
    using result_type = T;

    inline constexpr auto reset() -> void 
    {

    }

    inline constexpr auto result() -> T 
    {

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
        static_assert(false, "Not implemented yet");
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
        static_assert(false, "Not implemented yet");
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
 * @tparam SumPairwise
 * @tparam
 * @tparam std::plus
 */
template < typename T,
           template < typename... > class SumAlgo    = SumPairwise,
           NaNPolicy nan_p                           = NaNPolicy::propagate,
           template < typename... > class BinaryFunc = std::plus >
class Kurtosis final

{
  public:
    using result_type = T;

    inline auto
    operator()(const T& v) -> T
    {
        static_assert(false, "Not implemented yet");
    }

    template < typename Iter, typename Getter >
    inline auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter) -> T
    {
        static_assert(false, "Not implemented yet");
    }

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
 * @tparam SumPairwise
 * @tparam NaNPolicy::propagate
 * @tparam std::plus
 */
template < typename T,
           template < typename... > class SumAlgo    = SumPairwise,
           NaNPolicy nan_p                           = NaNPolicy::propagate,
           template < typename... > class BinaryFunc = std::plus >
class ExcessKurtosis final

{
    
    Kurtosis< T, SumAlgo, nan_p, BinaryFunc > _kurtosis public
        : using result_type = T;

public:

    inline auto reset() -> void 
    {   
        _kurtosis.reset();
    }


    inline auto result() -> T 
    {
        return _kurtosis.result();
    }

    inline auto
    operator()(const T& v) -> T
    {
        _kurtosis(v);
        return _kurtosis.result() - 3.;
    }

    template < typename Iter, typename Getter >
    inline auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter) -> T
    {
        _kurtosis(std::forward< Iter >(begin),
                  std::forward< Iter >(end),
                  std::forward< Getter >(getter));
        return _kurtosis.result() - 3.;
    }

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
 * @tparam NaNPolicy::propagate 
 * @tparam std::hash 
 * @tparam std::equal_to 
 */
template < typename T,
           NaNPolicy nan_p                        = NaNPolicy::propagate,
           template < typename... > class Hasher  = std::hash,
           template < typename... > class Compare = std::equal_to >
class Mode final
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
                                [](auto&& lhs, auto&& rhs) {
                                    return lhs.second < rhs.second;
                                })
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
        ++_counts[T];

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
        for(; begin != end; ++begin)
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
class Quantile final

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
        static_assert(false, "Not implemented yet");
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
        static_assert(false, "Not implemented yet");
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
    operator()(const T& v)
    {
        static_assert(false, "Not implemented yet");
    }

    template < typename Iter, typename Getter >
    inline auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter) -> result_type
    {
        static_assert(false, "Not implemented yet");
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
        static_assert(false, "Not implemented yet");
    }

    template < typename Iter, typename Getter >
    inline auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter) -> result_type
    {
        static_assert(false, "Not implemented yet");
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

template < NaNPolicy nan_p,
           ,
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

    template < typename T >
    inline auto
    opertator()(T&& v) -> T
    {
        static_assert(false, "Not implemented yet");
    }

    template < typename Iter, typename Getter >
    inline auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter) -> result_type
    {
        static_assert(false, "Not implemented yet");
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