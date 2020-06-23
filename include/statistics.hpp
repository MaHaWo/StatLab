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
template < typename T, NaNPolicy nan_p = NaNPolicy::propagate >
class SumPairwise final
{
  private:
    T _s;

  public:
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
    template < typename U >
    inline constexpr T
    operator()(U&& value)
    {
        static_assert(std::is_convertible< U, T >::value,
                      "Error, type T in struct 'sum' not comptabile with type "
                      "'U' used by call operator");

        if constexpr (nan_p == NaNPolicy::propagate)
        {
            _s += value;
        }
        else if constexpr (nan_p == NaNPolicy::skip)
        {
            if (not std::isnan(value) and not std::isinf(value))
            {
                _s += value;
            }
        }
        else
        {
            if (not std::isnan(value) and not std::isinf(value))
            {
                _s += value;
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
template < typename T, NaNPolicy nan_p = NaNPolicy::propagate >
class SumKahan final
{
  private:
    T _y;
    T _t;
    T _s;
    T _comp;

  public:
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
            _y = value - _comp;
            _t = _s + _y;
            _s = _t;
        }
        else if constexpr (nan_p == NaNPolicy::skip)
        {
            if (not std::isnan(value) and not std::isinf(value))
            {
                _y = value - _comp;
                _t = _s + _y;
                _s = _t;
            }
        }
        else
        {
            if (not std::isnan(value) and not std::isinf(value))
            {
                _y = value - _comp;
                _t = _s + _y;
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
template < typename T, NaNPolicy nan_p = NaNPolicy::propagate >
class Product final
{
  private:
    T _result;

  public:
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
           template < typename... > class SumAlgo = SumPairwise,
           NaNPolicy nan_p                        = NaNPolicy::propagate >
class Moment final
{
  private:
    SumAlgo< T, nan_p > _sum;
    T                   _n;

  public:
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
    operator()(T v) -> T
    {
        ++_n;
        _sum(v)
    }

    template < typename Iter, typename Getter >
    inline auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter) -> T
    {
        for (; begin != end; ++begin)
        {
            ++_n;
            _sum(getter(*begin));
        }

        return std::pow(_sum.result(), order) / _n;
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
 * @tparam order
 * @tparam SumPairwise
 * @tparam NaNPolicy::propagate
 */
template < typename T,
           std::size_t order,
           template < typename... > class SumAlgo = SumPairwise,
           NaNPolicy nan_p                        = NaNPolicy::propagate >
class CentralMoment final
{
  private:
    Moment< T, order, SumAlgo, nan_p > _moment;

  public:
    inline auto
    operator()(T v) -> T
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

template < typename T,
           template < typename... > class SumAlgo = SumPairwise,
           NaNPolicy nan_p                        = NaNPolicy::propagate >
class ArithmeticMean final

{
  public:
    inline auto
    operator()(T v) -> T
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

template < typename T,
           template < typename... > class SumAlgo = SumPairwise,
           NaNPolicy nan_p                        = NaNPolicy::propagate >
class HarmonicMean final

{
  public:
    inline auto
    operator()(T v) -> T
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

template < typename T,
           template < typename... > class SumAlgo = SumPairwise,
           NaNPolicy nan_p                        = NaNPolicy::propagate >
class GeometricMean final

{
  public:
    inline auto
    operator()(T v) -> T
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

template < typename T,
           template < typename... > class SumAlgo = SumPairwise,
           NaNPolicy nan_p                        = NaNPolicy::propagate >
class Variance final

{
  public:
    inline auto
    operator()(T v) -> T
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

template < typename T,
           template < typename... > class SumAlgo = SumPairwise,
           NaNPolicy nan_p                        = NaNPolicy::propagate >
class Std final

{
  public:
    inline auto
    operator()(T v) -> T
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

template < typename T,
           template < typename... > class SumAlgo = SumPairwise,
           NaNPolicy nan_p                        = NaNPolicy::propagate >
class Skewness final

{
  public:
    inline auto
    operator()(T v) -> T
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

template < typename T,
           template < typename... > class SumAlgo = SumPairwise,
           NaNPolicy nan_p                        = NaNPolicy::propagate >
class Kurtosis final

{
  public:
    inline auto
    operator()(T v) -> T
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

template < typename T,
           template < typename... > class SumAlgo = SumPairwise,
           NaNPolicy nan_p                        = NaNPolicy::propagate >
class ExcessKurtosis final

{
  public:
    inline auto
    operator()(T v) -> T
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

template < typename T, NaNPolicy nan_p = NaNPolicy::propagate >
class Mode final

{
  public:
    inline auto
    operator()(T v) -> T
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

template < typename T,
           std::size_t percent,
           NaNPolicy   nan_p = NaNPolicy::propagate >
class Quantile final

{
  public:
    inline auto
    operator()(T v) -> T
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

template < typename T,
           std::size_t percent,
           NaNPolicy   nan_p = NaNPolicy::propagate >
using Median = Quantile< T, 50, nan_p >;

}

#endif