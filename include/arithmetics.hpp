#ifndef STATLAB_ARITHMETICS_HPP
#define STATLAB_ARITHMETICS_HPP

#include <algorithm>
#include <cmath>
#include <concepts>
#include <functional>
#include <ranges>
#include <type_traits>

#include "concepts.hpp"
#include "utils.hpp"

namespace StatLab
{
namespace Arithmetics
{

/**
 * @brief General skeleton for arithmetic operations.
 *
 * @tparam T Type to process
 * @tparam nan_policy Policy on how to deal with non-finite values
 * @tparam BinaryFunc Function to merge two values into one, e.g. std::plus
 */
template < typename T,
           StatLab::Utils::NaNPolicy nan_policy,
           class BinaryFunc,
           class RangeAlgo >
class ArithmeticOperation final
{

    static_assert(StatLab::Concepts::BinaryCallable< BinaryFunc, T, T >);

  protected:
    T               _init;
    BinaryFunc _func;

  public:
    /**
     * @brief TODO
     *
     * @param init
     * @return auto
     */
    auto
    reset(T init = T{})
    {
        _init = init;
        _func = BinaryFunc();
    }

    /**
     * @brief get result attained until now
     * 
     * @return auto 
     */
    auto result() -> T
    {
        return _init;
    }

    /**
     * @brief Construct a new Arithmetic Operation object
     *
     */
    ArithmeticOperation() noexcept = default;

    /**
     * @brief Construct a new Arithmetic Operation object
     *
     */
    ArithmeticOperation(const ArithmeticOperation&) noexcept = default;

    /**
     * @brief Construct a new Arithmetic Operation object
     *
     */
    ArithmeticOperation(ArithmeticOperation&&) noexcept = default;

    /**
     * @brief Move-construct an ArithmeticOperation object
     *
     * @return ArithmeticOperation&
     */
    ArithmeticOperation&
    operator=(ArithmeticOperation&&) noexcept = default;

    /**
     * @brief Copy on ArithmeticOperation object
     *
     * @return ArithmeticOperation&
     */
    ArithmeticOperation&
    operator=(const ArithmeticOperation&) noexcept = default;

    /**
     * @brief Destroy the Arithmetic Operation object
     *
     */
    ~ArithmeticOperation() noexcept = default;

    /**
     * @brief Construct a new Arithmetic Operation object
     *
     * @param init
     * @param func
     */
    ArithmeticOperation(T init, BinaryFunc func) : _init(init), _func(func)
    {
    }

    /**
     * @brief Construct a new Arithmetic Operation object
     *
     * @param func
     */
    ArithmeticOperation(BinaryFunc func) : _init(T{}), _func(func)
    {
    }

    /**
     * @brief Process a single value
     *
     * @param value
     * @return T
     */
    auto
    operator()(const T& value) const -> T
        requires StatLab::Concepts::Scalar< T >
    {

        if constexpr (nan_policy == StatLab::Utils::NaNPolicy::propagate)
        {
            _init = _func(_init, value);
        }
        else if constexpr (nan_policy == StatLab::Utils::NaNPolicy::skip)
        {
            if (not std::isnan(value) and not std::isinf(value))
            {
                _init = _func(_init, value);
            }
        }
        else
        {
            if (not std::isnan(value) and not std::isinf(value))
            {
                _init = _func(_init, value);
            }
            else
            {
                throw std::invalid_argument("Error: non finite value found");
            }
        }

        return _init;
    }

    /**
     * @brief Process a range of values
     *
     * @tparam Iterable
     * @param values
     * @return T requires std::ranges::range< T >
     */
    template < std::ranges::range I >
    auto
    operator()(I&& values) const -> T
        requires std::ranges::range< T >
    {
        return RangeAlgo{}(
            std::forward< I >(values),
            [this](auto&& init, auto&& value)
            { return this->operator()(init, value); },
            _init);
    }

    /**
     * @brief Process a range fo values extracted from another range
     *
     * @tparam I
     * @tparam Getter
     * @param values
     * @param getter
     * @return T
     */
    template < std::ranges::range                       I,
               std::invocable< typename I::value_type > Getter >
    auto
    operator()(I&& values, Getter&& getter) -> T
    {
        return operator()(values | std::ranges::views::transform(getter));
    }
};

/**
 * @brief TODO
 *
 */
struct naive_iteration
{
    /**
     * @brief TODO
     *
     * @tparam T
     * @tparam V
     * @tparam Func
     * @param range
     * @param func
     * @param init
     * @return decltype(auto)
     */
    template < std::ranges::range T, typename V, std::invocable< T > Func >
    auto
    operator()(T&& range, Func&& func, V&& init) -> decltype(auto)
    {
        for (auto v : range)
        {
            init = func(init, v);
        }

        return init;
    }
};

/**
 * @brief TODO
 *
 */
struct pairwise_iteration
{
    /**
     * @brief Pairwise, i.e., 
     *
     * @tparam T
     * @tparam V
     * @tparam Func
     * @tparam splitsize
     * @param range
     * @param func
     * @param init
     * @return decltype(auto)
     */
    template < std::ranges::range T,
               typename V,
               std::invocable< T, T > Func,
               int                    splitsize = 2500 >
    auto
    operator()(T&& range, Func&& func, V&& init) -> decltype(auto)
    {
        auto size     = std::ranges::size(range);
        auto sum_algo = naive_iteration{};
        if (size <= splitsize)
        {
            return sum_algo(range, func);
        }
        else
        {
            return this->operator()(std::ranges::subrange(
                                        std::ranges::begin(range),
                                        std::ranges::begin(range) + splitsize),
                                    func,
                                    init) +
                   this->operator()(std::ranges::subrange(
                                        std::ranges::begin(range) + splitsize,
                                        std::ranges::end(range)),
                                    func,
                                    init);
        }
    }
};

/**
 * @brief Kahan's algorithm
 *
 * @tparam T
 */
template < Concepts::Scalar T >
struct kahan_sum_single
{
    T _y    = T{};
    T _t    = T{};
    T _comp = T{};

    /**
     * @brief TODO
     *
     * @param init
     * @param value
     * @return T
     */
    auto
    operator()(T init, T value) const -> T
    {
        _y    = value - _comp;
        _t    = init + _y;
        _comp = (_t - init) - _y;
        init  = _t;
        return init;
    }
};

/**
 * @brief Sum values using naive iteration
 *
 * @tparam T
 * @tparam nan_policy
 */
template < typename T, StatLab::Utils::NaNPolicy nan_policy >
using Sum = ArithmeticOperation< T, nan_policy, std::plus<T>, naive_iteration >;

/**
 * @brief Sum values using Pairwise summation
 *
 * @tparam T
 * @tparam nan_policy
 */
template < typename T, StatLab::Utils::NaNPolicy nan_policy >
using SumPaiwise =
    ArithmeticOperation< T, nan_policy, std::plus<T>, pairwise_iteration >;

/**
 * @brief Sum values using Kahan's algorithm
 *
 * @tparam T
 * @tparam nan_policy
 */
template < typename T, StatLab::Utils::NaNPolicy nan_policy >
using SumKahan =
    ArithmeticOperation< T, nan_policy, kahan_sum_single<T>, naive_iteration >;

// multiplication shit goes here

/**
 * @brief Multiply values
 * 
 * @tparam T 
 * @tparam nan_policy 
 */
template < typename T, StatLab::Utils::NaNPolicy nan_policy >
using Multiply =
    ArithmeticOperation< T, nan_policy, std::multiplies<T>, naive_iteration >;

}
}
#endif
