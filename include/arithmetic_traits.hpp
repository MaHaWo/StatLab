#ifndef STATLAB_ARITHMETIC_TRAITS_HPP
#define STATLAB_ARITHMETIC_TRAITS_HPP
#include <cmath>
#include <type_traits>
#include <utility>
#include <functional>

namespace StatLab {
namespace Traits {

/**
 * @brief Determine if a type T has a multiplicative inverse
 *        Prototype version
 *
 * @tparam T Type to check for multiplicative inverse existence
 */
template<typename T, typename V = std::void_t<>>
struct has_multiplicative_inverse : std::false_type
{
};

/**
 * @brief Determine if a type T has a multiplicative inverse
 *
 * @tparam T Type to check for multiplicative inverse existence
 */
template<typename T>
struct has_multiplicative_inverse<T, std::void_t<decltype(std::pow(std::declval<T>(), -1))>>
  : std::true_type
{
};

/**
 * @brief Shorthand for has_multiplicative_inverse<T>::value
 *
 * @tparam T Type to check if it has multiplicative inverse
 */
template<typename T>
inline constexpr bool has_multiplicative_inverse_v =
  has_multiplicative_inverse<T>::value;

/**
 * @brief Determine if a type T has an additive inverse
 *        Prototype version
 *
 * @tparam T Type to check if it has additive inverse
 */
template<typename T, typename V = std::void_t<>>
struct has_additive_inverse : std::false_type
{
};

/**
 * @brief Determine if a type T has an additive inverse
 *
 * @tparam T
 */
template<typename T>
struct has_additive_inverse<T, std::void_t<std::negate<T>>> : std::true_type
{
};

/**
 * @brief Shorthand for has_additive_inverse<T>::value
 *
 * @tparam T Type to check if it has additive inverse
 */
template<typename T>
inline constexpr bool has_additive_inverse_v = has_additive_inverse<T>::value;

/**
 * @brief Check if operation Op for type T is commutative, prototype
 *
 * @tparam T type to use Op on
 * @tparam Op Operation to use, e.g. std::plus
 * @tparam std::void_t<>
 */
template<typename T, template<typename> class Op, typename U = std::void_t<>>
struct is_commutative
{
    static constexpr bool value = Op<T>::is_commutative;
};

/**
 * @brief Check if an operation 'Op' on a type 'T' is commutative
 *
 * @tparam T type to use Op on
 * @tparam Op Operation to use, e.g. std::plus
 */
template<typename T, template<typename> class Op>
struct is_commutative<T, Op, std::void_t<decltype(Op<T>::is_commutative)>>
{
    static constexpr bool value = Op<T>::is_commutative;
};

/**
 * @brief Commutativity for std::plus with arithmetic types is true
 *
 * @tparam T Type to check for commutativity, only valid if it is arithmetic
 */
template<typename T>
struct is_commutative_plus
  : is_commutative<T, std::plus, std::enable_if_t<std::is_arithmetic_v<T>, int>>
{
    static constexpr bool value = true;
};

template<typename T> 
inline constexpr bool is_commutative_plus_v = is_commutative_plus<T>::value;

/**
 * @brief Commutativity for std::multiplies with arithmetic types is true
 *
 * @tparam T Type to check for commutativity, only valid if it is arithmetic
 */
template<typename T>
struct is_commutative_multiply
  : is_commutative<T,
                   std::multiplies,
                   std::enable_if_t<std::is_arithmetic_v<T>, int>>
{
    static constexpr bool value = true;
};

template<typename T> 
inline constexpr bool is_commutative_multiply_v = is_commutative_multiply<T>::value;

/**
 * @brief Shorthand for is_commutative<T>::value
 *
 * @tparam T type to check commutativity for
 */
template<typename T, template<typename> class Op>
inline constexpr bool is_commutative_v = is_commutative<T, Op>::value;

/**
 * @brief Identity function for type T, with value v of type T: f(v) -> v
 *
 * @tparam T type to define identity function for
 */
template<typename T>
struct Identity
{
    constexpr T operator()(T t) const noexcept
    {
        return t;
    }
};


/**
 * @brief TODO
 *
 * @tparam Condition
 * @tparam T
 */
template<template<typename...> class Condition, typename T>
using requires = std::enable_if_t<Condition<T>::value, T>;


/**
 * @brief TODO
 * 
 * @tparam Condition 
 * @tparam T 
 */
template<typename T, template<typename...> class... Conditions>
using requires_all_of = std::enable_if_t<std::conjunction_v<Conditions<T>...>, T>;



/**
 * @brief TODO
 * 
 * @tparam Condition 
 * @tparam T 
 */
/**
 * @brief TODO
 * 
 * @tparam Condition 
 * @tparam T 
 */
template<typename T, template<typename...> class... Conditions>
using requires_one_of = std::enable_if_t<std::disjunction_v<Conditions<T>...>, T>;


} // namespace Utils
} // namespace StatLab

#endif