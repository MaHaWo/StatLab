#ifndef STATLAB_CONCEPTS_HPP
#define STATLAB_CONCEPTS_HPP

#include "llvm/15.0.6/include/c++/v1/__concepts/convertible_to.h"
#include "llvm/15.0.6/include/c++/v1/__concepts/invocable.h"
#include <concepts>
#include <iterator>
#include <limits>
#include <type_traits>
#include <ranges> 

namespace StatLab
{
namespace Concepts
{

template < typename T > 
concept Scalar = std::is_scalar_v<T>;


template <typename T, typename Arg1, typename Arg2> 
concept BinaryCallable = requires(T object, Arg1 arg1, Arg2 arg2)
{
    std::invocable<T, Arg1, Arg2>;
    std::convertible_to<Arg1, Arg2>;
    std::convertible_to<Arg2, Arg1>;
    {object(arg1, arg2)} -> std::convertible_to<Arg1>;
};

/**
 * @brief Concept that names types that can be added together, i.e., have
 * operators '+' and '-'
 *
 * @tparam T
 */
template < typename T >
concept Additive = requires(T lhs, T rhs) {
                       {
                           lhs + rhs
                           } -> std::convertible_to< T >;
                       {
                           lhs - rhs
                           } -> std::convertible_to< T >;
                   };

/**
 * @brief Concept that names types that can be multiplied together, i.e., have
 * operators '*' and '/'
 *
 * @tparam T
 */
template < typename T >
concept Multiplicative = requires(T lhs, T rhs) {
                             {
                                 lhs* rhs
                                 } -> std::convertible_to< T >;
                             {
                                 lhs / rhs
                                 } -> std::convertible_to< T >;
                         };

}
}
#endif