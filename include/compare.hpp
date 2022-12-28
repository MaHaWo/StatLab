#ifndef STATLAB_COMPARE_HPP
#define STATLAB_COMPARE_HPP

#include <algorithm>
#include <concepts>
#include <ranges>

namespace StatLab
{
namespace Compare
{

/**
 * @brief General function for comparing things for lhs == rhs relation.
 *        Containers are compared elementwise, but only if they are equally long
 *        Unequal-length containers compare as fals
 *
 * @tparam T Automatically determined
 * @param lhs object to compare
 * @param rhs object to compare
 * @return true if (elementwise) lhs == rhs
 * @return false else
 */
struct IsEqual
{
    template < typename T >

    constexpr auto
    operator()(T&& lhs, T&& rhs, [[maybe_unused]] double tol = 2e-15) -> bool
    {
        if constexpr (std::is_floating_point_v< std::decay_t< T > >)
        {
            static_assert(std::equality_comparable< std::decay_t< T > >);

            if (std::abs(lhs) < tol && std::abs(rhs) < tol)
            {
                return true;
            }
            else
            {
                // // from 'The Art of computer programming, vol 2, p. 218'
                return std::abs(lhs - rhs) <=
                       ((std::abs(lhs) > std::abs(rhs) ? std::abs(rhs)
                                                       : std::abs(lhs)) *
                        tol);
            }
        }
        else if constexpr (std::ranges::range< std::decay_t< T > >)
        {
            static_assert(std::equality_comparable<
                          std::decay_t< decltype(*std::begin(lhs)) > >);

            bool equal = lhs.size() == rhs.size();
            if (equal)
            {
                for (std::size_t i = 0; i < lhs.size(); ++i)
                {
                    equal = operator()(*std::next(std::begin(lhs), i),
                                       *std::next(std::begin(rhs), i),
                                       tol);
                    if (!equal)
                    {
                        break;
                    }
                }
            }
            return equal;
        }
        else
        {
            static_assert(std::equality_comparable< std::decay_t< T > >);
            return lhs == rhs;
        }
    }
};

/**
 * @brief General function for comparing things for lhs > rhs relation.
 *        Containers are compared elementwise, but only if they are equally long
 *        Unequal-length containers compare as fals
 *
 * @tparam T Automatically determined
 * @param lhs object to compare
 * @param rhs object to compare
 * @return true if (elementwise) lhs > rhs
 * @return false else
 */
struct IsLess
{
    template < typename T >
        requires std::totally_ordered< T >
    constexpr auto
    operator()(T&& lhs, T&& rhs) -> bool
    {
        bool less = true;
        if constexpr (std::ranges::range< std::decay_t< T > >)
        {
            static_assert(std::totally_ordered<
                          std::decay_t< decltype(*std::begin(lhs)) > >);

            for (std::size_t i = 0; i < std::min(lhs.size(), rhs.size()); ++i)
            {
                less = *std::next(std::begin(lhs), i) <
                       *std::next(std::begin(rhs), i);
                if (!less)
                {
                    break;
                }
            }
        }
        else
        {
            static_assert(std::totally_ordered< std::decay_t< T > >);
            less = lhs < rhs;
        }
        return less;
    }
};

/**
 * @brief General function for comparing things for lhs < rhs relation.
 *        Containers are compared elementwise, but only if they are equally long
 *        Unequal-length containers compare as fals
 *
 * @tparam T Automatically determined
 * @param lhs object to compare
 * @param rhs object to compare
 * @return true if (elementwise) lhs < rhs
 * @return false else
 */
struct IsGreater
{
    template < typename T >
    constexpr auto
    operator()(T&& lhs, T&& rhs) -> bool
    {
        if constexpr (std::ranges::range< std::decay_t< T > >)
        {
            static_assert(std::totally_ordered<
                          std::decay_t< decltype(*std::begin(lhs)) > >);

            bool greater = true;

            for (std::size_t i = 0; i < std::min(lhs.size(), rhs.size()); ++i)
            {
                greater = *std::next(std::begin(lhs), i) >
                          *std::next(std::begin(rhs), i);
                if (!greater)
                {
                    greater = false;
                    break;
                }
            }
            return greater;
        }
        else
        {
            static_assert(std::totally_ordered< std::decay_t< T > >);
            return lhs > rhs;
        }
    }
};
}
}
#endif