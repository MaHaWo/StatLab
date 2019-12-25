#ifndef STATISTICS_UTILS_HH
#define STATISTICS_UTILS_HH

#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <type_traits>

namespace Statistics
{
namespace Utils
{


/**
 * @brief Meta-function for checking if a type is a stringtype, that is
 *        it is either std::string, const char* or char* - prototype
 *
 * @tparam T
 */
template < typename T >
struct is_string_helper : public std::false_type
{
};

/**
 * @brief Specialization of 'is_string_helper' for std::string
 *
 * @tparam
 */
template <>
struct is_string_helper< std::string > : public std::true_type
{
};

/**
 * @brief Metafunction for determining if a type is a string-like type, i.e.
 *        std::string, const char*, char*
 * @tparam T
 */
template < typename T >
struct is_string : public is_string_helper< std::decay_t< T > >
{
};

/**
 * @brief Overload of is_string for pure const char*, which would loose its
 *        pointer qualifier if this was not provided
 * @tparam
 */
template <>
struct is_string< const char* > : public std::true_type
{
};

/**
 * @brief Overload of is_string for pure const char*, which would loose its
 *        pointer qualifier if this was not provided
 * @tparam
 */
template <>
struct is_string< char* > : public std::true_type
{
};

/**
 * @brief Shorthand for 'is_string<T>::value'
 *
 */
template < typename T >
constexpr inline bool is_string_v = is_string< T >::value;

/**
 * @brief Metafunction for checking if a type is a containertype, which does
 *        not include string types - prototype
 *
 * @tparam T
 * @tparam std::void_t<>
 */
template < class T, class U = std::void_t<> >
struct is_container : public std::false_type
{
};

/**
 * @brief Metafunction for checking if a type is a containertype, which does
 *        not include string types
 *
 * @tparam T
 */
template < typename T >
struct is_container< T,
                     std::void_t< typename std::decay_t< T >::iterator,
                                  std::enable_if_t< !is_string_v< T >, int > > >
    : public std::true_type
{
};

/**
 * @brief Shorthand for 'is_container::value
 *
 * @tparam T
 */
template < typename T >
inline constexpr bool is_container_v = is_container< T >::value;

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
    constexpr bool
    operator()(const T& lhs, const T& rhs, [[maybe_unused]] double tol = 2e-15)
    {
        if constexpr (std::is_floating_point_v< T >)
        {
            if (std::abs(lhs) < tol && std::abs(rhs) < tol)
            {
                return true;
            }
            else
            {
                // // from 'The Art of computer programming, vol 2, p. 218'
                // return std::abs(lhs - rhs) <=
                //        ((std::abs(lhs) > std::abs(rhs) ? std::abs(rhs) :
                //        std::abs(lhs)) * tol);
                return std::abs(lhs - rhs) < tol;
            }
        }
        if constexpr (std::is_integral_v< T >)
        {
            return lhs == rhs;
        }
        if constexpr (is_container_v< T >)
        {
            bool equal = lhs.size() == rhs.size();
            if (equal)
            {
                for (std::size_t i = 0; i < lhs.size(); ++i)
                {
                    equal = operator()(lhs[i], rhs[i], tol);
                    if (!equal)
                    {
                        break;
                    }
                }
            }
            return equal;
        }
        if constexpr (std::is_pointer_v< T >)
        {
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
    constexpr bool
    operator()(const T& lhs, const T& rhs)
    {
        bool less = true;
        if constexpr (is_container_v< T >)
        {
            for (std::size_t i = 0; i < std::min(lhs.size(), rhs.size()); ++i)
            {
                less = (lhs[i] < rhs[i]);
                if (!less)
                {
                    break;
                }
            }
        }
        else
        {
            less = (lhs < rhs);
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
    constexpr bool
    operator()(const T& lhs, const T& rhs)
    {
        if constexpr (is_container_v< T >)
        {
            bool greater = true;

            for (std::size_t i = 0; i < std::min(lhs.size(), rhs.size()); ++i)
            {
                greater = (lhs[i] > rhs[i]);
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
            return (lhs > rhs);
        }
    }
};

/**
 * @brief Output operator for containers
 *
 * @tparam T
 * @tparam 0
 * @param out
 * @param container
 * @return std::ostream&
 */
template < typename T, std::enable_if_t< is_container_v< T >, int > = 0 >
std::ostream&
operator<<(std::ostream& out, const T& container)
{
    std::setprecision(16);
    out << "[";
    for (std::size_t i = 0; i < (container.size() - 1); ++i)
    {
        out << container[i] << ",";
    }
    out << container.back() << "]";

    return out;
}

template < typename T, typename U = void >
struct extract_basic_type
{
};

template < typename T >
struct extract_basic_type< T, std::enable_if_t< !is_container_v< T > > >
{
    using type = T;
};

template < typename T >
struct extract_basic_type< T, std::enable_if_t< is_container_v< T > > >
{
    using type = typename T::value_type;
};

} // namespace Utils
} // namespace Statistics

#endif
