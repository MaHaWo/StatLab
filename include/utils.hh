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
 * @brief Helper function for removing pointer qualifiers from a type recursivly
 *        - recursion base case which provides a type equal to T
 * @tparam T
 * @tparam 0
 */
template <typename T, typename U = std::void_t<>>
struct remove_pointer
{
    using type = T;
};

/**
 * @brief Helper function for removing pointer qualifiers from a type recursivly
 *        Provides a member type definition called 'type' which is equal to T
 *        if the first template argument is of type T* or T** or T***...
 * @tparam T
 * @tparam 0
 */
template <typename T>
struct remove_pointer<T, std::enable_if_t<std::is_pointer_v<T>, std::void_t<>>>
{
    using type = typename remove_pointer<std::remove_pointer_t<T>>::type;
};

/**
 * @brief Shorthand for 'typename remove_pointer<T>::type'
 *
 * @tparam T
 */
template <typename T>
using remove_pointer_t = typename remove_pointer<T>::type;

/**
 * @brief Oveload of 'remove_pointer' metafunction for array types (stack
 * allocated)
 *
 * @tparam T
 */
template <typename T>
struct remove_pointer<T, std::enable_if_t<std::is_array_v<T>, std::void_t<>>>
{
    using type = typename remove_pointer<std::remove_all_extents_t<T>>::type;
};

// remove qualifiers. FIXME: this is not optimal currently, because it does not
// work recursivly

/**
 * @brief Function for removing the qualifiers from
 *
 * @tparam T
 */
template <typename T>
struct remove_qualifier
{
    using type = std::remove_cv_t<remove_pointer_t<std::remove_reference_t<T>>>;
};

/**
 * @brief Shorthand for 'typename remove_qualifier::value'
 *
 * @tparam T
 */
template <typename T>
using remove_qualifier_t = typename remove_qualifier<T>::type;

/**
 * @brief Meta-function for checking if a type is a stringtype, that is
 *        it is either std::string, const char* or char* - prototype
 *
 * @tparam T
 */
template <typename T>
struct is_string_helper : public std::false_type
{
};

/**
 * @brief Specialization of 'is_string_helper' for std::string
 *
 * @tparam
 */
template <>
struct is_string_helper<std::string> : public std::true_type
{
};

/**
 * @brief Metafunction for determining if a type is a string-like type, i.e.
 *        std::string, const char*, char*
 * @tparam T
 */
template <typename T>
struct is_string : public is_string_helper<remove_qualifier_t<T>>
{
};

/**
 * @brief Overload of is_string for pure const char*, which would loose its
 *        pointer qualifier if this was not provided
 * @tparam
 */
template <>
struct is_string<const char*> : public std::true_type // is_string_helper<const char*>
{
};

/**
 * @brief Overload of is_string for pure const char*, which would loose its
 *        pointer qualifier if this was not provided
 * @tparam
 */
template <>
struct is_string<char*> : public std::true_type // public is_string_helper<char*>
{
};

/**
 * @brief Shorthand for 'is_string<T>::value'
 *

 */
template <typename T>
constexpr inline bool is_string_v = is_string<T>::value;

/**
 * @brief Metafunction for checking if a type is a containertype, which does
 *        not include string types - prototype
 *
 * @tparam T
 * @tparam std::void_t<>
 */
template <class T, class U = std::void_t<>>
struct is_container : public std::false_type
{
};

/**
 * @brief Metafunction for checking if a type is a containertype, which does
 *        not include string types
 *
 * @tparam T
 */
template <typename T>
struct is_container<T, std::void_t<typename remove_qualifier_t<T>::iterator, std::enable_if_t<!is_string_v<T>, int>>>
    : public std::true_type
{
};

/**
 * @brief Shorthand for 'is_container::value
 *
 * @tparam T
 */
template <typename T>
inline constexpr bool is_container_v = is_container<T>::value;

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
    template <typename T>
    constexpr bool operator()(const T& lhs, const T& rhs, [[maybe_unused]] double tol = 2e-15)
    {
        if constexpr (std::is_floating_point_v<T>)
        {
            if (std::abs(lhs) < tol && std::abs(rhs) < tol)
            {
                return true;
            }
            else
            {
                // // from 'The Art of computer programming, vol 2, p. 218'
                // return std::abs(lhs - rhs) <=
                //        ((std::abs(lhs) > std::abs(rhs) ? std::abs(rhs) : std::abs(lhs)) * tol);
                return std::abs(lhs - rhs) < tol;
            }
        }
        else if constexpr (std::is_integral_v<T>)
        {
            return lhs == rhs;
        }
        else if constexpr (is_container_v<T>)
        {
            bool equal = lhs.size() == rhs.size();
            if (equal)
            {
                for (std::size_t i = 0; i < lhs.size(); ++i)
                {
                    equal = IsEqual()(lhs[i], rhs[i], tol);
                    if (!equal)
                    {
                        break;
                    }
                }
            }
            return equal;
        }
        else if constexpr (std::is_pointer_v<T>)
        {
            return lhs == rhs;
        }
        else
        {
            throw std::runtime_error("Unknown objects to compare in is_equal");
            return false;
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
    template <typename T>
    constexpr bool operator()(const T& lhs, const T& rhs)
    {
        bool less = true;
        if constexpr (is_container_v<T>)
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
    template <typename T>
    constexpr bool operator()(const T& lhs, const T& rhs)
    {
        if constexpr (is_container_v<T>)
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
template <typename T, std::enable_if_t<is_container_v<T>, int> = 0>
std::ostream& operator<<(std::ostream& out, const T& container)
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

// namespace for implementing tuple for each
namespace visit_impl
{
// functions needed for implementing the visit from below
template <std::size_t idx, typename Function, typename... Tuplelike>
constexpr auto apply_to_index(Function&& f, Tuplelike&&... tuplelike)
{
    using std::get;
    return std::apply(f, std::forward_as_tuple(std::get<idx>(tuplelike)...));
}

// function for applying the above to each tuple- or array- or pair-element
template <class Returntype, typename Function, std::size_t... idxs, typename... Tuplelike>
constexpr auto apply_to_indices(Function&& f, std::index_sequence<idxs...>, Tuplelike&&... tuples)
{
    using std::get;

    if constexpr (std::is_same_v<Returntype, void>)
    {
        (apply_to_index<idxs>(std::forward<Function>(f), std::forward<Tuplelike>(tuples)...), ...);
    }
    else
    {
        return Returntype{apply_to_index<idxs>(
            std::forward<Function>(f), std::forward<Tuplelike>(tuples)...)...};
    }
}

} // namespace visit_impl

/**
 * @brief      Function which takes an arbitrary number of tuples, and reduces
 * them along equal indices, via the function f, which takes each nth element of
 * the tuples as arguments and returns a single value:
 * Example: std::tuple<int, double, int> a;
 *          std::tuple<double, double, int> b;
 *          std::tuple<double, double, double> res;
 *          res =  visit([](auto& i, auto& j){ return i*j;}, a, b);
 *
 * @param[in]  f function to reduce the tuples
 * @param[in]  tuples  Arbitrary many tuples or related objects to reduce
 *
 * @return     Returntype holding the results of the reduction operation
 */
template <class Returntype, typename Function, typename... Tuplelike>
constexpr auto visit(Function&& f, Tuplelike&&... tuples)
{
    using std::get;
    using std::tuple_size;

    constexpr std::size_t N =
        std::min({tuple_size<typename remove_qualifier<Tuplelike>::type>::value...});

    return visit_impl::apply_to_indices<Returntype>(
        std::forward<Function>(f), std::make_index_sequence<N>(),
        std::forward<Tuplelike>(tuples)...);
}

/**
 * @brief Functor for extracting basic type from collection.
 *
 * @tparam T
 * @tparam void
 */
template <typename T, typename U = void>
struct extract_basic_type
{
};

template <typename T>
struct extract_basic_type<T, std::enable_if_t<!is_container_v<T>>>
{
    using type = T;
};

/**
 * @brief
 *
 * @tparam T
 */
template <typename T>
struct extract_basic_type<T, std::enable_if_t<is_container_v<T>>>
{
    using type = typename T::value_type;
};

} // namespace Utils
} // namespace Statistics

#endif
