#ifndef STATISTICS_UTILS_HH
#define STATISTICS__UTILS_HH

#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <type_traits>

namespace Statistics{
namespace Utils
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
                return std::abs(lhs - rhs) / std::max(std::abs(lhs), std::abs(rhs)) < tol;
            }
        }
        else if constexpr (std::is_integral_v<T>)
        {
            return lhs == rhs;
        }
        else if constexpr (DataIO::is_container_v<T>)
        {
            bool equal = lhs.size() == rhs.size();
            if (equal)
            {
                for (std::size_t i = 0; i < lhs.size(); ++i)
                {
                  equal = IsEqual()(lhs[i], rhs[i], tol);
                  if (!equal) {
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
            less = lhs.size() == rhs.size();
            if (less)
            {
                for (std::size_t i = 0; i < lhs.size(); ++i)
                {
                    less = (lhs[i] < rhs[i]);
                    if (!less)
                    {
                        break;
                    }
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
        bool greater = true;
        if constexpr (DataIO::is_container_v<T>)
        {
            greater = lhs.size() == rhs.size();
            if (equal)
            {
                for (std::size_t i = 0; i < lhs.size(); ++i)
                {
                    greater = (lhs[i] > rhs[i]);
                    if (!greater)
                    {
                        break;
                    }
                }
            }
        }
        else
        {
            greater = (lhs > rhs);
        }
        return greater;
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
template <typename T, std::enable_if_t<DataIO::is_container_v<T>, int> = 0>
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
namespace tuple_reduce_impl
{
// functions needed for implementing the tuple_reduce from below
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
        (apply_to_index<idxs>(std::forward<Function>(f),
                                     std::forward<Tuplelike>(tuples)...), ...);
    }
    else
    {
        return Returntype{apply_to_index<idxs>(
            std::forward<Function>(f), std::forward<Tuplelike>(tuples)...)...};
    }
}

} // namespace tuple_reduce_impl

/**
 * @brief      Function which takes an arbitrary number of tuples, and reduces
 * them along equal indices, via the function f, which takes each nth element of
 * the tuples as arguments and returns a single value:
 * Example: std::tuple<int, double, int> a;
 *          std::tuple<double, double, int> b;
 *          std::tuple<double, double, double> res;
 *          res =  tuple_reduce([](auto& i, auto& j){ return i*j;}, a, b);
 *
 * @param[in]  f function to reduce the tuples
 * @param[in]  tuples  Arbitrary many tuples or related objects to reduce
 *
 * @return     Returntype holding the results of the reduction operation
 */
template <class Returntype, typename Function, typename... Tuplelike>
constexpr auto tuple_reduce(Function&& f, Tuplelike&&... tuples)
{
    using std::get;
    using std::tuple_size;

    constexpr std::size_t N = std::min(
        {tuple_size<typename DataIO::remove_qualifier<Tuplelike>::type>::value...});

    return tuple_reduce_impl::apply_to_indices<Returntype>(
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
struct ExtractBasicType
{
};

template <typename T>
struct ExtractBasicType<T, std::enable_if_t<!is_container_v<T>>>
{
    using type = T;
};

template <typename T>
struct ExtractBasicType<T, std::enable_if_t<is_container_v<T>>>
{
    using type = typename T::value_type;
};

} // namespace Utils
} // namespace Statistics

#endif
