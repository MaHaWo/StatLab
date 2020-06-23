#ifndef STATLAB_CORE_HPP
#define STATLAB_CORE_HPP

#include "utilities.hpp"

namespace StatLab
{

/**
 * @brief TODO
 *
 * @tparam T
 */
template < typename T >
class Product final
{
    T _result;

  public:
    /**
     * @brief TODO
     *
     * @param value
     * @return T
     */
    inline auto
    operator()(T value) -> T
    {
        _result *= T;
        return _result;
    }

    /**
     * @brief TODO
     *
     * @tparam Iter
     * @tparam Getter
     * @param begin
     * @param end
     * @param getter
     * @return auto
     */
    template < typename Iter, typename Getter >
    inline auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter)
    {
        std::for_each(std::forward< Iter >(begin),
                      std::forward< Iter >(end),
                      [&getter, &_result](auto&& v) { _result *= v; });

        return _result;
    }

    /**
     * @brief Construct a new Product object
     *
     */
    Product() : _result(T{})
    {
    }

    /**
     * @brief Construct a new Product object
     *
     */
    Product(const Product&) = default;

    /**
     * @brief Construct a new Product object
     *
     */
    Product(Product&&) = default;

    /**
     * @brief Assign Product object
     *
     * @return Product&
     */
    Product&
    operator=(const Product&) = default;

    /**
     * @brief Assign Product object
     *
     * @return Product&
     */
    Product&
    operator=(Product&&) = default;

    /**
     * @brief Destroy the Product object
     *
     */
    ~Product();
};

/**
 * @brief Sum up a range of (real) numbers using the kahan summation algorithm
 * to limit the accumulation of summation error:
 * https://en.wikipedia.org/wiki/Kahan_summation_algorithm
 *
 * @tparam T
 */
template < typename T >
class SumKahan final
{
    T _result;
    T _comp;
    T _y;
    T _t;

  public:
    // rule of five + custom constructor
    SumKahan() : _result(T{}), _comp(T{}), _y(T{}), _t(T{})
    {
    }
    SumKahan(const SumKahan&) = default;
    SumKahan(SumKahan&&)      = default;
    SumKahan&
    operator=(const SumKahan&) = default;
    SumKahan&
    operator=(SumKahan&&) = default;
    ~SumKahan()           = default;

    /**
     * @brief TODO
     *
     * @return T
     */
    inline constexpr T
    result()
    {
        return _result;
    }

    /**
     * @brief Reset internal variable to zero or equivalent for type T
     *
     */
    inline constexpr void
    reset()
    {
        _result = T{};
    }

    /**
     * @brief TODO
     *
     * @tparam BinaryFunc
     * @param value
     * @param plus
     * @param Minus
     * @return T
     */
    template < typename BinaryFunc >
    inline constexpr auto
    operator()(T value, BinaryFunc&& plus) -> T
    {
        _y      = plus(value, -1 * _comp);
        _t      = plus(_result, _y);
        _comp   = plus(plus(_t, -1 * _result), -1 * _y);
        _result = _t;
        return _result;
    }

    /**
     * @brief
     *
     * @param value
     * @return T
     */
    inline constexpr auto
    operator()(T value) -> T
    {
        return operator()(value, std::plus<>{});
    }

    /**
     * @brief TODO
     *
     * @tparam Iter
     * @tparam Getter
     * @param begin
     * @param end
     * @param getter
     * @return Utils::requires<Utils::is_addable, T>
     */
    template < typename Iter, typename Getter, typename BinaryFunc >
    inline auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter, BinaryFunc&& plus)
        -> T
    {
        for (; begin != end; ++begin)
        {
            this->operator()(getter(*begin), std::forward< BinaryFunc >(plus));
        }
        return _result;
    }

    /**
     * @brief TODO
     *
     * @tparam Iter
     * @tparam Getter
     * @param begin
     * @param end
     * @param getter
     * @return auto
     */
    template < typename Iter, typename Getter >
    inline auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter) -> T
    {
        return operator()(std::forward< Iter >(begin),
                          std::forward< Iter >(end),
                          std::forward< Getter >(getter),
                          std::plus<>{});
    }

    /**
     * @brief
     *
     * @tparam Iter
     * @param begin
     * @param end
     * @return auto
     */
    template < typename Iter >
    inline auto
    operator()(Iter&& begin, Iter&& end) -> T
    {
        return operator()(std::forward< Iter >(begin),
                          std::forward< Iter >(end),
                          Utils::Identity< T >{});
    }
};

/**
 * @brief Sum a range of numbers using pairwise summation algorithm
 *
 * @tparam T
 */
template < typename T >
class SumPairwise final
{
    T _result;

  public:
    // rule of file
    SumPairwise() : _result(T{})
    {
    }
    SumPairwise(const SumPairwise&) = default;
    SumPairwise(SumPairwise&&)      = default;
    SumPairwise&
    operator=(const SumPairwise&) = default;
    SumPairwise&
    operator=(SumPairwise&&) = default;
    ~SumPairwise()           = default;

    /**
     * @brief TODO
     *
     * @return T
     */
    inline constexpr T
    result()
    {
        return _result;
    }

    /**
     * @brief Reset the internal variables to zero or equivalent for type T
     *
     */
    inline constexpr void
    reset()
    {
        _result = T{};
    }

    /**
     * @brief Add a number 'value' of type T to the value stored in the caller
     *
     * @param value value to add
     * @return constexpr T
     */
    template < typename BinaryFunc >
    inline constexpr auto
    operator()(T value, BinaryFunc&& plus) -> T
    {
        _result = plus(_result, value);
        return _result;
    }

    /**
     * @brief TODO
     *
     * @param value
     * @return T
     */
    inline constexpr auto
    operator()(T value) -> T
    {
        return operator()(value, std::plus<>{});
    }

    /**
     * @brief Add a range [begin, end) of values to the value stored by the
     * caller. The range is represented by Iterators begin and end, and the
     * values to add are the result of applying the callable 'getter' to the
     * result of dereferencing each iterator.
     * Uses partition sum to limit numerical errors.
     * @tparam Iter Iterator type
     * @tparam Getter Callable type to use
     * @tparam base base size of range which tells the partition summation at
     * which range size to evaluate the sum. If the range is longer than base,
     *         recursion will be used with subranges until these are smaller or
     * equal to base.
     * @param begin Iterator to begin of range
     * @param end Iterator to end of range
     * @param getter Callable which is invoked on each result of dereferncing an
     * iterator
     * @param plus Function for adding two objects of type T
     * @return T sum over the range, with each element being passed through
     * 'getter' first.
     */
    template < typename Iter,
               typename Getter,
               typename BinaryFunc,
               std::size_t base = 1000 >
    inline auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter, BinaryFunc&& plus)
        -> T
    {

        using IVT =
            typename std::iterator_traits< std::decay_t< Iter > >::value_type;

        // check that the getter is usefully invocable
        static_assert(
            std::is_invocable_v< std::decay_t< Getter >, IVT >,
            "Error, Getter cannot be invoked with to value type of iterator");

        double result = 0;
        double size   = std::distance(begin, end);
        // use partition sum to limit numerical errors
        if (size <= base)
        {
            // recursion base case
            for (; begin != end; ++begin)
            {
                _result = plus(_result, getter(*begin));
            }
        }
        else
        {

            std::size_t half = std::round(size / 2.);

            // recursion
            this->operator()(std::next(std::forward< Iter >(begin), half),
                             std::forward< Iter >(end),
                             std::forward< Getter >(getter),
                             std::forward< BinaryFunc >(plus));

            this->operator()(std::forward< Iter >(begin),
                             std::next(std::forward< Iter >(begin), half),
                             std::forward< Getter >(getter),
                             std::forward< BinaryFunc >(plus));
        }

        return _result;
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
        return operator()(std::forward< Iter >(begin),
                          std::forward< Iter >(end),
                          std::forward< Getter >(getter),
                          std::plus<>{});
    }

    /**
     * @brief TODO
     *
     * @tparam Iter
     * @param begin
     * @param end
     * @return auto
     */
    template < typename Iter >
    inline auto
    operator()(Iter&& begin, Iter&& end)
    {
        return operator()(std::forward< Iter >(begin),
                          std::forward< Iter >(end),
                          Utils::Identity< T >{});
    }
};

//==============================================================================
//==============================================================================
//============================ Sum based statistics ============================
//==============================================================================
//==============================================================================

/**
 * @brief TODO
 *
 * @tparam T
 * @tparam SumAlgorithm
 * @tparam order
 */
template < typename T,
           template < typename >
           class SumAlgorithm,
           std::size_t order >
class Moment
{
    SumAlgorithm< T > _accumulator;
    T                 _size;

  public:
    Moment() : _accumulator(SumAlgorithm< T >{}), _size(T{})
    {
    }
    Moment(const Moment&) = default;
    Moment(Moment&&)      = default;
    Moment&
    operator=(const Moment&) = default;
    Moment&
    operator=(Moment&&) = default;
    ~Moment()           = default;

    /**
     * @brief Return the result stored in the functor.
     *
     * @return T
     */
    inline constexpr T
    result()
    {
        return _accumulator.result() / _size;
    }

    /**
     * @brief reset the functor, setting its internally stored sum value to zero
     *
     */
    inline constexpr void
    reset()
    {
        _accumulator.reset();
        _size = 0;
    }

    /**
     * @brief Add a single value to the sum stored within this functor
     *
     * @param value Value to add
     * @return current summed value stored within this functor
     */
    inline constexpr auto
    operator()(T value) -> T
    {
        _size += 1;
        _accumulator(std::pow(value, order));
        return _accumulator.result() / _size;
    }

    /**
     * @brief TODO
     *
     * @tparam Iter
     * @tparam Getter
     * @tparam BinaryFunc1
     * @tparam BinaryFunc2
     * @param begin
     * @param end
     * @param getter
     * @param plus
     * @return T
     */
    template < typename Iter, typename Getter, typename BinaryFunc >
    inline auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter, BinaryFunc&& plus)
        -> T
    {
        _size += std::distance(std::forward< Iter >(begin),
                               std::forward< Iter >(end));

        return _accumulator(
                   std::forward< Iter >(begin),
                   std::forward< Iter >(end),
                   [getter](auto&& v) { return std::pow(getter(v), order); },
                   std::forward< BinaryFunc >(plus)) /
               _size;
    }

    /**
     * @brief TODO
     *
     * @tparam Iter
     * @tparam Getter
     * @param begin
     * @param end
     * @param getter
     * @return auto
     */
    template < typename Iter, typename Getter >
    inline auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter) -> T
    {
        return operator()(std::forward< Iter >(begin),
                          std::forward< Iter >(end),
                          std::forward< Getter >(getter),
                          std::plus<>{});
    }

    /**
     * @brief
     *
     * @tparam Iter
     * @param begin
     * @param end
     * @return auto
     */
    template < typename Iter >
    inline auto
    operator()(Iter&& begin, Iter&& end) -> T
    {
        return operator()(std::forward< Iter >(begin),
                          std::forward< Iter >(end),
                          Utils::Identity< T >{});
    }
};

/**
 * @brief TODO
 *
 * @tparam T
 * @tparam SumAlgorithm
 */
template < typename T, template < typename > class SumAlgorithm >
class ArithmeticMean final : public Moment< T, SumAlgorithm, 1 >
{
};

/**
 * @brief TODO
 *
 * @tparam T
 * @tparam SumAlgorithm
 */
template < typename T, template < typename > class SumAlgorithm >
class HarmonicMean final
{
    Moment< T, SumAlgorithm, 1 > _moment;

  public:
    /**
     * @brief TODO
     *
     * @tparam BinaryFunc
     * @param value
     * @param plus
     * @return constexpr auto
     */
    template < typename BinaryFunc >
    inline constexpr auto
    operator()(T value, BinaryFunc&& plus)
    {
        _moment(static_cast< T >(1) / value, plus);
        return static_cast< T >(1) / _moment.result();
    }

    /**
     * @brief TODO
     *
     * @tparam Iter
     * @tparam Getter
     * @tparam BinaryFunc
     * @param begin
     * @param end
     * @param getter
     * @param plus
     * @return constexpr auto
     */
    template < typename Iter, typename Getter, typename BinaryFunc >
    inline constexpr auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter, BinaryFunc&& plus)
        -> T
    {
        return static_cast< T >(1) / _moment(
                                         std::forward< Iter >(begin),
                                         std::forward< Iter >(end),
                                         [&getter](auto&& v) {
                                             return static_cast< T >(1) /
                                                    getter(v);
                                         },
                                         plus);
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
    inline constexpr auto
    operator()(Iter&& begin, Iter&& end, Getter&& getter) -> T
    {
        return operator()(std::forward< Iter >(begin),
                          std::forward< Iter >(end),
                          std::forward< Getter >(getter),
                          std::plus<>{});
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
    inline constexpr auto
    operator()(Iter&& begin, Iter&& end) -> T
    {
        return operator()(std::forward< Iter >(begin),
                          std::forward< Iter >(end),
                          Utils::Identity< T >{});
    }

    HarmonicMean() : _moment(Moment< T, SumAlgorithm, 1 >{})
    {
    }
    HarmonicMean(const HarmonicMean&) = default;
    HarmonicMean(HarmonicMean&&)      = default;
    HarmonicMean&
    operator=(const HarmonicMean&) = default;
    HarmonicMean&
    operator=(HarmonicMean&&) = default;
    ~HarmonicMean()           = default;
};

/**
 * @brief
 *
 * @tparam T
 * @tparam SumAlgorithm
 */
template < typename T >
class GeometricMean final
{

  public:
    GeometricMean()                     = default;
    GeometricMean(const GeometricMean&) = default;
    GeometricMean(GeometricMean&&)      = default;
    GeometricMean&
    operator=(const GeometricMean&) = default;
    GeometricMean&
    operator=(GeometricMean&&) = default;
    ~GeometricMean()           = default;
};

} // namespace StatLab

#endif