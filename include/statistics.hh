/** @file statistics.hh
 *  @brief Implements functions for computing statistics.
 *
 *  @author Harald Mack, harald.mack@iup.uni-heidelberg.de
 *  @bug Quantiles return only exact quantiles, that is only
 *       values which are contained in the container.
 *  @todo Optimize, perhaps use template metaprogramming where appropriate
 *  @todo Implement interpolation for quantiles, since they are defined in the
 *  @todo Remove Base class
 * respective way!
 *
 */
#ifndef UTOPIA_MODELS_AMEEMULTI_STATISTICS_HH
#define UTOPIA_MODELS_AMEEMULTI_STATISTICS_HH

#include "utils.hh"
#include <algorithm>
#include <numeric>
#include <type_traits>
#include <unordered_map>



using namespace Statistics::Utils;
namespace Statistics
{
/**
 * @brief Functor for pairwise summation
 *
 * @tparam T return type
 */
template <typename T>
struct Sum
{
private:
    T _s;

public:
    using value_type = T;
    /**
     * @brief Construct a new Sum object
     *
     */
    Sum() : _s(T())
    {
    }
    Sum(T v) : _s(v)
    {
    }

    Sum(const Sum& other) = default;
    Sum(Sum&& other) = default;
    Sum& operator=(const Sum& other) = default;
    Sum& operator=(Sum&& other) = default;
    virtual ~Sum() = default;

    /**
     * @brief Reset the internal state
     *
     */
    void reset()
    {
        _s = T();
    }

    /**
     * @brief get result of computation
     *
     * @return T
     */
    T result()
    {
        return _s;
    }

    /**
     * @brief
     *
     *             Ignores NaN and Inf values
     *
     * @param[in]  begin          The begin iterator
     * @param[in]  end            The end iterator
     * @param[in]  f              function to apply to the values to
     * accumulate.
     *
     *
     * @return     sum over (f(it)) for 'it' in [begin, end)
     */

    template <typename InputIterator, typename Getter>
    inline value_type operator()(InputIterator&& begin, InputIterator&& end, Getter&& getter)
    {
        if (begin != end)
        {
            if (std::distance(begin, end) <= 1000)
            {
                value_type s = std::accumulate(
                    begin, end, value_type(),
                    [&](auto&& lhs, auto&& rhs) { return lhs + getter(rhs); });

                _s += s;
            }
            else
            {
                std::size_t t = std::floor(double(std::distance(begin, end) / 2));

                this->operator()(std::forward<InputIterator>(begin),
                                 std::next(std::forward<InputIterator>(begin), t),
                                 std::forward<Getter>(getter));

                this->operator()(std::next(std::forward<InputIterator>(begin), t),
                                 std::forward<InputIterator>(end),
                                 std::forward<Getter>(getter));
            }
        }
        return _s;
    }

    /**
     * @brief
     *
     * @tparam InputIterator
     * @param begin
     * @param end
     * @return value_type
     */
    template <typename InputIterator>
    inline value_type operator()(InputIterator&& begin, InputIterator&& end)
    {
        if (begin != end)
        {
            this->operator()(std::forward<InputIterator>(begin),
                             std::forward<InputIterator>(end),
                             [](auto& value) { return value; });
        }
        return _s;
    }

    /**
     * @brief updating the sum with a single value
     *
     * @param value
     * @return T
     */
    inline T operator()(T value)
    {
        _s += value;
        return _s;
    }
};

/**
 * @brief Functor for the kahan summation algorithm
 *
 * @tparam T
 */
template <typename T>
struct SumKahan
{
private:
    T _y;
    T _t;
    T _s;
    T _comp;

public:
    using value_type = T;
    /**
     * @brief Construct a new Sum Kahan object
     *
     */
    SumKahan() : _y(T()), _t(T()), _s(T()), _comp(T())
    {
    }

    /**
     * @brief Construct a new Sum Kahan object
     *
     * @param v initializer
     */
    SumKahan(T v) : _y(v), _t(v), _s(v), _comp(v)
    {
    }

    /**
     * @brief Construct a new Sum Kahan object
     *
     * @param other
     */
    SumKahan(const SumKahan& other) = default;

    /**
     * @brief Construct a new Sum Kahan object
     *
     * @param other
     */
    SumKahan(SumKahan&& other) = default;

    /**
     * @brief
     *
     * @param other
     * @return SumKahan&
     */
    SumKahan& operator=(const SumKahan& other) = default;

    /**
     * @brief
     *
     * @param other
     * @return SumKahan&
     */
    SumKahan& operator=(SumKahan&& other) = default;

    /**
     * @brief Destroy the Sum Kahan object
     *
     */
    virtual ~SumKahan() = default;
    /**
     * @brief reset internal state
     *
     */
    void reset()
    {
        _s = T();
        _y = T();
        _t = T();
        _comp = T();
    }

    /**
     * @brief get result of summation
     *
     * @return double
     */
    T result()
    {
        return _s;
    }
    /**
     * @brief      Summation of the range [begin, end) using Kahan's
     * algorithm.
     *
     * @param[in]  begin          The begin iterator
     * @param[in]  end            The end iterator
     * @param[in]  f              function to apply to the values to
     * accumulate.
     * @param[in]  args           additional arguments to f
     *
     *
     * @return     sum(f(it)) for it in [star, end)
     */

    template <typename InputIterator, typename Getter>
    inline T operator()(InputIterator&& begin, InputIterator&& end, Getter&& getter)
    {
        if (begin != end)
        {
            _comp = 0.;
            for (; begin != end; ++begin)
            {
                _y = getter(*begin) - _comp;
                _t = _s + _y;
                _comp = (_t - _s) - _y;
                _s = _t;
            }
        }
        return _s;
    }

    /**
     * @brief Summation of the range [begin, end) using Kahan's algorithm.
     *
     * @tparam InputIterator
     * @param begin begin iterator
     * @param end   end iterator
     * @return T
     */
    template <typename InputIterator>
    inline T operator()(InputIterator&& begin, InputIterator&& end)
    {
        return this->operator()(std::forward<InputIterator>(begin),
                                std::forward<InputIterator>(end),
                                [](auto&& value) { return value; });
    }

    /**
     * @brief Update the sum with a new value.
     *
     * @tparam Value Some type convertible to double
     * @param value
     * @return double
     */
    inline T operator()(T value)
    {
        _y = value - _comp;
        _t = _s + _y;
        _s = _t;
        return _s;
    }
};

template <typename T, int order, template <typename...> class Sum>
class Moment
{
protected:
    T _size;
    Sum<T> _sum;

public:
    inline T result()
    {
        return _sum.result() / _size;
    }

    inline void operator()(T value)
    {
        _size += 1;
        _sum(std::pow(value, order));
    }

    template <class InputIter, class Getter>
    T operator()(InputIter&& begin, InputIter&& end, Getter&& getter)
    {
        _size += std::distance(begin, end);
        _sum(std::forward<InputIter>(begin), std::forward<InputIter>(end),
             [&](auto&& value) { return std::pow(getter(value), order); });

        return result();
    }

    template <class InputIter>
    T operator()(InputIter&& begin, InputIter&& end)
    {
        return operator()(std::forward<InputIter>(begin), std::forward<InputIter>(end),
                          [](auto&& v) { return v; });
    }

    Moment() = default;
    Moment(const Moment& other) = default;
    Moment(Moment&& other) = default;
    Moment& operator=(Moment&& other) = default;
    Moment& operator=(const Moment& other) = default;
    virtual ~Moment() = default;
};

/**
 * Functor for computing the arithmetic mean. using a given summation algorithm.
 * @tparam: T type to represent the mean value.
 * @tparam: Sum Functor implementing a summation algorithm.
 */
template <typename T, template <typename...> class Sum>
using ArithmeticMean = Moment<T, 1, Sum>;

template <typename T, template <typename...> class Sum>
class HarmonicMean
{
protected:
    Sum<T> _mean;
    T _size;

public:
    void reset()
    {
        _mean.reset();
    }

    T result()
    {
        return _size / _mean.result();
    }

    T operator()(T value)
    {
        _size += 1;
        _mean(1. / value);
    }

    template <typename Iter, typename Getter>
    T operator()(Iter&& begin, Iter&& end, Getter&& getter)
    {
        _size += std::distance(std::forward<Iter>(begin), std::forward<Iter>(end));
        std::for_each(std::forward<Iter>(begin), std::forward<Iter>(end),
                      [&](auto&& value) { return 1. / getter(value); });

        return result();
    }

    template <typename Iter>
    T operator()(Iter&& begin, Iter&& end)
    {
        return operator()(std::forward<Iter>(begin), std::forward<Iter>(end),
                          [](auto&& v) { return v; });
    }

    HarmonicMean() = default;
    HarmonicMean(const HarmonicMean& other) = default;
    HarmonicMean(HarmonicMean&& other) = default;
    HarmonicMean& operator=(const HarmonicMean&& other) = default;
    HarmonicMean& operator=(HarmonicMean&& other) = default;
    virtual ~HarmonicMean() = default;
};

template <typename T, template <typename...> class Sum>
class Variance
{
protected:
    T _n;
    Sum<T> _mean;
    T _d;
    Sum<T> _M2;

public:
    inline T result()
    {
        return _M2.result() / (_n - 1.);
    }

    inline void reset()
    {
        _n = 0;
        _mean.reset();
        _d = 0;
        _M2.reset();
    }

    inline void operator()(T value)
    {
        _n += 1;

        _d = value - _mean.result();
        _mean(_d / _n);
        _M2(_d * (value - _mean.result()));
    }

    template <typename InputIter, typename Getter>
    T operator()(InputIter&& begin, InputIter&& end, Getter&& getter)
    {
        std::for_each(std::forward<InputIter>(begin), std::forward<InputIter>(end),
                      [&](auto&& value) { operator()(getter(value)); });
        return result();
    }

    template <typename InputIter>
    T operator()(InputIter&& begin, InputIter&& end)
    {
        return operator()(std::forward<InputIter>(begin), std::forward<InputIter>(end),
                          [](auto&& v) { return v; });
    }

    Variance() : _n(0), _d(0), _mean(Sum<T>()), _M2(Sum<T>())
    {
    }
    Variance(const Variance& other) = default;
    Variance(Variance&& other) = default;
    Variance& operator=(Variance&& other) = default;
    Variance& operator=(const Variance& other) = default;
    virtual ~Variance() = default;
};

template <typename T, template <typename...> class Sum>
class Standarddeviation
{
protected:
    Variance<T, Sum> _variance;

public:
    T result()
    {
        return std::sqrt(_variance.result());
    }

    void reset()
    {
        _variance.reset();
    }

    void operator()(T value)
    {
        _variance(value);
    }

    template <typename InputIter, typename Getter>
    T operator()(InputIter&& begin, InputIter&& end, Getter&& getter)
    {
        return std::sqrt(_variance(std::forward<InputIter>(begin),
                                   std::forward<InputIter>(end),
                                   std::forward<Getter>(getter)));
    }

    template <typename InputIter>
    T operator()(InputIter&& begin, InputIter&& end)
    {
        return operator()(std::forward<InputIter>(begin), std::forward<InputIter>(end),
                          [](auto&& v) { return v; });
    }

    Standarddeviation()
    {
    }
    Standarddeviation(const Standarddeviation& other) = default;
    Standarddeviation(Standarddeviation&& other) = default;
    Standarddeviation& operator=(Standarddeviation&& other) = default;
    Standarddeviation& operator=(const Standarddeviation& other) = default;
    virtual ~Standarddeviation() = default;
};

template <typename T, template <typename...> class Sum>
class Skewness
{
protected:
    T _n;
    Sum<T> _mean;
    Sum<T> _M2;
    Sum<T> _M3;
    T _n1;
    T _delta;
    T _delta_n;
    T _term1;

public:
    inline void reset()
    {
        _n = 0;
        _mean.reset();
        _M2.reset();
        _M3.reset();
        _n1 = 0;
        _delta = 0;
        _delta_n = 0;
        _term1 = 0;
    }

    inline T result()
    {
        return std::pow(_n, 0.5) * _M3.result() / std::pow(_M2.result(), 1.5);
    }

    void operator()(T value)
    {
        _n1 = _n;
        _n += 1;
        _delta = value - _mean.result();
        _delta_n = _delta / _n;
        _term1 = _delta * _delta_n * _n1;
        _mean(_delta_n);
        _M3(_term1 * _delta_n * (_n - 2) - 3 * _delta_n * _M2.result());
        _M2(_term1);
    }

    template <typename InputIter, typename Getter>
    T operator()(InputIter&& begin, InputIter&& end, Getter&& getter)
    {
        std::for_each(std::forward<InputIter>(begin), std::forward<InputIter>(end),
                      [&](auto&& value) { operator()(getter(value)); });
        return result();
    }

    template <typename InputIter>
    T operator()(InputIter&& begin, InputIter&& end)
    {
        return operator()(std::forward<InputIter>(begin), std::forward<InputIter>(end),
                          [](auto&& v) { return v; });
    }

    Skewness()
        : _n(0),
          _mean(Sum<T>()),
          _M2(Sum<T>()),
          _M3(Sum<T>()),
          _n1(0),
          _delta(0),
          _delta_n(0),
          _term1(0)
    {
    }
    Skewness(const Skewness& other) = default;
    Skewness(Skewness&& other) = default;
    Skewness& operator=(Skewness&& other) = default;
    Skewness& operator=(const Skewness& other) = default;
    virtual ~Skewness() = default;
};

template <typename T, template <typename...> class Sum>
class Kurtosis
{
protected:
    T _n;
    Sum<T> _mean;
    Sum<T> _M2;
    Sum<T> _M3;
    Sum<T> _M4;
    T _n1;
    T _delta;
    T _delta_n;
    T _delta_n2;
    T _term1;

public:
    /**
     * [result description]
     * @method result
     * @return [description]
     */
    inline T result()
    {
        return _n * _M4.result() / (std::pow(_M2.result(), 2)) - 3.;
    }

    inline void reset()
    {
        _n = 0;
        _mean = Sum<T>();
        _M2 = Sum<T>();
        _M3 = Sum<T>();
        _M4 = Sum<T>();
        _n1 = 0;
        _delta = 0;
        _delta_n = 0;
        _delta_n2 = 0;
        _term1 = 0;
    }

    /**
     * [operator description]
     * @method operator
     */
    void operator()(T value)
    {
        _n1 = _n;
        _n += 1;
        _delta = value - _mean.result();
        _delta_n = _delta / _n;
        _delta_n2 = _delta_n * _delta_n;
        _term1 = _delta * _delta_n * _n1;
        _mean(_delta_n);
        _M4(_term1 * _delta_n2 * (_n * _n - 3 * _n + 3) +
            6 * _delta_n2 * _M2.result() - 4 * _delta_n * _M3.result());
        _M3(_term1 * _delta_n * (_n - 2) - 3 * _delta_n * _M2.result());
        _M2(_term1);
    }

    /**
     * [operator description]
     * @method operator
     * @return [description]
     */
    template <typename InputIter, typename Getter>
    T operator()(InputIter&& begin, InputIter&& end, Getter&& getter)
    {
        std::for_each(std::forward<InputIter>(begin), std::forward<InputIter>(end),
                      [&](auto&& value) { operator()(getter(value)); });
        return result();
    }

    template <typename InputIter>
    T operator()(InputIter&& begin, InputIter&& end)
    {
        return operator()(std::forward<InputIter>(begin),
                          std::forward<InputIter>(end), [](auto v) { return v; });
    }

    /**
     * [Kurtosis description]
     * @method Kurtosis
     */
    Kurtosis()
        : _n(0),
          _mean(Sum<T>()),
          _M2(Sum<T>()),
          _M3(Sum<T>()),
          _M4(Sum<T>()),
          _n1(0),
          _delta(0),
          _delta_n(0),
          _delta_n2(0),
          _term1(0)
    {
    }

    /**
     * Kurstosis Copy constructor
     * @method Kurtosis
     * @param  other    Object to copy from.
     */
    Kurtosis(const Kurtosis& other) = default;

    /**
     * Move Constructor
     * @method Kurtosis
     * @param  other    Object to move from
     */
    Kurtosis(Kurtosis&& other) = default;

    /**
     * [Kurtosis description]
     * @method Kurtosis
     */
    Kurtosis& operator=(Kurtosis&& other) = default;

    /**
     * [Kurtosis description]
     * @method Kurtosis
     */
    Kurtosis& operator=(const Kurtosis& other) = default;

    /**
     * Kurtosis destructor
     */
    virtual ~Kurtosis() = default;
};

/**
 * @brief Functor for computing Mode (most common value) of a stream of
 *        values. Be careful with large streams, because the records are
 *        stored in here.
 *
 */
template <typename T>
class Mode
{
protected:
    std::unordered_map<T, std::size_t> _counter;
    T _mode;

public:
    inline void reset()
    {
        _counter.clear();
        _mode = 0;
    }

    inline T result()
    {
        return *std::max_element(_counter.begin(), _counter.end(),
                                 [](auto&& value) { return value.second; });
    }

    void operator()(T value)
    {
        ++_counter[value];
    }

    template <typename Iter, typename Getter>
    T operator()(Iter&& begin, Iter&& end, Getter&& getter)
    {
        std::for_each(std::forward<Iter>(begin), std::forward<Iter>(end),
                      [&getter](auto&& value) { return getter(value); });
        return result();
    }

    template <typename Iter>
    T operator()(Iter&& begin, Iter&& end)
    {
        return operator()(std::forward<Iter>(begin), std::forward<Iter>(end),
                          [](auto&& value) { return value; });
    }

    Mode() : _counter(std::unordered_map<T, std::size_t>()), _mode(0)
    {
    }
    Mode(const Mode& other) = default;
    Mode(Mode&& other) = default;
    Mode& operator=(const Mode& other) = default;
    Mode& operator=(Mode&& other) = default;
    virtual ~Mode() = default;
};

/**
 * @brief Functor for computing the 'percent'-th quantile of a data stream.
 *        Be careful with large data, because the record is stored in the class.
 *        Partially sorts the array, hence best used if multiple quantiles of
 *        the same data are needed.
 *
 * @tparam T Datatype to use
 * @tparam percent percentage giving the quantile to compute, e.g. 30 -> 30% quantile
 * @tparam reservesize the number of elements the internal buffer shall allocate memory for.
 *                     defaults to 10000.
 */
template <typename T, std::size_t percent, std::size_t reservesize = 10000>
class Quantile
{
protected:
    std::vector<T> _data;

public:
    inline T result()
    {
        // FIXME: Interpolation method needed here!
        return *std::nth_element(
            _data.begin(),
            _data.begin() + std::floor(_data.size() * (static_cast<double>(percent) / 100)),
            _data.end());
    }

    inline void reset()
    {
        _data.clear();
    }

    T operator()(T value)
    {
        _data.push_back(value);
    }

    template <typename Iter, typename Getter>
    T operator()(Iter&& begin, Iter&& end, Getter&& getter)
    {
        std::for_each(std::forward<Iter>(begin), std::forward<Iter>(end),
                      [&getter, this](auto&& v) { _data.push_back(getter(v)); });
        return result();
    }

    template <typename Iter>
    T operator()(Iter&& begin, Iter&& end)
    {
        _data.insert(_data.end(), std::forward<Iter>(begin), std::forward<Iter>(end));
        return result();
    }

    Quantile() : _data(std::vector<T>())
    {
        _data.reserve(reservesize);
    }
    Quantile(const Quantile& other) = default;
    Quantile(Quantile&& other) = default;
    Quantile& operator=(const Quantile& other) = default;
    Quantile& operator=(Quantile&& other) = default;
    virtual ~Quantile() = default;
};

template <typename T, std::size_t... quantiles>
class OnlineQuantile
{
protected:
    std::array<std::size_t, sizeof...(quantiles)> _quantiles;

public:
    OnlineQuantile() = default;
    OnlineQuantile(const OnlineQuantile& other) = default;
    OnlineQuantile(OnlineQuantile&& other) = default;
    OnlineQuantile& operator=(const OnlineQuantile& other) = default;
    OnlineQuantile& operator=(OnlineQuantile&& other) = default;
    virtual ~OnlineQuantile() = default;
};

} // namespace Statistics

#endif
