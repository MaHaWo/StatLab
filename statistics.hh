#ifndef UTOPIA_MODELS_AMEEMULTI_STATISTICS_HH
#define UTOPIA_MODELS_AMEEMULTI_STATISTICS_HH
/** @file statistics.hh
 *  @brief Implements functions for computing statistics.
 *
 *  @author Harald Mack, harald.mack@iup.uni-heidelberg.de
 *  @bug Quantiles return only exact quantiles, that is only
 *       values which are contained in the container.
 *  @todo Optimize, perhaps use template metaprogramming where appropriate
 *  @todo Use partition sums to make numerical errors smaller, test this!
 *  @todo Implement interpolation for quantiles, since they are defined in the
 * respective way!
 *
 */

#include "utils.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <random>
#include <type_traits>
#include <unordered_map>
#include <valarray>
#include <vector>

using namespace Statistics::Utils;
namespace Statistics {
template <typename T> struct SumNaive {
  private:
    T _s;

  public:
    using value_type = T;
    /**
     * @brief Construct a new Sum object
     *
     */
    SumNaive() : _s(T()) {}
    SumNaive(T&& v) : _s(v) {}
    SumNaive(const SumNaive& other) = default;
    SumNaive(SumNaive&& other) = default;
    SumNaive& operator=(const SumNaive& other) = default;
    SumNaive& operator=(SumNaive&& other) = default;
    virtual ~SumNaive() = default;

    /**
     * @brief Reset the internal state
     *
     */
    void reset() { _s = T(); }

    /**
     * @brief get result of computation
     *
     * @return T
     */
    T result() { return _s; }

    /**
     * @brief
     *
     *             Ignores NaN and Inf values
     *
     * @param[in]  begin          The begin iterator
     * @param[in]  end            The end iterator
     * @param[in]  f              function to apply to the values to accumulate.
     *
     *
     * @return     sum over (f(it)) for 'it' in [begin, end)
     */

    template <typename InputIterator, typename Getter>
    inline T operator()(InputIterator&& begin, InputIterator&& end,
                        Getter&& getter) {
        _s += std::accumulate(
            begin, end, value_type(),
            [&](auto&& lhs, auto&& rhs) { return lhs + getter(rhs); });
        return _s;
    }

    /**
     * @brief
     *
     * @tparam InputIterator
     * @param begin
     * @param end
     * @return T
     */
    template <typename InputIterator>
    inline T operator()(InputIterator&& begin, InputIterator&& end) {
        this->operator()(std::forward<InputIterator>(begin),
                         std::forward<InputIterator>(end),
                         [](auto&& value) { return value; });
        return _s;
    }

    /**
     * @brief updating the sum with a single value
     *
     * @param value
     * @return T
     */
    T operator()(T&& value) {
        _s += value;
        return _s;
    }
};
/**
 * @brief Functor for pairwise summation
 *
 * @tparam T return type
 */
template <typename T> struct SumPairwise {
  private:
    T _s;

  public:
    using value_type = T;
    /**
     * @brief Construct a new Sum object
     *
     */
    SumPairwise() : _s(T()) {}
    SumPairwise(T&& v) : _s(v) {}

    SumPairwise(const SumPairwise& other) = default;
    SumPairwise(SumPairwise&& other) = default;
    SumPairwise& operator=(const SumPairwise& other) = default;
    SumPairwise& operator=(SumPairwise&& other) = default;
    virtual ~SumPairwise() = default;

    /**
     * @brief Reset the internal state
     *
     */
    void reset() { _s = T(); }

    /**
     * @brief get result of computation
     *
     * @return T
     */
    T result() { return _s; }

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
    inline value_type operator()(InputIterator&& begin, InputIterator&& end,
                                 Getter&& getter) {
        if (begin != end) {
            if (std::distance(begin, end) <= 1000) {
                value_type s = std::accumulate(
                    begin, end, value_type(),
                    [&](auto&& lhs, auto&& rhs) { return lhs + getter(rhs); });

                _s += s;
            } else {
                std::size_t t =
                    std::floor(double(std::distance(begin, end) / 2));

                this->operator()(
                    std::forward<InputIterator>(begin),
                    std::next(std::forward<InputIterator>(begin), t),
                    std::forward<Getter>(getter));

                this->operator()(
                    std::next(std::forward<InputIterator>(begin), t),
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
    inline value_type operator()(InputIterator&& begin, InputIterator&& end) {
        if (begin != end) {
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
    inline T operator()(T&& value) {
        _s += value;
        return _s;
    }
};

/**
 * @brief Functor for the kahan summation algorithm
 *
 * @tparam T
 */
template <typename T> struct SumKahan {
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
    SumKahan() : _y(T()), _t(T()), _s(T()), _comp(T()) {}

    /**
     * @brief Construct a new Sum Kahan object
     *
     * @param v initializer
     */
    SumKahan(T&& v) : _y(v), _t(v), _s(v), _comp(v) {}

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
    void reset() {
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
    T result() { return _s; }
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
    inline T operator()(InputIterator&& begin, InputIterator&& end,
                        Getter&& getter) {
        if (begin != end) {
            _comp = 0.;
            for (; begin != end; ++begin) {
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
    inline T operator()(InputIterator&& begin, InputIterator&& end) {
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
    inline T operator()(T&& value) {
        _y = value - _comp;
        _t = _s + _y;
        _s = _t;
        return _s;
    }
};

/**
 * @brief Functor for computing 'order'-th moments wrt zero
 *
 * @tparam order Order of the desired moment
 * @tparam Summation Summation algorithm, defaults to 'Sum', but can also be
 * 'SumKahan'
 */
template <typename T, int order,
          template <typename> class Summation = SumPairwise>
struct Moment {
  private:
    T _mmnt;
    T _n;
    Summation<T> _sum;

  public:
    using value_type = T;
    Moment() : _mmnt(T()), _n(T()), _sum(Summation<T>()) {}
    virtual ~Moment() = default;
    Moment(const Moment& other) = default;
    Moment(Moment&& other) = default;
    Moment& operator=(const Moment& other) = default;
    Moment& operator=(Moment&& other) = default;

    /**
     * @brief get current value of the  order-th moment wrt zero.
     *
     * @return double
     */
    double result() { return _mmnt / _n; }

    /**
     * @brief reset all internal members including result
     *
     */
    void reset() {
        _mmnt = 0;
        _n = 0;
        _sum.reset();
    }

    /**
     * @brief      Computes the 'order'-th moment of the range [begin, end)
     *
     * @param[in]  begin          The begin iterator
     * @param[in]  end            The end iterator
     * @param[in]  getter         A function taking an argument of type
     * InputIterator::value_type and returning a member of it or modifying
     * it in another way.
     *
     *
     * @return    order-th moment wrt zero.
     */
    template <typename InputIterator, typename Getter>
    inline T operator()(InputIterator&& begin, InputIterator&& end,
                        Getter&& getter) {
        if (begin != end) {
            _n += std::distance(begin, end);
            _mmnt += _sum(std::forward<InputIterator>(begin),
                          std::forward<InputIterator>(end),
                          [&getter](auto&& value) -> value_type {
                              return std::pow(getter(value), order);
                          });
        }
        return _mmnt / _n;
    }

    /**
     * @brief      Computes the 'order'-th moment of the range [begin, end)
     *
     * @param[in]  order          The order
     * @param[in]  begin          The begin iterator
     * @param[in]  end            The end iterator
     *
     *
     *
     * @return    order-th moment wrt zero.
     */
    template <typename InputIterator>
    inline T operator()(InputIterator&& begin, InputIterator&& end) {
        return this->operator()(std::forward<InputIterator>(begin),
                                std::forward<InputIterator>(end),
                                [](auto&& value) { return value; });
    }

    /**
     * @brief update the current value of the moment with the value given in
     * the argument
     *
     * @tparam Value Some type convertible to double
     * @param value Value to update the moment with
     * @return double updated  order-th moment wrt zero.
     */
    inline T operator()(T&& value) {
        ++_n;
        _mmnt = _sum(std::pow(value, order));
        return _mmnt / _n;
    }
};

/**
 * @brief Functor for computing arithmetic mean
 *
 */
template <typename T, template <typename> class Summation = SumPairwise>
struct ArithmeticMean : Moment<T, 1, Summation> {
    ArithmeticMean() = default;
    virtual ~ArithmeticMean() = default;
    ArithmeticMean(const ArithmeticMean& other) = default;
    ArithmeticMean(ArithmeticMean&& other) = default;
    ArithmeticMean& operator=(const ArithmeticMean& other) = default;
    ArithmeticMean& operator=(ArithmeticMean&& other) = default;
};

/**
 * @brief Functor for computing Harmonic mean
 *
 */
template <typename T, template <typename> class Summation = SumPairwise>
struct HarmonicMean {
  private:
    T _hm;
    T _n;
    Summation<T> _sum;

  public:
    using value_type = T;
    HarmonicMean() : _hm(T()), _n(T()), _sum(Summation<T>()) {}
    virtual ~HarmonicMean() = default;
    HarmonicMean(const HarmonicMean& other) = default;
    HarmonicMean(HarmonicMean&& other) = default;
    HarmonicMean& operator=(const HarmonicMean& other) = default;
    HarmonicMean& operator=(HarmonicMean&& other) = default;

    void reset() {
        _hm = 0;
        _n = 0;
        _sum.reset();
    }

    inline T result() { return _n / _hm; }
    /**
     * @brief      Function computes the harmonic mean of the range [begin,
     * end).
     *
     * @param[in]  begin          The begin iterator
     * @param[in]  end            The end iterator
     * @param[in]  getter         The getter function extracting an element
     * from the object the iterator points to
     *
     *
     * @return     The harmonic mean of the range begin,end
     */

    template <typename InputIterator, typename Getter>
    inline T operator()(InputIterator&& begin, InputIterator&& end,
                        Getter&& getter) {
        if (begin != end) {
            _n += std::distance(begin, end);
            _hm += _sum(std::forward<InputIterator>(begin),
                        std::forward<InputIterator>(end),
                        [&](T&& value) { return 1. / getter(*value); });
        }
        return _n / _hm;
    }

    /**
     * @brief      Computes the harmonic mean of the values in the range
     * [begin, end)
     *
     * @param[in]  begin          The begin iterator of the range
     * @param[in]  end            The end iterator of the range
     *
     * @tparam     InputIterator  Automatically determined template
     * parameter
     *
     * @return     harmonic mean of [begin, end)
     */
    template <typename InputIterator>
    inline T operator()(InputIterator&& begin, InputIterator&& end) {
        return this->operator()(std::forward<InputIterator>(begin),
                                std::forward<InputIterator>(end),
                                [](auto& value) { return value; });
    }

    /**
     * @brief
     *
     * @tparam Value
     * @param value
     * @return double
     */
    inline T operator()(T&& value) {
        ++_n;
        _hm += _sum(1. / value);
        return _n / _hm;
    }
};

// FIXME: create online algorithm for this!
template <int order> struct CentralMoment {
    // FIXME: put in online algorithm here here
};

/**
 * @brief Functor for computing variance
 *
 * @tparam T
 * @tparam Sum
 */
template <typename T, template <typename> class Summation = SumPairwise>
struct Variance {
  private:
    // using V = std::valarray<T>;
    T _n;
    T _mean;
    T _M2;
    T _d;
    // V _values;
    Summation<T> _sum;
    // enum { n = 0, mean = 1, M2 = 2 };

  public:
    using value_type = T;
    Variance()
        : _n(T()), _mean(T()), _M2(T()), _d(T()), _sum(Summation<T>(T())) {}
    virtual ~Variance() = default;
    Variance(const Variance& other) = default;
    Variance(Variance&& other) = default;
    Variance& operator=(const Variance& other) = default;
    Variance& operator=(Variance&& other) = default;

    inline T result() { return _M2 / (_n - 1.); }

    void reset() {
        _n = T();
        _mean = T();
        _M2 = T();
        _d = T();
        _sum.reset();
    }
    /**
     * @brief      Computes the sample-variance of the range [begin, end).
     *
     * @param[in]  begin      The iterator pointing to the first element in
     *                        [begin, end)
     * @param[in]  end        The iterator pointing to the element after the
     * last element in the range [begin, end)
     * @param[in]  getter     Function converting an iterator to a value.
     *
     *
     * @return    The sample-variance of the range [begin, end)
     */
    template <typename InputIterator, typename Getter>
    inline T operator()(InputIterator&& begin, InputIterator&& end,
                        Getter&& getter) {
        // this should be changed to returning a numeric vector [_n,
        // _mean, _M2] or even tuple
        _M2 += _sum(std::forward<InputIterator>(begin),
                    std::forward<InputIterator>(end), [&](auto& value) -> T {
                        _n += 1;
                        _d = getter(value) - _mean;
                        _mean += _d / _n;
                        return _d * (getter(value) - _mean); // m2
                    });

        return _M2 / (_n - 1);
    }

    /**
     * @brief      Computes the sample-variance of the range [begin, end).
     *
     * @param[in]  begin      The iterator pointing to the first element in
     *                        [begin, end)
     * @param[in]  end        The iterator pointing to the element after the
     * last element in the range [begin, end)
     * @param[in]  getter     Function converting an iterator to a value.
     *
     *
     * @return    The sample-variance of the range [begin, end)
     */
    template <typename InputIterator>
    inline T operator()(InputIterator&& begin, InputIterator&& end) {
        return operator()(std::forward<InputIterator>(begin),
                          std::forward<InputIterator>(end),
                          [](auto& value) { return value; });
    }

    /**
     * @brief For updating the current value of the valriance
     *
     * @tparam Value Some type which is convertible to double
     * @param value actual value to add
     * @return double updated mean
     */
    inline T operator()(T&& value) {
        _n += 1;
        _d = value - _mean;
        _mean += _d / _n;
        _M2 += _d * (value - _mean); // m2
        return _M2 / (_n - 1);
    }
};

/**
 * @brief Functor for computing the standard deviation
 *
 * @tparam T
 * @tparam Summation
 */
template <typename T, template <typename> class Summation = SumPairwise>
struct Stddev : Variance<T, Summation> {
  public:
    using value_type = T;
    using Base = Variance<T, Summation>;
    Stddev() : Base() {}
    virtual ~Stddev() = default;
    Stddev(const Stddev& other) = default;
    Stddev(Stddev&& other) = default;
    Stddev& operator=(const Stddev& other) = default;
    Stddev& operator=(Stddev&& other) = default;

    inline T result() {
        return std::pow(static_cast<Base&>(*this).result(), 0.5);
    }

    /**
     * @brief      Computes the sample-variance of the range [begin, end).
     *
     * @param[in]  begin      The iterator pointing to the first element in
     *                        [begin, end)
     * @param[in]  end        The iterator pointing to the element after the
     * last element in the range [begin, end)
     * @param[in]  getter     Function converting an iterator to a value.
     *
     *
     * @return    The sample-variance of the range [begin, end)
     */
    // FIXME: use partition algorithm for the summation. Possible??
    template <typename InputIterator, typename Getter>
    inline T operator()(InputIterator&& begin, InputIterator&& end,
                        Getter&& getter) {
        // this should be changed to returning a numeric vector [_n,
        // _mean, _M2] or even tuple
        return std::pow(static_cast<Base&>(*this).operator()(
                            std::forward<InputIterator>(begin),
                            std::forward<InputIterator>(end),
                            std::forward<Getter>(getter)),
                        0.5);
    }

    /**
     * @brief      Computes the sample-variance of the range [begin, end).
     *
     * @param[in]  begin      The iterator pointing to the first element in
     *                        [begin, end)
     * @param[in]  end        The iterator pointing to the element after the
     * last element in the range [begin, end)
     * @param[in]  getter     Function converting an iterator to a value.
     *
     *
     * @return    The sample-variance of the range [begin, end)
     */
    // FIXME: use partition algorithm for the summation.
    template <typename InputIterator>
    inline T operator()(InputIterator&& begin, InputIterator&& end) {
        return this->operator()(std::forward<InputIterator>(begin),
                                std::forward<InputIterator>(end),
                                [](auto& value) { return value; });
    }

    /**
     * @brief For updating the current value of the valriance
     *
     * @tparam Value Some type which is convertible to double
     * @param value actual value to add
     * @return double updated mean
     */
    inline T operator()(T&& value) {
        return std::pow(
            static_cast<Base&>(*this).operator()(std::forward<T>(value)), 0.5);
    }
};

/**
 * @brief      Computes the sample-skewness of the range [begin, end).
 *
 * @param[in]  begin      The iterator pointing to the first element in
 *                        [begin, end)
 * @param[in]  end        The iterator pointing to the element after the
 * last element in the range [begin, end)
 * @param[in]  getter     Function converting an iterator to a value.
 *
 *
 *
 * @return     The sample-skewness of [begin, end)
 */
template <typename T, template <typename> class Summation = SumPairwise>
struct Skewness {
  private:
    T _n;
    T _mean;
    T _M2;
    T _M3;
    T _n1;
    T _delta;
    T _delta_n;
    T _term1;
    Summation<T> _sum;

  public:
    Skewness()
        : _n(T()), _mean(T()), _M2(T()), _M3(T()), _n1(T()), _delta(T()),
          _delta_n(T()), _term1(T()), _sum(Summation<T>()) {}
    virtual ~Skewness() = default;
    Skewness(const Skewness& other) = default;
    Skewness(Skewness&& other) = default;
    Skewness& operator=(const Skewness& other) = default;
    Skewness& operator=(Skewness&& other) = default;

    void reset() {
        _n = 0;
        _mean = 0;
        _M2 = 0;
        _M3 = 0;
        _n1 = 0;
        _delta = 0;
        _delta_n = 0;
        _term1 = 0;
        _sum.reset();
    }

    T result() { return std::pow(_n, 0.5) * _M3 / std::pow(_M2, 1.5); }

    template <typename InputIterator, typename Elementgetter>
    T operator()(InputIterator&& begin, InputIterator&& end,
                 Elementgetter&& getter) {
        _M3 += _sum(std::forward<InputIterator>(begin),
                    std::forward<InputIterator>(end), [&](auto& value) {
                        _n1 = _n;
                        _n += 1;
                        _delta = getter(value) - _mean;
                        _delta_n = _delta / _n;
                        _term1 = _delta * _delta_n * _n1;
                        _mean += _delta_n;

                        double M3 =
                            _term1 * _delta_n * (_n - 2) - 3 * _delta_n * _M2;
                        _M2 += _term1;
                        return M3;
                    });

        return std::pow(_n, 0.5) * _M3 / std::pow(_M2, 1.5);
    }

    /**
     * @brief      Computes the sample-skewness of the range [begin, end).
     *
     * @param[in]  begin      The iterator pointing to the first element in
     *                        [begin, end)
     * @param[in]  end        The iterator pointing to the element after the
     * last element in the range [begin, end)
     * @param[in]  getter     Function converting an iterator to a value.
     *
     *
     *
     * @return     The sample-skewness of [begin, end)
     */
    template <typename InputIterator>
    T operator()(InputIterator&& begin, InputIterator&& end) {
        return this->operator()(std::forward<InputIterator>(begin),
                                std::forward<InputIterator>(end),
                                [](auto& value) { return value; });
    }

    /**
     * @brief Adding a single value to the current stuff
     *
     * @param value
     * @return T
     */
    T operator()(T&& value) {
        _n1 = _n;
        _n += 1;
        _delta = value - _mean;
        _delta_n = _delta / _n;
        _term1 = _delta * _delta_n * _n1;
        _mean += _delta_n;

        _M3 += _term1 * _delta_n * (_n - 2) - 3 * _delta_n * _M2;
        _M2 += _term1;
        return std::pow(_n, 0.5) * _M3 / std::pow(_M2, 1.5);
    }
};

/**
 * @brief      Computes the sample-kurtosis of the range [begin, end).
 *
 * @param[in]  begin      The iterator pointing to the first element in
 *                        [begin, end)
 * @param[in]  end        The iterator pointing to the element after the
 * last element in the range [begin, end)
 * @param[in]  getter     Function converting an iterator to a value.
 *
 *
 *
 * @return     The sample-kurtosis [begin, end)
 */

template <typename T, template <typename> class Summation = SumPairwise>
struct Kurtosis {
  private:
    T _n;
    T _mean;
    T _M2;
    T _M3;
    T _M4;
    T _n1 = 0;
    T _delta;
    T _delta_n;
    T _delta_n2;
    T _term1;
    Summation<T> _sum;

  public:
    Kurtosis()
        : _n(T()), _mean(T()), _M2(T()), _M3(T()), _M4(T()), _n1(T()),
          _delta(T()), _delta_n(T()), _delta_n2(T()), _term1(T()),
          _sum(Summation<T>()) {}

    virtual ~Kurtosis() = default;
    Kurtosis(const Kurtosis& other) = default;
    Kurtosis(Kurtosis&& other) = default;
    Kurtosis& operator=(const Kurtosis& other) = default;
    Kurtosis& operator=(Kurtosis&& other) = default;

    void reset() {
        _n = 0;
        _mean = 0;
        _M2 = 0;
        _M3 = 0;
        _M4 = 0;
        _n1 = 0;
        _delta = 0;
        _delta_n = 0;
        _delta_n2 = 0;
        _term1 = 0;
        _sum.reset();
    }

    T result() { return _n * _M4 / (std::pow(_M2, 2)); }

    template <typename InputIterator, typename Elementgetter>
    double operator()(InputIterator&& begin, InputIterator&& end,
                      Elementgetter&& getter) {
        _M4 = _sum(
            std::forward<InputIterator>(begin),
            std::forward<InputIterator>(end), [&](auto&& value) {
                _n1 = _n;
                _n += 1;
                _delta = getter(value) - _mean;
                _delta_n = _delta / _n;
                _delta_n2 = _delta_n * _delta_n;
                _term1 = _delta * _delta_n * _n1;
                _mean += _delta_n;
                double M4 = _M4 + _term1 * _delta_n2 * (_n * _n - 3 * _n + 3) +
                            6 * _delta_n2 * _M2 - 4 * _delta_n * _M3;
                _M3 = _M3 + _term1 * _delta_n * (_n - 2) - 3 * _delta_n * _M2;
                _M2 = _M2 + _term1;
                return M4;
            });

        return _n * _M4 / (std::pow(_M2, 2));
    }

    /**
     * @brief      Computes the sample-kurtosis of the range [begin, end).
     *
     * @param[in]  begin      The iterator pointing to the first element in
     *                        [begin, end)
     * @param[in]  end        The iterator pointing to the element after the
     * last element in the range [begin, end)
     * @param[in]  getter     Function converting an iterator to a value.
     *
     *
     *
     * @return     The sample-kurtosis [begin, end)
     */
    template <typename InputIterator>
    double operator()(InputIterator&& begin, InputIterator&& end) {
        return operator()(std::forward<InputIterator>(begin),
                          std::forward<InputIterator>(end),
                          [](auto& value) { return value; });
    }

    /**
     * @brief Adding a single value to the kurtosis
     *
     * @param value
     * @return T
     */
    T operator()(T&& value) {
        _n1 = _n;
        _n += 1;
        _delta = value - _mean;
        _delta_n = _delta / _n;
        _delta_n2 = _delta_n * _delta_n;
        _term1 = _delta * _delta_n * _n1;
        _mean += _delta_n;
        _M4 += _term1 * _delta_n2 * (_n * _n - 3 * _n + 3) +
               6 * _delta_n2 * _M2 - 4 * _delta_n * _M3;
        _M3 += _term1 * _delta_n * (_n - 2) - 3 * _delta_n * _M2;
        _M2 += _term1;
        return _n * _M4 / (std::pow(_M2, 2));
    }
};

/**
 * @brief Functor for getting the Excess kurtosis
 *
 * @tparam T
 * @tparam Summation
 */
template <typename T, template <typename> class Summation = SumPairwise>
struct ExcessKurtosis : Kurtosis<T, Summation> {
  public:
    using Base = Kurtosis<T, Summation>;
    /**
     * @brief      Computes the excess-kurtosis of [begin, end)
     *
     * @param[in]  begin      The iterator pointing to the first element in
     *                        [begin, end)
     * @param[in]  end        The iterator pointing to the element after the
     * last element in the range [begin, end)
     * @param[in]  getter     Function converting an iterator to a value.
     *
     *
     *
     * @return     Excess-kurtosis of [begin, end), i.e.,
     * sample-kurtosis-3.0.
     */
    template <typename InputIterator, typename Elementgetter>
    double operator()(InputIterator&& begin, InputIterator&& end,
                      Elementgetter&& getter) {
        return static_cast<Base&>(*this).operator()(
                   std::forward<InputIterator>(begin),
                   std::forward<InputIterator>(end),
                   std::forward<Elementgetter>(getter)) -
               3.;
    }

    /**
     * @brief      Computes the excess-kurtosis of [begin, end)
     *
     * @param[in]  begin      The iterator pointing to the first element in
     *                        [begin, end)
     * @param[in]  end        The iterator pointing to the element after the
     * last element in the range [begin, end)
     * @param[in]  getter     Function converting an iterator to a value.
     *
     *
     *
     * @return     Excess-kurtosis of [begin, end), i.e.,
     * sample-kurtosis-3.0.
     */
    template <typename InputIterator>
    T operator()(InputIterator&& begin, InputIterator&& end) {
        return static_cast<Base&>(*this).operator()(std::forward<T>(begin),
                                                    std::forward<T>(end)) -
               3.;
    }

    /**
     * @brief Update excess kurtosis with a single value
     *
     * @param value
     * @return T
     */
    T operator()(T&& value) {
        return static_cast<Base&>(*this).operator()(std::forward<T>(value)) -
               3.;
    }
};

/**
 * @brief Standardizes the values in the range [begin,end), i.e. replacing the
 * values by their respective zscores. Is there a online possibility for this?
 *
 */
template <typename T, template <typename> class Summation = SumPairwise>
struct Standardize {
  private:
    std::vector<T> _zscore;

  public:
    auto& result() { return _zscore; }

    void reset() { return _zscore.clear(); }

    template <typename InputIterator, typename Getter>
    auto operator()(InputIterator&& begin, InputIterator&& end,
                    Getter&& getter) {
        T avg = ArithmeticMean<T, Summation>()(begin, end, getter);
        T std = Stddev<T, Summation>()(begin, end, getter);

        _zscore.resize(std::distance(begin, end));
        auto vbegin = _zscore.begin();
        for (; begin != end; ++begin, ++vbegin) {
            *vbegin = (static_cast<long double>(getter(*begin)) - avg) / std;
        }

        return _zscore;
    }

    template <typename InputIterator>
    auto operator()(InputIterator&& begin, InputIterator&& end) {
        return this->operator()(std::forward<InputIterator>(begin),
                                std::forward<InputIterator>(end),
                                [](auto& value) { return value; });
    }
};

/**
 * @brief      Computes the covariance between the two ranges
 *             [begin_first, end_first) and [begin_second, end_second).
 *             If one range is shorter than another the covariance is
 *             computed only for the overlab of the two ranges.
 * @param[in]  begin_first   The iterator pointing to the begin of the first
 *                           range.
 * @param[in]  end_first     The iterator pointing to the element after the
 *                            end  of the first range
 * @param[in]  begin_second  The iterator pointing to the begin of the second
 *                           range.
 * @param[in]  end_second    The iterator pointing to the element after the
 *                           end  of the secodn range
 * @param[in]  getter        Function for converting an iterator into a value.
 *
 *
 * @return     The covariance of the two ranges:
 *              cov([begin_first, end_first), [begin_second, end_second))
 */
template <typename T, template <typename> class Summation = SumPairwise>
struct Covariance {
  private:
    T _n;
    T _mean1;
    T _mean2;
    T _M12;
    T _d1;
    T _d2;
    Summation<T> _sum;

  public:
    void reset() {
        _n = 0;
        _mean1 = 0;
        _mean2 = 0;
        _M12 = 0;
        _d1 = 0;
        _d2 = 0;
        _sum.reset();
    }

    T result() { return _M12 * (_n / (_n - 1)); }

    template <typename InputIterator1, typename InputIterator2,
              typename Getter1, typename Getter2>
    double operator()(InputIterator1&& begin_first, InputIterator1&& end_first,
                      InputIterator2&& begin_second,
                      InputIterator2&& end_second, Getter1&& getter1,
                      Getter2&& getter2) {
        // runs over min(std::distance(begin_first, end_first),
        // begin_second, end_second))

        if (std::distance(begin_first, end_first) !=
            std::distance(begin_second, end_second)) {
            throw std::invalid_argument("Iterator ranges have to be of equal "
                                        "length for computing covariance");
        }

        auto it1 = begin_first;
        auto it2 = begin_second;

        _M12 += _sum(begin_first, end_first, [&, this](auto& valuefirst) {
            _n += 1;
            _d1 = (getter1(valuefirst) - _mean1) / _n;
            _mean1 += _d1;
            _d2 = (getter2(*it2) - _mean2) / _n;
            _mean2 += _d2;
            T M = (_n - 1) * _d1 * _d2 - _M12 / _n;
            ++it2;
            return M;
        });

        return _M12 * (_n / (_n - 1));
    }

    /**
     * @brief      Computes the covariance between the two ranges
     *             [begin_first, end_first) and [begin_second, end_second).
     *             If one range is shorter than another the covariance is
     *             computed only for the overlab of the two ranges.
     * @param[in]  begin_first   The iterator pointing to the begin of the
     * first range.
     * @param[in]  end_first     The iterator pointing to the element after
     * the end  of the first range
     * @param[in]  begin_second  The iterator pointing to the begin of the
     * second range.
     * @param[in]  end_second    The iterator pointing to the element after
     * the end  of the secodn range
     * @param[in]  getter        Function for converting an iterator into a
     * value.
     *
     *
     * @return     The covariance of the two ranges:
     *              cov([begin_first, end_first), [begin_second,
     * end_second))
     */
    template <typename InputIterator1, typename InputIterator2>
    double operator()(InputIterator1 begin_first, InputIterator1 end_first,
                      InputIterator2 begin_second, InputIterator2 end_second) {
        // runs over min(std::distance(begin_first, end_first),
        // begin_second, end_second))
        return operator()(begin_first, end_first, begin_second, end_second,
                          [](auto& value) { return value; },
                          [](auto& value) { return value; });
    }

    /**
     * @brief updating covariance with single values
     *
     * @param first
     * @param second
     * @return T
     */
    T operator()(T&& first, T&& second) {
        _n += 1;
        _d1 = (first - _mean1) / _n;
        _mean1 += _d1;
        _d2 = (second - _mean2) / _n;
        _mean2 += _d2;
        _M12 += _sum((_n - 1) * _d1 * _d2 - _M12 / _n);

        return _M12 * (_n / (_n - 1));
    }
};

/**
 * @brief Struct for computing the quantile
 *
 * @tparam T
 */

// FIXME: need to implement P2 algorithm here!
template <typename T, int percent, typename Comp = IsLess> struct Quantile {
  private:
    T _res;
    Comp _comp;
    static constexpr double _percent = percent;
    // all kinds of other stuff here
  public:
    Quantile() = default;
    Quantile(const Quantile& other) = default;
    Quantile(Quantile&& other) = default;
    Quantile& operator=(Quantile&& other) = default;
    Quantile& operator=(const Quantile& other) = default;

    void reset() { _res = 0; }

    T result() { return _res; }
    /**
     * @brief      Computes the 'percent'-th quantile of the distribution.
     *
     * @param[in]  percent    The percentage to compute the quantile for.
     *
     * @param[in]  begin      The iterator pointing to the first element in
     *                        [begin, end)
     * @param[in]  end        The iterator pointing to the element after the
     * last element in the range [begin, end)
     * @param[in]  comp       The comparison-function to use, defaults to
     *                        operator<(const T&, const T&).
     *
     *
     * @return    The 'percentage'-th quantile of the range [begin, end).
     */
    template <typename InputIterator, typename Getter>
    T operator()(InputIterator&& begin, InputIterator&& end, Getter&& getter) {
        std::vector<T> temp;
        T size = std::distance(begin, end);

        temp.reserve(size);
        for (auto it = begin; it != end; ++it) {
            temp.emplace_back(getter(*it));
        }
        T idx = size * (_percent / 100.);

        std::nth_element(temp.begin(), temp.begin() + idx, temp.end(), _comp);
        return temp[std::floor(idx)];
    }

    /**
     * @brief Computes the 'percentage'-th quantile of the distribution.
     *
     * @tparam InputIterator
     * @param percent
     * @param begin
     * @param end
     * @return double
     */
    template <typename InputIterator>
    T operator()(InputIterator&& begin, InputIterator&& end) {
        return operator()(std::forward<InputIterator>(begin),
                          std::forward<InputIterator>(end),
                          [](const auto& value1, const auto& value2) {
                              return value1 < value2;
                          },
                          [](auto& value) { return value; });
    }

    T operator()(T&& value) { return 0; }
};

/**
 * @brief Functor for computing the median of a collection of numbers
 *
 * @tparam T
 * @tparam IsLess
 */
template <typename T, typename Comp = IsLess>
struct Median : Quantile<T, 50, Comp> {
    using Base = Quantile<T, 50, Comp>;

    Median() = default;
    Median(const Median& other) = default;
    Median(Median&& other) = default;
    Median& operator=(Median&& other) = default;
    Median& operator=(const Median& other) = default;

    /**
     * @brief Update with a range [begin, end)
     *
     * @tparam InputIterator
     * @tparam Getter
     * @param begin
     * @param end
     * @param getter
     * @return T
     */
    template <typename InputIterator, typename Getter>
    T operator()(InputIterator&& begin, InputIterator&& end, Getter&& getter) {
        return static_cast<Base&>(*this).operator()(
            std::forward<InputIterator>(begin),
            std::forward<InputIterator>(end), std::forward<Getter>(getter));
    }

    /**
     * @brief Update with a range [begin, end)
     *
     * @tparam InputIterator
     * @tparam Getter
     * @param begin
     * @param end
     * @param getter
     * @return T
     */
    template <typename InputIterator>
    T operator()(InputIterator&& begin, InputIterator&& end) {
        return static_cast<Base&>(*this).operator()(
            std::forward<InputIterator>(begin),
            std::forward<InputIterator>(end));
    }

    /**
     * @brief update with single value
     *
     * @param value
     * @return T
     */
    T operator()(T&& value) {
        return static_cast<Base&>(*this).operator()(
            std::forward<T>(std::forward<T>(value)));
    }
};

/**
 * @brief Functor for computing Mode
 *
 */
template <typename T, class KeyEqual = IsEqual> struct Mode {
  private:
    std::pair<T, std::size_t> _modecount;
    std::unordered_map<T, std::size_t, std::hash<T>, KeyEqual> _counter;
    T _maxelement;
    std::size_t _maxcount;

  public:
    /**
     * @brief
     *
     */
    void reset() {
        _modecount.swap(std::pair<T, std::size_t>());
        _counter.swap(
            std::unordered_map<T, std::size_t, std::hash<T>, KeyEqual>());
        _maxelement = 0;
        _maxcount = 0;
    }

    /**
     * @brief
     *
     * @return auto
     */
    auto result() { return _maxelement; }

    auto get_mode_and_count() {
        return std::make_pair(_maxelement, _modecount);
    }
    /**
     * @brief      Computes the mode of the range [begin, end), i.e., the
     * most frequent value in the range.
     *
     * @param[in]  begin      The iterator pointing to the first element in
     *                        [begin, end)
     * @param[in]  end        The iterator pointing to the element after the
     * last element in the range [begin, end)
     * @param[in]  getter     Function converting an iterator to a value.
     *
     *
     * @return     The mode of the range [begin, end) as pair (element,
     * count)
     */
    template <typename InputIterator, typename Getter>
    T operator()(InputIterator&& begin, InputIterator&& end, Getter&& getter) {
        for (; begin != end; ++begin) {
            ++_counter[getter(*begin)];
            if (_counter[getter(*begin)] > _maxcount) {
                _maxelement = getter(*begin);
                _maxcount = _counter[getter(*begin)];
            }
        }

        return _maxelement;
    }

    /**
     * @brief      Computes the mode of the range [begin, end), i.e., the
     * most frequent value in the range.
     *
     * @param[in]  begin      The iterator pointing to the first element in
     *                        [begin, end)
     * @param[in]  end        The iterator pointing to the element after the
     * last element in the range [begin, end)
     * @param[in]  getter     Function converting an iterator to a value.
     *
     * @param[in]  comp       The comparison-function to use, defaults to
     *                        operator<(const T&, const T&) for the second
     * of the elements of a pair (T, int).
     *
     *
     * @return     The mode of the range [begin, end)
     */
    template <typename InputIterator>
    T operator()(InputIterator&& begin, InputIterator&& end) {
        return operator()(std::forward<InputIterator>(begin),
                          std::forward<InputIterator>(end),
                          [](auto& value) { return value; });
    }

    /**
     * @brief update the mode with a single value
     *
     * @param value
     * @return std::pair<T, inta>
     */
    T operator()(T&& value) {
        _counter[value] += 1;
        if (_counter[value] > _maxcount) {
            _maxelement = value;
            _maxcount = _counter[value];
        }

        return _maxelement;
    }
};

/**
 * @brief
 *
 * @tparam T
 * @tparam IsLess
 */
template <typename T, typename Comp = IsLess> struct Minimum {
  private:
    T _min;
    Comp _comp;

  public:
    T result() { return _min; }

    void reset() { _min = T(); }
    /**
     * @brief Computes the minimum of the range [begin, end).
     *
     * @tparam InputIterator
     * @tparam Getter
     * @param begin
     * @param end
     * @param getter
     * @return auto
     */
    template <typename InputIterator, typename Getter>
    auto operator()(InputIterator&& begin, InputIterator&& end,
                    Getter&& getter) {
        auto value = getter(*begin);
        for (; begin != end; ++begin) {
            if (_comp(value, _min)) {
                _min = value;
            }
        }
        return _min;
    }

    /**
     * @brief Computes the minimum of the range [begin, end).
     *
     * @tparam InputIterator
     * @param begin
     * @param end
     * @param comp
     * @return auto
     */
    template <typename InputIterator>
    auto operator()(InputIterator&& begin, InputIterator&& end) {
        return *std::min_element(std::forward<InputIterator>(begin),
                                 std::forward<InputIterator>(end), _comp);
    }

    /**
     * @brief
     *
     * @param value
     * @return T
     */
    T operator()(T&& value) {
        if (_comp(value, _min)) {
            _min = value;
        }

        return _min;
    }
};

/**
 * @brief Functor for computing the maximum
 *
 * @tparam T
 * @tparam IsLess
 */
template <typename T, typename Comp = IsLess> struct Maximum {
  private:
    T _max;
    Comp _comp;

  public:
    T result() { return _max; }

    void reset() { _max = T(); }
    /**
     * @brief Computes the maximum of the range [begin, end).
     *
     * @tparam InputIterator
     * @tparam std::function<bool(const T &, const T &)>
     * @param begin
     * @param end
     * @param adaptor
     * @return auto
     */
    template <typename InputIterator, typename Getter>
    T operator()(InputIterator&& begin, InputIterator&& end, Getter&& getter) {
        for (; begin != end; ++begin) {
            if (_comp(_max, getter(*begin))) {
                _max = getter(*begin);
            }
        }
        return _max;
    }

    /**
     * @brief Computes the maximum of the range [begin, end).
     *
     * @tparam InputIterator
     * @tparam std::function<bool(const T &, const T &)>
     * @param begin
     * @param end
     * @param comp
     * @return auto
     */
    template <typename InputIterator>
    T operator()(InputIterator&& begin, InputIterator&& end) {
        return this->operator()(std::forward<InputIterator>(begin),
                                std::forward<InputIterator>(end),
                                [](auto& value) { return value; });
    }

    /**
     * @brief update with single value
     *
     * @param value
     * @return T
     */
    T operator()(T&& value) {
        if (_comp(_max, value)) {
            _max = value;
        }
        return _max;
    }
};

/**
 * @brief
 *
 */
struct Sample {
  private:
    double _share;

  public:
    /**
     * @brief      Function to generate a random sample of relative
     *              size 'percent' from a range [begin, end).
     *
     * @param[in]  percent      The sample size in percent.
     * @param[in]  begin      Iterator to the first element to consider
     * @param[in]  end        Iterator pointing to the element after the
     * last element to consider.
     *
     *
     * @return     Sample containing 'percent'-percent of the elements in
     * [begin,end) with the elemenst contained being randomly chosen. The
     * output type is the same as the container which provides the input
     * iterators.
     */
    template <typename InputIterator, typename Getter>
    auto operator()(InputIterator&& begin, InputIterator&& end,
                    Getter&& getter) {
        std::size_t size = std::distance(begin, end);
        std::size_t shr = (_share / 100.) * size;
        std::vector<typename std::result_of<Getter>::type> smpl(shr);
        std::vector<std::size_t> indices(size);
        std::random_device rd;
        std::minstd_rand0 generator(rd());
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), generator);

        auto smpl_it = smpl.begin();
        for (std::size_t i = 0; i < shr; ++i) {
            *smpl_it = getter(std::next(begin, indices[i]));
            ++smpl_it;
        }
        return smpl;
    }

    /**
     * @brief      Function to generate a random sample of relative
     *              size 'percent' from a range [begin, end).
     *
     * @param[in]  percent      The sample size in percent.
     * @param[in]  begin      Iterator to the first element to consider
     * @param[in]  end        Iterator pointing to the element after the
     * last element to consider.
     *
     *
     * @return     Sample containing 'percent'-percent of the elements in
     * [begin,end) with the elemenst contained being randomly chosen. The
     * output type is the same as the container which provides the input
     * iterators.
     */
    template <typename InputIterator>
    auto operator()(InputIterator begin, InputIterator end) {
        return operator()(std::forward<InputIterator>(begin),
                          std::forward<InputIterator>(end),
                          [](auto& value) { return value; });
    }
};

/**
 * @brief
 *
 * @tparam T
 * @tparam Funcs
 */
// FIXME: this needs some way to convey the template argumets for the functions,
// such that I do not need them to be given explicitly all the time.
template <typename... Funcs> struct Statistician {
  public:
    using value_type =
        typename std::tuple_element_t<0, std::tuple<Funcs...>>::value_type;
    using ResArr = std::array<value_type, sizeof...(Funcs)>;

  private:
    std::tuple<Funcs...> _funcs;
    ResArr _res;

  public:
    /**
     * @brief
     *
     * @return auto
     */
    auto result() {
        _res = tuple_reduce<ResArr>(
            [&](auto& f) -> value_type { return f.result(); }, _funcs);
        return _res;
    }

    /**
     * @brief
     *
     */
    void reset() {
        tuple_reduce<void>([&](auto& f) { f.reset(); }, _funcs);
    }

    /**
     * @brief Construct a new Statistician object
     *
     */
    Statistician() : _funcs(std::make_tuple(Funcs()...)), _res(ResArr()) {}

    /**
     * @brief Construct a new Statistician object
     *
     * @param other
     */
    Statistician(const Statistician& other) = default;

    /**
     * @brief Construct a new Statistician object
     *
     * @param other
     */
    Statistician(Statistician&& other) = default;

    /**
     * @brief
     *
     * @param other
     * @return Statistician&
     */
    Statistician& operator=(const Statistician& other) = default;

    /**
     * @brief
     *
     * @param other
     * @return Statistician&
     */
    Statistician& operator=(Statistician&& other) = default;

    /**
     * @brief Destroy the Statistician object
     *
     */
    virtual ~Statistician() = default;

    /**
     * @brief make array of desired functions from Funcs
     *
     * @tparam InputIterator
     * @tparam Getter
     * @param begin
     * @param end
     * @param getter
     * @return auto
     */
    template <typename InputIterator, typename Getter>
    auto operator()(InputIterator&& begin, InputIterator&& end,
                    Getter&& getter) {
        for (; begin != end; ++begin) {
            tuple_reduce<void>([&](auto&& f) { f(getter(*begin)); }, _funcs);
        }

        _res = tuple_reduce<ResArr>(
            [](auto& f) -> value_type { return f.result(); }, _funcs);
        return _res;
    }

    /**
     * @brief make array of desired functions from Funcs
     *
     * @tparam InputIterator
     * @param begin
     * @param end
     * @return auto
     */
    template <typename InputIterator>
    auto operator()(InputIterator&& begin, InputIterator&& end) {
        return this->operator()(std::forward<InputIterator>(begin),
                                std::forward<InputIterator>(end),
                                [](auto& value) { return value; });
    }

    /**
     * @brief
     *
     * @tparam InputIterator1
     * @tparam InputIterator2
     * @tparam Getter1
     * @tparam Getter2
     * @param begin1
     * @param end1
     * @param begin2
     * @param end2
     * @param getter1
     * @param getter2
     * @return auto
     */
    template <typename InputIterator1, typename InputIterator2,
              typename Getter1, typename Getter2>
    auto operator()(InputIterator1&& begin1, InputIterator1&& end1,
                    InputIterator2&& begin2, InputIterator2&& end2,
                    Getter1&& getter1, Getter2&& getter2) {
        for (; begin1 != end1 && begin2 != end2; ++begin1, ++begin2) {
            tuple_reduce<void>(
                [&](auto& f) { f(getter1(*begin1), getter2(*begin2)); },
                _funcs);
        }

        _res = tuple_reduce<ResArr>(
            [](auto& f) -> value_type { return f.result(); }, _funcs);
        return _res;
    }

    /**
     * @brief
     *
     * @tparam InputIterator1
     * @tparam InputIterator2
     * @param begin1
     * @param end1
     * @param begin2
     * @param end2
     * @return auto
     */
    template <typename InputIterator1, typename InputIterator2>
    auto operator()(InputIterator1&& begin1, InputIterator1&& end1,
                    InputIterator2&& begin2, InputIterator2&& end2) {
        return this->operator()(std::forward<InputIterator1>(begin1),
                                std::forward<InputIterator1>(end1),
                                std::forward<InputIterator2>(begin2),
                                std::forward<InputIterator2>(end2),
                                [](auto& v) { return v; },
                                [](auto& v) { return v; });
    }

    /**
     * @brief Update all member functors with a single value
     *
     * @tparam T automatically determined
     * @param value Single value to update the functors with
     * @return auto
     */
    auto operator()(value_type&& value1, value_type&& value2) {
        tuple_reduce(
            [&](auto& f) {
                f(std::forward<value_type>(value1),
                  std::forward<value_type>(value2));
            },
            _funcs);
        _res = tuple_reduce<std::array>([](auto& f) { return f.result(); },
                                        _funcs);
        return _res;
    }
};

} // namespace Statistics

#endif
