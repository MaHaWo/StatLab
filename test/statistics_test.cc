#include "../statistics.hh"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <thread>
#include <vector>
using namespace Statistics;
using namespace std::literals::chrono_literals;

template <typename T>
std::vector<T> generate_data(bool read, std::size_t size = 1000000, bool write = false)
{
    if (!read)
    {
        std::default_random_engine rng(347285);
        std::normal_distribution<T> n(0., 1.);
        std::vector<T> v(size);
        std::generate(v.begin(), v.end(), [&]() { return n(rng); });
        if (write)
        {
            std::ofstream file(
                "data.bin", std::ios::out | std::ios::binary | std::ios::trunc);
            file.write(reinterpret_cast<char*>(&size), sizeof(std::size_t));
            file.write(reinterpret_cast<char*>(v.data()), sizeof(T) * size);
        }
        return v;
    }
    else
    {
        std::vector<T> data;
        std::ifstream file("data.bin", std::ios::in | std::ios::binary);
        std::size_t s = 0;
        file.read(reinterpret_cast<char*>(&s), sizeof(std::size_t));
        data.resize(s);
        file.read(reinterpret_cast<char*>(data.data()), s * sizeof(T));
        return data;
    }
}

template <typename T>
auto read_test_results(std::string filename, std::size_t size)
{
    std::ifstream file(filename, std::ios::in | std::ios::binary);
    std::vector<T> data(size);
    file.read(reinterpret_cast<char*>(data.data()), size * sizeof(T));
    return data;
}

template <template <typename> class Sum>
using Mean = ArithmeticMean<double, Sum>;
template <template <typename> class Sum>
using Std = Stddev<double, Sum>;
using Min = Minimum<double>;
using Max = Maximum<double>;
template <template <typename> class Sum>
using Skew = Skewness<double, Sum>;
template <template <typename> class Sum>
using EKur = ExcessKurtosis<double, Sum>;

int main()
{
    // std::this_thread::sleep_for(20s);
    auto v = generate_data<double>(false, 1000000, true);
    auto ddsv = read_test_results<double>("description.bin", 10);
    ddsv = decltype(ddsv){ddsv[1], ddsv[2], ddsv[3], ddsv[7], ddsv[8], ddsv[9]};
    std::vector<std::string> names{"mean", "std",      "min",
                                   "max",  "skewness", "excess_kurtosis"};
    Statistician<Mean<SumNaive>, Std<SumNaive>, Min, Max, Skew<SumNaive>, EKur<SumNaive>> StatN;
    Statistician<Mean<SumPairwise>, Std<SumPairwise>, Min, Max, Skew<SumPairwise>, EKur<SumPairwise>> StatP;
    Statistician<Mean<SumKahan>, Std<SumKahan>, Min, Max, Skew<SumKahan>, EKur<SumKahan>> StatK;

    double m = std::accumulate(v.begin(), v.end(), 0.) / v.size();
    double s = std::accumulate(
                   v.begin(), v.end(), 0.,
                   [&m](auto& a, auto& b) { return a + std::pow(b - m, 2.); }) /
               (v.size() - 1);

    double min = *std::min_element(v.begin(), v.end());
    double max = *std::max_element(v.begin(), v.end());
    double sk = 0.;
    double ek = 0.;

    std::vector<double> simplevalues{m, s, min, max, sk, ek};
    StatN(v.begin(), v.end());
    StatP(v.begin(), v.end());
    StatK(v.begin(), v.end());

    std::cout << std::setw(25) << "name" << std::setw(25) << "Simple"
              << std::setw(25) << "Pandas" << std::setw(25) << "Naive"
              << std::setw(25) << "Pairwise" << std::setw(25) << "Kahan" << std::endl;

    for (std::size_t i = 0; i < ddsv.size(); ++i)
    {
        std::cout << std::setw(25) << names[i] << std::setprecision(16)
                  << std::setw(25) << simplevalues[i] << std::setw(25)
                  << ddsv[i] << std::setw(25) << StatN.result()[i]
                  << std::setw(25) << StatP.result()[i] << std::setw(25)
                  << StatK.result()[i] << std::endl;
    }

    Mean<SumNaive> mn;
    Mean<SumPairwise> mp;
    Mean<SumKahan> mk;

    Std<SumNaive> stdn;
    Std<SumPairwise> stdp;
    Std<SumKahan> stdk;

    mn(v.begin(), v.end());
    mp(v.begin(), v.end());
    mk(v.begin(), v.end());

    stdn(v.begin(), v.end());
    stdp(v.begin(), v.end());
    stdk(v.begin(), v.end());

    std::cout << std::setw(25) << "Simple" << std::setw(25) << "Pandas"
              << std::setw(25) << "Naive" << std::setw(25) << "Pairwise"
              << std::setw(25) << "Kahan" << std::endl;
    std::cout << std::setw(25) << m << std::setw(25) << std::setw(25) << ddsv[0]
              << std::setw(25) << mn.result() << std::setw(25) << mp.result()
              << std::setw(25) << mk.result() << std::endl;
    std::cout << std::setw(25) << s << std::setw(25) << ddsv[1] << std::setw(25)
              << stdn.result() << std::setw(25) << stdp.result()
              << std::setw(25) << stdk.result() << std::endl;
    return 0;
}
