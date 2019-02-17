/**
 * @todo - device test case for floating point accuracy
 * @todo - check pairwise summation -> something feels wrong there
 * @todo - get a way to automatically determine significant digits
           and adjust tolerance of comparison accordingly
 * @todo - compare with boost -> is there anyway
 */

#include "../include/statistics.hh"
#include "../include/test_utils.hh"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

using namespace Statistics;
using namespace std::literals::chrono_literals;

template <typename Distribution>
auto generate_sample(std::size_t samplesize, std::size_t seed, Distribution&& dist, bool write)
{
    using T = typename Distribution::result_type;
    std::mt19937 generator(seed);
    std::vector<T> data(samplesize);
    std::generate(data.begin(), data.end(), [&]() { return dist(generator); });

    if (write)
    {
        std::fstream file("data.bin", std::ios::out | std::ios::binary);
        file.write(reinterpret_cast<char*>(data.data()),
                   sizeof(typename decltype(data)::value_type) * samplesize);
    }

    return data;
}

int main()
{
    auto data = generate_sample(10000, 75843, std::normal_distribution<double>(0, 5), true);

    // test summation
    Sum<double> sum_pairwise;
    SumKahan<double> sum_kahan;

    auto pairwisesum = sum_pairwise(data.begin(), data.end());
    auto kahansum = sum_kahan(data.begin(), data.end());
    double pythonsum = 429.451407544296; // put fucking floating point to 1.....
    double pythonkahansum = 429.4514075442959; // put fucking floating point to 1.....

    ASSERT_EQ_CUSTOM(pairwisesum, pythonsum, 2.4e-13);
    ASSERT_EQ_CUSTOM(kahansum, pythonkahansum, 1e-16);

    // test mean
    ArithmeticMean<double, Sum> pairwise_mean;
    ArithmeticMean<double, SumKahan> kahan_mean;

    double pairwise_arithmeticmean = pairwise_mean(data.begin(), data.end());
    double kahan_arithmeticmean = kahan_mean(data.begin(), data.end());
    double python_arithmeticmean = 0.0429451407544296;
    double python_kahanmean = 0.04294514075442959;
    double test_mean = std::accumulate(data.begin(), data.end(), 0.) / data.size();
    ASSERT_EQ_CUSTOM(pairwise_arithmeticmean, python_arithmeticmean, 1e-16);
    ASSERT_EQ_CUSTOM(kahan_arithmeticmean, python_kahanmean, 1e-16);

    // Test Variance
    Variance<double, SumKahan> variance;
    double var = variance(data.begin(), data.end());

    double pythonvariance = 25.636281809383632;
    double testvariance = std::accumulate(data.begin(), data.end(), 0.,
                                          [&](auto&& r, auto&& v) {
                                              return r + (v - test_mean) * (v - test_mean);
                                          }) /
                          (data.size() - 1);

    // FIXME: find reason for the large differences here!
    ASSERT_EQ_CUSTOM(var, pythonvariance, 2.5e-14);
    ASSERT_EQ_CUSTOM(var, testvariance, 5.4e-14);

    // Test skewness
    Skewness<double, SumKahan> skewness;
    double pythonskewness = -0.0057856842695107645;
    double skew = skewness(data.begin(), data.end());
    ASSERT_EQ_CUSTOM(skew, pythonskewness, 5e-16);

    // Test kurtosis
    Kurtosis<double, Sum> kurtosis_pairwise;
    Kurtosis<double, SumKahan> kurtosis;

    double pythonkurtosis = 0.01342924877361451;
    double kur = kurtosis(data.begin(), data.end());
    double pkur = kurtosis_pairwise(data.begin(), data.end());
    ASSERT_EQ_CUSTOM(kur, pythonkurtosis, 4e-15);
    ASSERT_EQ_CUSTOM(pkur, pythonkurtosis, 3e-14);

    return 0;
}
