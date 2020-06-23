#define BOOST_TEST_MODULE StatLab_test

// this library
#include "../include/statistics.hpp"

// boost deps
#include <boost/test/unit_test.hpp>

// stl deps
#include <fstream>
#include <numeric>
#include <vector>

namespace utf = boost::unit_test;

struct Fixture
{
    std::vector< double > data;

    Fixture()
    {
        data.reserve(10000);

        std::ifstream normal_file("normal_0_1_data.csv", std::ios::in);
        double        x;
        while (file >> x)
        {
            data.push_back(x);
        }
        data.shrink_to_fit();
    }

    ~Fixture() = default;
};

BOOST_FIXTURE_TEST_CASE(Summation_test, Fixture)
{
    // desired value, derived from Python numpy
    double python_sum = -88.28237087766323;

    double sum_kahan_iter = StatLab::SumKahan< double >()(
        data.begin(), data.end(), StatLab::Utils::Identity< double >());

    double sum_pairwise_iter = StatLab::SumPairwise< double >()(
        data.begin(), data.end(), StatLab::Utils::Identity< double >());

    StatLab::SumPairwise< double > p_acc;
    StatLab::SumKahan< double >    k_acc;
    double                         sum_pairwise = 0;
    for (auto v : data)
    {
        sum_pairwise = p_acc(v);
    }

    double sum_kahan = 0;
    for (auto v : data)
    {
        sum_kahan = k_acc(v);
    }

    BOOST_TEST(python_sum == sum_kahan_iter,
               boost::test_tools::tolerance(7e-16));
    BOOST_TEST(python_sum == sum_kahan, boost::test_tools::tolerance(7e-16));

    BOOST_TEST(python_sum == sum_pairwise,
               boost::test_tools::tolerance(7.2e-15));
    BOOST_TEST(python_sum == sum_pairwise_iter,
               boost::test_tools::tolerance(2e-15));
}

BOOST_FIXTURE_TEST_CASE(Product_test, Fixture)
{
}

BOOST_FIXTURE_TEST_CASE(Moment_test, Fixture)
{
    // results from python
    double first_moment_python  = -0.008828237087766323;
    double second_moment_python = 1.0027303051670418;

    // first and second moments using kahan summation algorithm
    double first_moment_kahan_iter =
        StatLab::Moment< double, StatLab::SumKahan, 1 >()(
            data.begin(), data.end(), StatLab::Utils::Identity< double >());
    double second_moment_kahan_iter =
        StatLab::Moment< double, StatLab::SumKahan, 2 >()(
            data.begin(), data.end(), StatLab::Utils::Identity< double >());

    // first and second moment using pairwise summation
    double first_moment_pairwise_iter =
        StatLab::Moment< double, StatLab::SumPairwise, 1 >()(
            data.begin(), data.end(), StatLab::Utils::Identity< double >());
    double second_moment_pairwise_iter =
        StatLab::Moment< double, StatLab::SumPairwise, 2 >()(
            data.begin(), data.end(), StatLab::Utils::Identity< double >());

    // FIXME: these are wrong!

    std::cout << "datasize: " << data.size() << std::endl;

    BOOST_TEST(first_moment_python == first_moment_kahan_iter,
               boost::test_tools::tolerance(2e-15));
    BOOST_TEST(second_moment_python == second_moment_kahan_iter,
               boost::test_tools::tolerance(2e-15));

    BOOST_TEST(first_moment_python == first_moment_pairwise_iter,
               boost::test_tools::tolerance(2e-15));
    BOOST_TEST(second_moment_python == second_moment_pairwise_iter,
               boost::test_tools::tolerance(4.9e-15));

    // needed for testing elementwise computation
    StatLab::Moment< double, StatLab::SumKahan, 1 >    first_moment_kahan;
    StatLab::Moment< double, StatLab::SumKahan, 2 >    second_moment_kahan;
    StatLab::Moment< double, StatLab::SumPairwise, 1 > first_moment_pairwise;
    StatLab::Moment< double, StatLab::SumPairwise, 2 > second_moment_pairwise;

    std::size_t i = 0;
    for (auto& v : data)
    {
        first_moment_kahan(v);
        second_moment_kahan(v);
        first_moment_pairwise(v);
        second_moment_pairwise(v);
    }

    BOOST_TEST(first_moment_python == first_moment_kahan.result(),
               boost::test_tools::tolerance(5.9e-16));

    BOOST_TEST(second_moment_python == second_moment_kahan.result(),
               boost::test_tools::tolerance(2e-16));

    BOOST_TEST(first_moment_python == first_moment_pairwise.result(),
               boost::test_tools::tolerance(7.1e-15));

    BOOST_TEST(second_moment_python == second_moment_pairwise.result(),
               boost::test_tools::tolerance(2e-15)); // FIXME: this is wrong!!
}
