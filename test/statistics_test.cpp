
#define BOOST_TEST_MODULE UNARY_TESTS

#include <boost/test/included/unit_test.hpp>

#include <exception>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../include/statistics.hpp"

using namespace Statistics;
using namespace std::literals;

struct Thing
{
    double      v;
    std::string s;
};

struct Fix
{
    std::vector< double > u         = std::vector< double >(100000);
    std::vector< double > n         = std::vector< double >(100000);
    std::vector< double > u_inf_nan = std::vector< double >(100000);
    std::vector< Thing >  strct     = std::vector< Thing >(100000);

    Fix()
    {
        std::string       line = "";
        std::stringstream s;

        // read stuff
        for (auto&& [vec, name] : { std::make_pair(&n, "normal.csv"),
                                    std::make_pair(&u, "uniform.csv") })
        {
            std::ifstream stream(name, std::ios::in);

            for (std::size_t i = 0; std::getline(stream, line); ++i)
            {
                std::stringstream s(line);
                s >> (*vec)[i];
                strct[i].s = line;
                strct[i].v = (*vec)[i];
            }
        }

        u_inf_nan = u;
        u_inf_nan.push_back(std::numeric_limits< double >::quiet_NaN());
        u_inf_nan.push_back(std::numeric_limits< double >::infinity());
    }

    ~Fix()          = default;
    Fix(const Fix&) = default;
    Fix(Fix&&)      = default;
    Fix&
    operator=(const Fix&) = default;
    Fix&
    operator=(Fix&&) = default;
};

BOOST_TEST_DECORATOR(*boost::unit_test::tolerance(2e-13));
BOOST_FIXTURE_TEST_CASE(sum_test, Fix)
{
    // numpy values for sum, obtained via partial summation algorithm
    double normal = 856.8614085806706;

    SumPairwise< double, NaNPolicy::propagate > sum_propagate;
    SumPairwise< double, NaNPolicy::skip >      sum_skip;
    SumPairwise< double, NaNPolicy::error >     sum_throw;

}

BOOST_TEST_DECORATOR(*boost::unit_test::tolerance(2e-15));
BOOST_FIXTURE_TEST_CASE(sum_kahan_test, Fix)
{
    // numpy values for sum, obtained via kahan summation algorithm
    // stolen from rosetta code: 
    // 
    // def kahansum(input): 
    //     summ = c = 0 
    //     for num in input: 
    //         y = num - c 
    //         t = summ + y 
    //         c = (t - summ) - y 
    //         summ = t 
    //     return summ 
    // 
    double normal = 856.8614085806706;

    SumKahan< double, NaNPolicy::propagate > sum_propagate;
    SumKahan< double, NaNPolicy::skip >      sum_skip;
    SumKahan< double, NaNPolicy::error >     sum_throw;
}