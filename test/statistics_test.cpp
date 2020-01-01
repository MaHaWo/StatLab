#include <stdexcept>
#define BOOST_TEST_MODULE UNARY_TESTS
#include <boost/test/included/unit_test.hpp>
#include <sstream>

#include "../include/statistics.hpp"

#include <fstream>
#include <random>
#include <string>
#include <vector>

using namespace Statistics;
using namespace std::literals::chrono_literals;

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

        u_inf_nan      = u;
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

BOOST_TEST_DECORATOR(*boost::unit_test::tolerance(2e-14));
BOOST_FIXTURE_TEST_CASE(sum_test, Fix)
{
    // numpy values for sum, obtained via partial summation algorithm
    // uniform: -44701.74983180444
    // normal: -53313.558700174326

    Sum< double, NaNPolicy::propagate > sum_propagate;
    Sum< double, NaNPolicy::skip >      sum_skip;
    Sum< double, NaNPolicy::error >     sum_throw;

    sum_propagate.reset();
    BOOST_TEST(sum_propagate.result() == 0);

    // value based operator
    for (auto&& it = u.begin(); it != u.end(); ++it)
    {
        sum_propagate(*it);
    }

    BOOST_TEST(sum_propagate.result() == -44701.74983180444);

    // partition sum
    sum_propagate.reset();
    sum_propagate(u.begin(), u.end());
    BOOST_TEST(sum_propagate.result() == -44701.74983180444);

    sum_propagate.reset();
    sum_propagate(u_inf_nan.begin(), u_inf_nan.end());
    BOOST_TEST(std::isnan(sum_propagate.result()));

    sum_skip(u.begin(), u.end());
    BOOST_TEST(sum_skip.result() == -44701.74983180444);

    try{
        sum_throw(u_inf_nan.begin(), u_inf_nan.end());
    }
    catch(std::invalid_argument& e){
        BOOST_TEST(e.what() == "std::invalid_argument: Error, nan or inf found: nan");
    }
}

BOOST_TEST_DECORATOR(*boost::unit_test::tolerance(2e-14));
BOOST_FIXTURE_TEST_CASE(sum_kahan_test, Fix)
{
    // numpy values for sum, obtained via partial summation algorithm
    // uniform: -44701.74983180444
    // normal: -53313.558700174326

    SumKahan< double, NaNPolicy::propagate > sum_propagate;
    SumKahan< double, NaNPolicy::skip >      sum_skip;
    SumKahan< double, NaNPolicy::error >     sum_throw;

    sum_propagate.reset();
    BOOST_TEST(sum_propagate.result() == 0);

    // value based operator
    for (auto&& it = u.begin(); it != u.end(); ++it)
    {
        sum_propagate(*it);
    }

    BOOST_TEST(sum_propagate.result() == -44701.74983180444);
    
    sum_propagate.reset();
    BOOST_TEST(sum_propagate.result() == 0);
}