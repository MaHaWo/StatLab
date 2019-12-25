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

struct Fix
{

    std::vector< double > u = std::vector< double >(100000);
    std::vector< double > n = std::vector< double >(100000);

    Fix()
    {
        std::string       line = "";
        std::stringstream s;

        // read stuff
        for (auto&& [vec, name] : { std::make_pair(&n, "normal.csv"),
                                    std::make_pair(&u, "uniform.csv") })
        {
            std::ifstream stream(name, std::ios::in);

            std::size_t i = 0;
            // FIXME: this is wrong
            while (std::getline(stream, line))
            {
                if(i > 100){
                    break;
                }
                double f = 0;
                s << line;
                std::cout << s.str() << ", ";
                s >> (*vec)[i++];
                s>> f;
                std::cout << (*vec)[i] << ", " << f << std::endl;
                s.str() = "";
            }
        }
    }

    ~Fix()          = default;
    Fix(const Fix&) = default;
    Fix(Fix&&)      = default;
    Fix&
    operator=(const Fix&) = default;
    Fix&
    operator=(Fix&&) = default;
};

BOOST_TEST_DECORATOR(*boost::unit_test::tolerance(2e-15));
BOOST_FIXTURE_TEST_CASE(sum_test, Fix)
{
    // numpy values
    // uniform: -44701.74983180444
    // normal: -53313.558700174326
    Sum< double, NaNPolicy::propagate> sum;


    double s = sum(u.begin(), u.begin()+100);
    BOOST_TEST(s == -44701.74983180444);

    s = sum(u.begin(), u.end());
    BOOST_TEST(s == 2 * -44701.74983180444);

    sum.reset();

    BOOST_TEST(sum.result() == 0.);


    // s = sum(n.begin(), n.end());
    // BOOST_TEST(s == -53313.558700174326);

    // s = sum(n.begin(), n.end());
    // BOOST_TEST(s == 2 * -53313.558700174326);
}