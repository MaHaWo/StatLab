#include "../statistics.hh"
#include <iomanip>
#include <iostream>
#include <random>
#include <thread>
#include <vector>
using namespace Statistics;
using namespace std::literals::chrono_literals;
int main() {
    // std::this_thread::sleep_for(20s);
    std::default_random_engine rng(347285);
    std::normal_distribution<double> n(0., 1.);
    std::vector<double> v(1000000);
    std::generate(v.begin(), v.end(), [&]() { return n(rng); });

    SumNaive<double> sn;
    SumPairwise<double> sp;
    SumKahan<double> sk;

    std::cout << std::setprecision(16) << std::setw(18)
              << " sn: " << sn(v.begin(), v.end()) << "\n"
              << std::setprecision(16) << std::setw(18)
              << " sp: " << sp(v.begin(), v.end()) << " \n"
              << std::setprecision(16) << std::setw(18)
              << " sk: " << sk(v.begin(), v.end()) << std::endl;

    ArithmeticMean<double, SumNaive> an;
    ArithmeticMean<double, SumPairwise> ap;
    ArithmeticMean<double, SumKahan> ak;
    std::cout << std::setprecision(16) << std::setw(18)
              << " an: " << an(v.begin(), v.end()) << "\n"
              << std::setprecision(16) << std::setw(18)
              << " ap: " << ap(v.begin(), v.end()) << " \n"
              << std::setprecision(16) << std::setw(18)
              << " ak: " << ak(v.begin(), v.end()) << std::endl;

    Stddev<double, SumNaive> stdn;
    Stddev<double, SumPairwise> stdp;
    Stddev<double, SumKahan> stdk;

    std::cout << std::setprecision(16) << std::setw(18)
              << " stdn: " << stdn(v.begin(), v.end()) << "\n"
              << std::setprecision(16) << std::setw(18)
              << " stdp: " << stdp(v.begin(), v.end()) << " \n"
              << std::setprecision(16) << std::setw(18)
              << " stdk: " << stdk(v.begin(), v.end()) << std::endl;

    Skewness<double, SumNaive> skn;
    Skewness<double, SumPairwise> skp;
    Skewness<double, SumKahan> skk;

    std::cout << std::setprecision(16) << std::setw(18)
              << " skn: " << skn(v.begin(), v.end()) << "\n"
              << std::setprecision(16) << std::setw(18)
              << " skp: " << skp(v.begin(), v.end()) << " \n"
              << std::setprecision(16) << std::setw(18)
              << " skk: " << skk(v.begin(), v.end()) << std::endl;

    Kurtosis<double, SumNaive> kn;
    Kurtosis<double, SumPairwise> kp;
    Kurtosis<double, SumKahan> kk;

    std::cout << std::setprecision(16) << std::setw(18)
              << " kn: " << kn(v.begin(), v.end()) << "\n"
              << std::setprecision(16) << std::setw(18)
              << " kp: " << kp(v.begin(), v.end()) << " \n"
              << std::setprecision(16) << std::setw(18)
              << " kk: " << kk(v.begin(), v.end()) << std::endl;

    return 0;
}