#include "../statistics.hh"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <thread>
#include <vector>
using namespace Statistics;
using namespace std::literals::chrono_literals;

template <typename T>
void to_csv(std::string filename, std::string header,
            std::vector<T>& container) {
    std::ofstream file(filename, std::ios::out | std::ios::trunc);
    file << header << "\n";
    for (auto& v : container) {
        file << std::setprecision(16) << v << "\n";
    }
}

template <typename T> std::vector<T> from_csv(std::string filename) {
    std::ifstream file(filename, std::ios::in);
    // get number of lines in file
    std::string line;
    std::vector<T> data;
    std::size_t number_of_lines = 0;
    while (std::getline(file, line)) {
        ++number_of_lines;
    }
    file.clear();
    file.seekg(0, std::ios::beg);
    data.resize(number_of_lines - 1); // header ignore
    std::getline(file, line);         // header ignore
    for (std::size_t i = 0; i < data.size(); ++i) {
        std::getline(file, line);
        std::istringstream s(line);
        s >> data[i];
    }
    return data;
}

template <typename Iter> auto moment(int order, Iter&& begin, Iter&& end) {
    double m = std::accumulate(begin, end, 0., [&](auto& lhs, auto& rhs) {
        return lhs + std::pow(rhs, order);
    });
    return m / double(std::distance(begin, end));
}

template <typename Iter> auto mean(Iter&& begin, Iter&& end) {
    return moment(1, std::forward<Iter>(begin), std::forward<Iter>(end));
}

template <typename Iter> auto variance(Iter&& begin, Iter&& end) {
    double m = mean(std::forward<Iter>(begin), std::forward<Iter>(end));
    double v = std::accumulate(begin, end, 0., [&](auto& lhs, auto& rhs) {
        return lhs + std::pow(rhs - m, 2);
    });
    return v / double(std::distance(begin, end));
}

int main() {
    // std::this_thread::sleep_for(20s);
    std::default_random_engine rng(347285);
    std::normal_distribution<double> n(0., 1.);
    std::vector<double> v(1000000);
    std::generate(v.begin(), v.end(), [&]() { return n(rng); });
    to_csv("data.csv", "value", v);
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
              << " ak: " << ak(v.begin(), v.end()) << " \n"
              << std::setprecision(16) << std::setw(18)
              << " naive: " << mean(v.begin(), v.end()) << std::endl;

    Stddev<double, SumNaive> stdn;
    Stddev<double, SumPairwise> stdp;
    Stddev<double, SumKahan> stdk;

    std::cout << std::setprecision(16) << std::setw(18)
              << " stdn: " << stdn(v.begin(), v.end()) << "\n"
              << std::setprecision(16) << std::setw(18)
              << " stdp: " << stdp(v.begin(), v.end()) << " \n"
              << std::setprecision(16) << std::setw(18)
              << " stdk: " << stdk(v.begin(), v.end()) << " \n"
              << std::setprecision(16) << std::setw(18)
              << " naive: " << std::pow(variance(v.begin(), v.end()), 0.5)
              << std::endl;

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

    ArithmeticMean<double, SumPairwise> amp;
    for (auto& value : v) {
        amp(std::forward<double>(value));
    }
    std::cout << "arithmetic mean: " << amp.result() << std::endl;
    amp.reset();
    std::cout << "arithmetic mean using range: " << amp(v.begin(), v.end())
              << std::endl;
    Statistician<ArithmeticMean<double, SumPairwise>,
                 Variance<double, SumPairwise>>
        statistician;

    auto mean_var = statistician(v.begin(), v.end());
    std::cout << "Arithmetic Mean: " << std::endl;
    std::cout << " singular    : "
              << ArithmeticMean<double, SumPairwise>()(v.begin(), v.end())
              << "\n"
              << " statistician: " << mean_var[0]
              << "\n naive       : " << mean(v.begin(), v.end()) << std::endl;

    std::cout << "Variance: " << std::endl;
    std::cout << " singular    : "
              << Variance<double, SumPairwise>()(v.begin(), v.end()) << "\n"
              << " statistician: " << mean_var[1]
              << "\n naive       : " << variance(v.begin(), v.end())
              << std::endl;

    using T = double;
    ArithmeticMean<T, SumNaive> am;
    Variance<T, SumNaive> vm;
    Skewness<T, SumNaive> sm;
    Kurtosis<T, SumNaive> km;

    Statistician<ArithmeticMean<T, SumPairwise>, Variance<T, SumPairwise>,
                 Skewness<T, SumPairwise>, Kurtosis<T, SumPairwise>>
        new_statistician;

    auto statres = new_statistician(v.begin(), v.end());
    auto statres_i =
        std::vector<T>{am(v.begin(), v.end()), vm(v.begin(), v.end()),
                       sm(v.begin(), v.end()), km(v.begin(), v.end())};

    std::array<std::string, 4> functionnames{
        {"mean", "variance", "skewness", "kurtosis"}};

    std::cout << "Test Statistician against singular functions" << std::endl;

    for (std::size_t i = 0; i < statres.size(); ++i) {
        std::cout << functionnames[i] << "\n"
                  << std::setprecision(16) << std::setw(16)
                  << "statistician: " << statres[i] << "\n"
                  << std::setw(16) << " singular: " << statres_i[i]
                  << std::endl;
        std::cout << std::endl;
    }
    return 0;
}
