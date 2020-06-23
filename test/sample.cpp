#include <fstream>
#include <iomanip>
#include <random>
#include <vector>

int main() {
  std::mt19937_64 rng(6384239);
  std::uniform_real_distribution<double> dist(-1000, 1000);
  std::normal_distribution<double> ndist(0, 250);
    std::vector<double> x(100000);
    std::vector<double> n(100000);
  std::generate(x.begin(), x.end(), [&]() { return dist(rng); });
  std::generate(n.begin(), n.end(), [&]() { return ndist(rng); });

  std::ofstream uniform("uniform.csv", std::ios::trunc);
  std::ofstream normal("normal.csv", std::ios::trunc);

  for(std::size_t i = 0; i < 100000; ++i){
      uniform << std::setprecision(16) << x[i] << std::endl;
      normal << std::setprecision(16) << n[i] << std::endl;
  }

uniform.close();
normal.close();

  return 0;
}