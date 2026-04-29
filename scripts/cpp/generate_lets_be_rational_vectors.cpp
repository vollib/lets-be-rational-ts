#include <cmath>
#include <cfloat>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "lets_be_rational.h"
#include "normaldistribution.h"

namespace {

struct NormalisedCase {
  double x;
  double s;
  double q;
};

struct BlackCase {
  double F;
  double K;
  double sigma;
  double T;
  double q;
};

std::string decimal(double value) {
  if (std::isnan(value)) return "null";
  if (std::isinf(value)) return value > 0 ? "1e9999" : "-1e9999";
  std::ostringstream out;
  out << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
  return out.str();
}

std::string hex(double value) {
  std::ostringstream out;
  out << std::hexfloat << value;
  return out.str();
}

void field(const char* name, double value, bool last = false) {
  std::cout << "      \"" << name << "\": " << decimal(value) << (last ? "\n" : ",\n");
}

void hex_field(const char* name, double value, bool last = false) {
  std::cout << "      \"" << name << "\": \"" << hex(value) << "\"" << (last ? "\n" : ",\n");
}

} // namespace

int main() {
  const std::vector<double> zs = {-10.0, -8.0, -4.0, -1.0, -0.1, 0.0, 0.1, 1.0, 4.0, 8.0, 10.0};
  const std::vector<NormalisedCase> normalised_cases = {
    {-10.0, 0.05, 1.0},
    {-4.0, 0.5, 1.0},
    {-1.0, 0.25, 1.0},
    {-0.25, 0.05, 1.0},
    {-0.01, 0.001, 1.0},
    {0.0, 0.001, 1.0},
    {0.0, 0.2, 1.0},
    {0.01, 0.001, 1.0},
    {0.25, 0.05, 1.0},
    {1.0, 0.25, 1.0},
    {4.0, 0.5, 1.0},
    {10.0, 0.05, 1.0},
    {-10.0, 0.05, -1.0},
    {-4.0, 0.5, -1.0},
    {-1.0, 0.25, -1.0},
    {-0.25, 0.05, -1.0},
    {-0.01, 0.001, -1.0},
    {0.0, 0.001, -1.0},
    {0.0, 0.2, -1.0},
    {0.01, 0.001, -1.0},
    {0.25, 0.05, -1.0},
    {1.0, 0.25, -1.0},
    {4.0, 0.5, -1.0},
    {10.0, 0.05, -1.0},
    {-0.5, 4.0, 1.0},
    {0.5, 4.0, -1.0},
    {-64.0, 8.0, 1.0},
    {64.0, 8.0, -1.0}
  };
  const std::vector<BlackCase> black_cases = {
    {100.0, 100.0, 0.2, 0.5, 1.0},
    {100.0, 100.0, 0.2, 0.5, -1.0},
    {100.0, 95.0, 0.3, 0.5, 1.0},
    {100.0, 95.0, 0.3, 0.5, -1.0},
    {95.0, 100.0, 0.3, 0.5, 1.0},
    {95.0, 100.0, 0.3, 0.5, -1.0},
    {1000.0, 100.0, 0.1, 2.0, 1.0},
    {100.0, 1000.0, 0.1, 2.0, -1.0},
    {1e-8, 1e-8, 2.0, 0.01, 1.0},
    {1e8, 1e8, 0.01, 30.0, -1.0}
  };

  std::cout << "{\n";
  std::cout << "  \"provenance\": {\n";
  std::cout << "    \"source\": \"Peter Jaeckel LetsBeRational C++\",\n";
  std::cout << "    \"decimalPrecision\": \"std::numeric_limits<double>::max_digits10\",\n";
  std::cout << "    \"hexFloatIncluded\": true\n";
  std::cout << "  },\n";

  std::cout << "  \"normalDistribution\": [\n";
  for (std::size_t i = 0; i < zs.size(); ++i) {
    const double z = zs[i];
    std::cout << "    {\n";
    field("z", z);
    field("normCdf", norm_cdf(z));
    hex_field("normCdfHex", norm_cdf(z), true);
    std::cout << "    }" << (i + 1 == zs.size() ? "\n" : ",\n");
  }
  std::cout << "  ],\n";

  std::cout << "  \"normalisedBlack\": [\n";
  for (std::size_t i = 0; i < normalised_cases.size(); ++i) {
    const NormalisedCase c = normalised_cases[i];
    const double beta = NormalisedBlack(c.x, c.s, c.q);
    std::cout << "    {\n";
    std::cout << "      \"input\": {\n";
    field("x", c.x);
    field("s", c.s);
    field("q", c.q, true);
    std::cout << "      },\n";
    std::cout << "      \"output\": {\n";
    field("normalisedBlack", beta);
    field("normalisedImpliedVolatility", NormalisedImpliedBlackVolatility(beta, c.x, c.q));
    field("complementaryNormalisedBlack", ComplementaryNormalisedBlack(c.x, c.s));
    field("normalisedVega", NormalisedVega(c.x, c.s));
    hex_field("normalisedBlackHex", beta, true);
    std::cout << "      }\n";
    std::cout << "    }" << (i + 1 == normalised_cases.size() ? "\n" : ",\n");
  }
  std::cout << "  ],\n";

  std::cout << "  \"black\": [\n";
  for (std::size_t i = 0; i < black_cases.size(); ++i) {
    const BlackCase c = black_cases[i];
    const double price = Black(c.F, c.K, c.sigma, c.T, c.q);
    std::cout << "    {\n";
    std::cout << "      \"input\": {\n";
    field("F", c.F);
    field("K", c.K);
    field("sigma", c.sigma);
    field("T", c.T);
    field("q", c.q, true);
    std::cout << "      },\n";
    std::cout << "      \"output\": {\n";
    field("black", price);
    field("impliedBlackVolatility", ImpliedBlackVolatility(price, c.F, c.K, c.T, c.q));
    field("vega", Vega(c.F, c.K, c.sigma, c.T));
    hex_field("blackHex", price, true);
    std::cout << "      }\n";
    std::cout << "    }" << (i + 1 == black_cases.size() ? "\n" : ",\n");
  }
  std::cout << "  ]\n";
  std::cout << "}\n";

  return 0;
}
