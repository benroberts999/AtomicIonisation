#include "Matrix.hpp"
#include "PhysicsConstants.hpp"
#include "StandardHaloModel.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace Scattering {

//==============================================================================
std::tuple<std::vector<double>, std::vector<double>, Matrix<double>>
read_in_Kion(const std::string &in_filename) {

  // This is very sensitive to input file format; not robust, but it works
  std::ifstream in_file{in_filename};
  if (!in_file.is_open()) {
    std::cout << "Failed to open " << in_filename << '\n';
    return {{}, {}, {}};
  }

  std::vector<std::string> lines;
  std::string line;
  while (std::getline(in_file, line)) {
    if (line.empty() || line.front() == '#' || line.front() == '\n')
      continue;
    lines.push_back(line);
  }

  std::vector<double> q_grid, E_grid;

  // Read in the E grid (assumed to be in atomic units):
  std::stringstream ss_E{lines.at(0)};
  double temp_E;
  while (ss_E >> temp_E) {
    E_grid.push_back(temp_E);
  }

  // Read in the q grid (assumed to be in atomic units):
  std::stringstream ss_q{lines.at(1)};
  double temp_q;
  while (ss_q >> temp_q) {
    q_grid.push_back(temp_q);
  }

  // first two 'lines' are q and E grids
  assert(lines.size() - 2 == E_grid.size());

  Matrix<double> Kion(E_grid.size(), q_grid.size());

  // Read in the Kion matrix:
  for (std::size_t iE = 0; iE < E_grid.size(); ++iE) {
    std::stringstream ss_K{lines.at(iE + 2)}; // first two are q and E grids
    std::size_t iq = 0;
    while (iq < q_grid.size() && ss_K >> Kion(iE, iq)) {
      ++iq;
    }
    assert(iq == q_grid.size());
  }

  return {E_grid, q_grid, Kion};
}

//==============================================================================
double dsdE_E(const Matrix<double> &Kion, std::size_t iE, double E, double v,
              double mx, const std::function<double(double)> &Fx2,
              const std::vector<double> &q_grid) {
  assert(iE < Kion.rows());
  assert(Kion.cols() == q_grid.size());

  const auto factor = 1.0 / (2.0 * v * v);

  // assume q grid is logarithmic:
  // dq = (dq/dt) * dt
  // (dq/dt) = q // logarithmic
  // dt = log(maxq / minq) / (num_q - 1);
  const auto maxq = *std::max_element(q_grid.begin(), q_grid.end());
  const auto minq = *std::min_element(q_grid.begin(), q_grid.end());
  const auto dt = std::log(maxq / minq) / double(q_grid.size() - 1);

  const auto arg = std::pow(mx * v, 2) - 2.0 * mx * E;
  // conservation of energy:
  if (arg < 0.0)
    return 0.0;

  const auto qminus = mx * v - std::sqrt(arg);
  const auto qplus = mx * v + std::sqrt(arg);
  // conservation of momentum:
  if (qminus > maxq || qplus < minq)
    return 0.0;

  // Integrate over q:
  double dsdE = 0;
  for (std::size_t iq = 0; iq < q_grid.size(); iq++) {
    const auto q = q_grid.at(iq);
    // conservation of momentum:
    if (q < qminus || q > qplus)
      continue;
    const auto dqdt = q; // for log grid only
    dsdE += q * dqdt * Fx2(q) * Kion(iE, iq);
  } // q int
  dsdE *= factor * dt;

  return dsdE;
}

//==============================================================================
double dsvdE_E(const Matrix<double> &Kion, std::size_t iE, double E, double mx,
               const std::function<double(double)> &Fx2,
               const std::vector<double> &q_grid,
               const std::function<double(double)> &Fv, int num_v) {
  assert(iE < Kion.rows());
  assert(Kion.cols() == q_grid.size());

  // max v (for integration), in atomic units:
  const auto vmax =
      1.2 * Astro::Constants::v_max / Physics::Conversions::V_to_kms;
  // factor 1.2 just for safety: to ensure we encompass max escape velocity,
  // even for strange velocity distributions
  const double dv = vmax / num_v;
  double tmp_dsvdE_E = 0.0;
  for (int iv = 0; iv < num_v; ++iv) {
    const double v = (iv + 1) * dv;
    tmp_dsvdE_E += v * Fv(v) * dsdE_E(Kion, iE, E, v, mx, Fx2, q_grid);
  }
  return tmp_dsvdE_E * dv;
}

//==============================================================================
std::vector<double> dsvdE(const Matrix<double> Kion, double mx,
                          const std::function<double(double)> &Fx2,
                          const std::vector<double> &E_grid,
                          const std::vector<double> &q_grid,
                          const std::function<double(double)> &Fv, int num_v) {
  assert(Kion.rows() == E_grid.size());
  assert(Kion.cols() == q_grid.size());

  std::vector<double> dsvdE_vec;
  dsvdE_vec.reserve(E_grid.size());
  for (std::size_t iE = 0; iE < E_grid.size(); ++iE) {
    const auto E = E_grid.at(iE);
    const auto tmp = dsvdE_E(Kion, iE, E, mx, Fx2, q_grid, Fv, num_v);
    dsvdE_vec.push_back(tmp);
  }
  return dsvdE_vec;
}

//==============================================================================
double gaussian(double sigma, double x) {
  const auto root_2pi = std::sqrt(2.0 * M_PI);
  const auto arg = (x / sigma) * (x / sigma);
  return std::exp(-0.5 * arg) / (sigma * root_2pi);
}

//==============================================================================
std::vector<double>
convolvedRate(const std::vector<double> &dRdE,
              const std::function<double(double, double)> &f,
              const std::vector<double> &E_grid) {
  assert(dRdE.size() == E_grid.size());

  std::vector<double> conv;
  conv.reserve(E_grid.size());

  const auto Emin = *std::min_element(E_grid.begin(), E_grid.end());
  const auto Emax = *std::max_element(E_grid.begin(), E_grid.end());
  // This assumes a logarithmic E grid:
  const auto dtE = std::log(Emax / Emin) / double(E_grid.size() - 1);

  for (std::size_t iE = 0; iE < E_grid.size(); iE++) {
    const auto E = E_grid.at(iE);

    double tmp = 0.0;
    for (std::size_t iEp = 0; iEp < E_grid.size(); iEp++) {
      const auto Ep = E_grid.at(iEp);
      const auto dEdt = E; // log grid
      tmp += f(E, Ep) * dRdE[iEp] * dEdt;
    }
    tmp *= dtE;
    conv.push_back(tmp);
  }
  return conv;
}

//==============================================================================
double binnedCount(const std::vector<double> &dSdE, double Ea, double Eb,
                   const std::vector<double> &E_grid) {
  assert(dSdE.size() == E_grid.size());

  const auto Emin = *std::min_element(E_grid.begin(), E_grid.end());
  const auto Emax = *std::max_element(E_grid.begin(), E_grid.end());
  // This assumes a logarithmic E grid:
  const auto dt = std::log(Emax / Emin) / double(E_grid.size() - 1);

  double tmp = 0.0;
  for (std::size_t iE = 0; iE < E_grid.size(); iE++) {
    // nb: for sparse E grids, this is unstable (particularly at large E)!
    const auto E = E_grid.at(iE);
    if (E < Ea || E > Eb)
      continue;
    const auto dEdt = E; // log grid
    tmp += dSdE.at(iE) * dEdt;
  }
  return tmp * dt / (Eb - Ea);
}

//==============================================================================
double sigtot_E(const Matrix<double> &Kion, const std::vector<double> &E_grid,
                const std::vector<double> &q_grid, double Ei) {
  assert(Kion.rows() == E_grid.size());
  assert(Kion.cols() == q_grid.size());

  double sigtot = 0.0;

  const auto Emin = *std::min_element(E_grid.begin(), E_grid.end());
  const auto Emax = *std::max_element(E_grid.begin(), E_grid.end());
  const auto qmin = *std::min_element(q_grid.begin(), q_grid.end());
  const auto qmax = *std::max_element(q_grid.begin(), q_grid.end());

  const auto dtE = std::log(Emax / Emin) / double(E_grid.size() - 1);
  const auto dtQ = std::log(qmax / qmin) / double(q_grid.size() - 1);

  for (std::size_t iE = 0; iE < E_grid.size(); iE++) {
    auto E = E_grid.at(iE);
    if (E >= Ei)
      continue;
    const double qminus = std::sqrt(2.0 * Ei) - std::sqrt(2.0 * (Ei - E));
    const double qplus = std::sqrt(2.0 * Ei) + std::sqrt(2.0 * (Ei - E));
    for (std::size_t iq = 0; iq < q_grid.size(); iq++) {
      const double q = q_grid.at(iq);
      if (q < qminus || q > qplus)
        continue;
      const double dEdq = E * q;
      const double FX2 = 1.0 / (q * q * q);
      sigtot += dEdq * FX2 * Kion(iE, iq);
    }
  }
  sigtot *= dtE * dtQ;
  return (4.0 * M_PI / Ei) * sigtot;
}

//==============================================================================
std::vector<double> sigtot(const Matrix<double> &Kion,
                           const std::vector<double> &E_grid,
                           const std::vector<double> &q_grid) {
  assert(Kion.rows() == E_grid.size());
  assert(Kion.cols() == q_grid.size());

  std::vector<double> sig_vec;
  sig_vec.reserve(E_grid.size());
  for (const auto &E : E_grid) {
    sig_vec.push_back(sigtot_E(Kion, E_grid, q_grid, E));
  }
  return sig_vec;
}
} // namespace Scattering