#include "StandardHaloModel.hpp"
#include <LinAlg/Matrix.hpp>
#include <Physics/PhysConst_constants.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

const double V_to_kms = (PhysConst::c_SI / PhysConst::c) / 1000.0;
const double M_to_GeV = PhysConst::m_e_MeV / 1000.0;
const double au_to_eV = PhysConst::Hartree_eV;

//==============================================================================
//! Reads E grid, q grid, and Kion(E,q) from '_mat' formatted text file. Note:
//! assumes it was in ATOMIC units.
std::tuple<std::vector<double>, std::vector<double>, LinAlg::Matrix<double>>
read_in_Kion(const std::string &in_filename) {
  std::ifstream in_file{in_filename};
  std::vector<std::string> lines;
  std::string line;
  while (std::getline(in_file, line)) {
    if (line.empty() || line.front() == '#' || line.front() == '\n')
      continue;
    lines.push_back(line);
  }

  std::vector<double> q_grid, E_grid;

  std::stringstream ss_E{lines.at(0)};
  double temp_E;
  while (ss_E >> temp_E) {
    E_grid.push_back(temp_E);
  }

  std::stringstream ss_q{lines.at(1)};
  double temp_q;
  while (ss_q >> temp_q) {
    q_grid.push_back(temp_q);
  }

  // first two 'lines' are q and E grids
  assert(lines.size() - 2 == E_grid.size());

  LinAlg::Matrix<double> Kion(E_grid.size(), q_grid.size());

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
//! Calculates ds/dE for given energy deposition, E, and velocity, v,
//! in units of sigma-bar_e (v and E in atomic units)
double dsdE_E(const LinAlg::Matrix<double> &Kion, std::size_t iE, double E,
              double v, double mx, const std::function<double(double)> &Fx2,
              const std::vector<double> &q_grid) {

  assert(Kion.rows() > iE);
  assert(Kion.cols() == q_grid.size());

  const auto factor = 1.0 / (2.0 * v * v);

  // assume q grid is logarithmic:
  // dq = (dq/dt) * dt
  // (dq/dt) = q // logarithmic
  // du = log(maxq / minq) / (num_q - 1);
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
//! Calculates <ds.v>/dE for given energy deposition, E,
//! in units of sigma-bar_e (E in atomic units)
double dsvdE_E(const LinAlg::Matrix<double> &Kion, std::size_t iE, double E,
               double mx, const std::function<double(double)> &Fx2,
               const std::vector<double> &q_grid,
               const std::function<double(double)> &Fv, int num_v) {
  const auto vmax = AstroConsts::v_max / V_to_kms; // in au
  const double dv = vmax / num_v;
  double tmp_dsvdE_E = 0.0;
  for (int iv = 0; iv < num_v; ++iv) {
    double v = (iv + 1) * dv;
    const auto tmp_dsde = dsdE_E(Kion, iE, E, v, mx, Fx2, q_grid);
    tmp_dsvdE_E += v * Fv(v) * tmp_dsde;
  }
  return tmp_dsvdE_E * dv;
}

//==============================================================================
//! Calculates <ds.v>/dE as function of energy deposition, E,
//! in units of sigma-bar_e (E in atomic units)
std::vector<double> dsvdE(const LinAlg::Matrix<double> Kion, double mx,
                          const std::function<double(double)> &Fx2,
                          const std::vector<double> &E_grid,
                          const std::vector<double> &q_grid,
                          const std::function<double(double)> &Fv, int num_v) {
  assert(Kion.rows() == E_grid.size());
  assert(Kion.cols() == q_grid.size());

  std::vector<double> dsvdE_vec;
  dsvdE_vec.reserve(E_grid.size());
  for (std::size_t iE = 0; iE < E_grid.size(); ++iE) {
    auto E = E_grid.at(iE);
    double tmp = dsvdE_E(Kion, iE, E, mx, Fx2, q_grid, Fv, num_v);
    dsvdE_vec.push_back(tmp);
  }
  return dsvdE_vec;
}

//==============================================================================
//! Calculates total cross-section for given electron energy, Ei
//! in units of 10^-15 cm^2
double sigtot_E(const LinAlg::Matrix<double> Kion,
                const std::vector<double> &E_grid,
                const std::vector<double> &q_grid, double Ei) {
  double sigtot = 0.0;

  auto Emin = E_grid.at(0);
  auto Emax = E_grid.at(E_grid.size() - 1);
  auto qmin = q_grid.at(0);
  auto qmax = q_grid.at(q_grid.size() - 1);

  auto duE = std::log((Emax / Emin) / (double)(E_grid.size() - 1));
  auto duQ = std::log((qmax / qmin) / (double)(q_grid.size() - 1));

  for (std::size_t iE = 0; iE < E_grid.size(); iE++) {
    auto E = E_grid.at(iE);
    if (E >= Ei)
      continue;
    double qminus = std::sqrt(2.0 * Ei) - std::sqrt(2.0 * (Ei - E));
    double qplus = std::sqrt(2.0 * Ei) + std::sqrt(2.0 * (Ei - E));
    for (std::size_t iq = 0; iq < q_grid.size(); iq++) {
      double q = q_grid.at(iq);
      if (q < qminus || q > qplus)
        continue;
      double dEdq = E * q * duE * duQ;
      double FX2 = 1.0 / (q * q * q);
      sigtot += dEdq * FX2 * Kion(iE, iq);
    }
  }
  // double a02 = 0.0279841;
  double a02 = PhysConst::aB_cm * PhysConst::aB_cm;
  return (4.0 * M_PI * a02 / Ei) * sigtot;
}

//==============================================================================
//! Calculates total cross-section as function of incoming electron energy
//! in units of
std::vector<double> sigtot(const LinAlg::Matrix<double> Kion,
                           const std::vector<double> &E_grid,
                           const std::vector<double> &q_grid) {
  assert(Kion.rows() == E_grid.size());
  assert(Kion.cols() == q_grid.size());

  std::vector<double> sig_vec;
  sig_vec.reserve(E_grid.size());
  for (std::size_t iE = 0; iE < E_grid.size(); iE++) {
    auto Ei = E_grid.at(iE);
    double tmp = sigtot_E(Kion, E_grid, q_grid, Ei);
    sig_vec.push_back(tmp);
  }
  return sig_vec;
}

//==============================================================================
//==============================================================================
//==============================================================================
int main(int argc, char *argv[]) {
  const std::vector<std::string> a;

  const std::string in_filename = argc > 1 ? argv[1] : "";
  if (in_filename.empty()) {
    std::cout << "No file given (pass as command-line arg)\n";
    return 1;
  }

  // Other options for mass etc can be passed like this
  // nb: can put all the input options in a text file,
  // and run using ./xsection << text_file.txt
  std::string option_2 = argc > 2 ? argv[2] : "";
  std::cout << "Option_2 = " << option_2 << "\n";

  // Read in Kion
  const auto [E_grid, q_grid, Kion] = read_in_Kion(in_filename);

  // Check formats correct:
  assert(E_grid.size() == Kion.rows());
  assert(q_grid.size() == Kion.cols());

  // Define DM form-factor (here, heavy mediator)
  const auto Fx2 = [](double q) {
    (void)q; // suppress unused variable q warning
    return 1.0;
  };

  // Construct SHM object: Note: this uses km/s units
  StandardHaloModel fSHM{};

  // Define 'usable' f(f) function, in atomic units
  // Before: [f] = [1/(km/s)], since Integral(f)=1
  const auto Fv = [&](double v) {
    const auto v_kms = v * V_to_kms;
    const auto fv = fSHM.fv(v_kms);
    return fv / (1.0 / V_to_kms);
  };

  // manually set DM mass - should be input
  const auto mx = 10.0 / M_to_GeV;

  // calculate <ds.v>/dE
  const auto dsvde = dsvdE(Kion, mx, Fx2, E_grid, q_grid, Fv, 1000);

  // // For now, just output to screen
  // for (std::size_t i = 0; i < dsvde.size(); ++i) {
  //   std::cout << E_grid[i] << " " << dsvde[i] << "\n";
  // }

  // calculate sig(E) for e-e impact ion
  const auto sig_impact = sigtot(Kion, E_grid, q_grid);
  for (std::size_t i = 0; i < dsvde.size(); ++i) {
    std::cout << E_grid[i] << " " << sig_impact[i] << "\n";
  }

  // output to file
  std::ofstream ofile_sig;
  ofile_sig.open("impact_" + in_filename);
  std::cout << sig_impact.size() << std::endl;
  for (std::size_t i = 0; i < sig_impact.size(); i++) {
    ofile_sig << E_grid[i] * au_to_eV << " " << sig_impact[i] << "\n";
  }
  ofile_sig.close();

  // Can copy functions over from old dmXsection file
  // There, K was stored as 3D vector of vector...
  // Here, we drop the 'state' index (already summed), and use LinAlg::Matrix
  // instead

  // 1. Read in K(E,q) matrix
  // 2. Calculate ds/dE
  // 3. Calculate <ds.v>/dE
  // 4. Calculate sigma(E)
  // 5. Smear with Gausian
}
