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
const double E_to_keV = PhysConst::Hartree_eV / 1000.;

const double V_to_cms =
    (PhysConst::c_SI * PhysConst::alpha) * 100.;     // au -> cm/s
const double V_to_cmday = V_to_cms * (24 * 60 * 60); // cm/s -> cm/day

// DM energy density (in GeV/cm^3)
const double rhoDM_GeVcm3 = 0.4; // GeV/cm^3

// Default value used for bar-sigma_e + conversions from a.u.
const double sbe_1e37_cm2 = 1.e-37;
const double dsdE_to_cm2keV = sbe_1e37_cm2 / E_to_keV; // au -> cm^2/keV
const double dsvdE_to_cm3keVday = dsdE_to_cm2keV * V_to_cmday;

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

// //==============================================================================
// //! Quick exponentional, good to parts in 1.e-5
// double quickExp(double x) {
//   if (x < 0)
//     return 1.0 / quickExp(-x);
//   if (x < 0.05)
//     return 1.0 + x + 0.5 * x * x;
//   if (x < 0.25)
//     return 1.0 + x + 0.5 * x * x + 0.16666667 * x * x * x +
//            0.041666667 * x * x * x * x;
//   return std::exp(x);
// }

//==============================================================================
//! Simple Gaussian, g(x) with standard dev of sigma
double gaussian(double sigma, double x) {
  double a = 0.398942 / sigma;
  double y = (x / sigma) * (x / sigma);
  // return a * quickExp(-0.5 * y);
  return a * std::exp(-0.5 * y);
}

//==============================================================================
//! Convolution across energy grid
std::vector<double> convolvedRate(const std::vector<double> &dRdE,
                                  const std::vector<double> &E_grid,
                                  double units) {
  std::vector<double> conv;
  conv.reserve(E_grid.size());

  const auto Emin = *std::min_element(E_grid.begin(), E_grid.end());
  const auto Emax = *std::max_element(E_grid.begin(), E_grid.end());
  const auto duE = std::log(Emax / Emin) / double(E_grid.size() - 1);

  double alpha = 0.310;
  double beta = 0.0037;

  std::vector<std::vector<double>> gausVec(E_grid.size());

  for (std::size_t i = 0; i < E_grid.size(); i++) {
    double Eobs = E_grid[i];
    double sigma =
        (alpha * std::sqrt(Eobs * E_to_keV) + beta * (Eobs * E_to_keV)) /
        E_to_keV;
    for (std::size_t j = 0; j < E_grid.size(); j++) {
      double Eer = E_grid[j];
      double g = gaussian(sigma, Eobs - Eer);

      if (Eer < E_grid[0] || Eobs < E_grid[0])
        g = 0.0;

      gausVec[i].push_back(g);
    }
  }

  for (std::size_t iE = 0; iE < E_grid.size(); iE++) {
    auto E = E_grid.at(iE);

    // double sigma =
    //     (alpha * std::sqrt(E * E_to_keV) + beta * (E * E_to_keV)) / E_to_keV;

    double tmp = 0.0;
    for (std::size_t iEp = 0; iEp < E_grid.size(); iEp++) {
      auto Ep = E_grid.at(iEp);
      double dEdt = E;
      // double g = gaussian(sigma, Ep - E);
      tmp += gausVec[iE][iEp] * dRdE[iEp] * dEdt;
    }
    tmp *= duE * units;
    conv.push_back(tmp);
  }
  return conv;
}

//==============================================================================
//! Calculates total cross-section for given electron energy, Ei
//! in units of 10^-15 cm^2
double sigtot_E(const LinAlg::Matrix<double> Kion,
                const std::vector<double> &E_grid,
                const std::vector<double> &q_grid, double Ei) {
  double sigtot = 0.0;

  const auto Emin = *std::min_element(E_grid.begin(), E_grid.end());
  const auto Emax = *std::max_element(E_grid.begin(), E_grid.end());
  const auto qmin = *std::min_element(q_grid.begin(), q_grid.end());
  const auto qmax = *std::max_element(q_grid.begin(), q_grid.end());

  const auto duE = std::log(Emax / Emin) / double(E_grid.size() - 1);
  const auto duQ = std::log(qmax / qmin) / double(q_grid.size() - 1);

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
      double dEdq = E * q;
      double FX2 = 1.0 / (q * q * q);
      sigtot += dEdq * FX2 * Kion(iE, iq);
    }
  }
  sigtot *= duE * duQ;
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

  double mn_xe = 131. * (PhysConst::u_NMU * PhysConst::m_e_kg);
  // Convert units to counts/day/kg/keV
  double rho_on_mxc2 = rhoDM_GeVcm3 / (mx * M_to_GeV);
  double rate_units = dsvdE_to_cm3keVday * rho_on_mxc2 / mn_xe;

  // // For now, just output to screen
  // for (std::size_t i = 0; i < dsvde.size(); ++i) {
  //   std::cout << E_grid[i] << " " << dsvde[i] << "\n";
  // }
  std::ofstream ofile_dsvde;
  ofile_dsvde.open("dsvde_" + in_filename);
  for (std::size_t i = 0; i < dsvde.size(); i++) {
    ofile_dsvde << E_grid[i] * E_to_keV << " "
                << dsvde[i] * rate_units * 1000.0 * 365.0 << "\n";
  }

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

  // // Total atomic mass for xenon in kg
  // double mn_xe = 131. * (PhysConst::u_NMU * PhysConst::m_e_kg);

  // Efficiency as function of deposited energy of Xe1T detector
  auto effic = [](double energy) {
    // if (energy < 0.5 / E_to_keV)
    //   return 0.0;
    // else
    return (0.876 - 7.39e-04 * energy * E_to_keV) /
           std::pow(
               (1.0 + 0.104 * std::exp(-(energy * E_to_keV - 1.98) / 0.360)),
               2.03);
  };

  // // Convert units to counts/day/kg/keV
  // double rho_on_mxc2 = rhoDM_GeVcm3 / (mx * M_to_GeV);
  // double rate_units = dsvdE_to_cm3keVday * rho_on_mxc2 / mn_xe;

  std::vector<double> dSdE_E;
  dSdE_E.resize(E_grid.size());

  dSdE_E = convolvedRate(dsvde, E_grid, rate_units);
  // for (std::size_t iE = 0; iE < E_grid.size(); iE++) {
  //   double E = E_grid[iE];
  //   dSdE_E[iE] *= effic(E);
  // }

  std::ofstream ofile_dsde;
  ofile_dsde.open("obsratenewconv_noeffic_" + in_filename);
  for (std::size_t i = 0; i < dSdE_E.size(); i++) {
    ofile_dsde << E_grid[i] * E_to_keV << " " << dSdE_E[i] * 1000.0 * 365.0
               << "\n";
  }
  ofile_dsde.close();
}
