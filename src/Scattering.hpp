#pragma once
#include "Matrix.hpp"
#include "PhysicsConstants.hpp"
#include "StandardHaloModel.hpp"
#include <functional>
#include <string>
#include <tuple>
#include <vector>

//! Formulas for calculating ionisation cross-sections and rates
namespace Scattering {

//! Reads E grid, q grid, and Kion(E,q) from '_mat' formatted text file. Note:
//! assumes it was in ATOMIC units.
std::tuple<std::vector<double>, std::vector<double>, Matrix<double>>
read_in_Kion(const std::string &in_filename);

//! Calculates dsigma/dE for given energy deposition, E, and velocity, v,
//! in units of sigma-bar_e (v and E are in atomic units)
double dsdE_E(const Matrix<double> &Kion, std::size_t iE, double E, double v,
              double mx, const std::function<double(double)> &Fx2,
              const std::vector<double> &q_grid);

//! Calculates <dsigma.v>/dE for given energy deposition, E, in units of
//! sigma-bar_e (E in atomic units)
double dsvdE_E(const Matrix<double> &Kion, std::size_t iE, double E, double mx,
               const std::function<double(double)> &Fx2,
               const std::vector<double> &q_grid,
               const std::function<double(double)> &Fv, int num_v);

//! Calculates <ds.v>/dE as function of energy deposition, E,
//! in units of sigma-bar_e (E in atomic units)
std::vector<double> dsvdE(const Matrix<double> Kion, double mx,
                          const std::function<double(double)> &Fx2,
                          const std::vector<double> &E_grid,
                          const std::vector<double> &q_grid,
                          const std::function<double(double)> &Fv, int num_v);

//! Simple Gaussian, g(x) with standard deviation sigma
double gaussian(double sigma, double x);

//! Convolution across energy grid: dSdE(e) = Int[ dRdE(e') f(e,e') de']
std::vector<double>
convolvedRate(const std::vector<double> &dRdE,
              const std::function<double(double, double)> &f,
              const std::vector<double> &E_grid);

//! Returns binned counts, (1/(Eb-Ea))*Int_Ea^Eb dSdE dE. Ea and Eb must be
//! given in atomic units. Result will have the same units as dSdE.
double binnedCount(const std::vector<double> &dSdE, double Ea, double Eb,
                   const std::vector<double> &E_grid);

//! Calculates total electron-impact ionisation cross-section for given
//! incoming electron energy, Ei, in units of a0^2 (i.e., atomic units)
double sigtot_E(const Matrix<double> &Kion, const std::vector<double> &E_grid,
                const std::vector<double> &q_grid, double Ei);

//! Calculates total electron-impact ionisation cross-section as function of
//! incoming electron energy, in units of a0^2 (i.e., atomic units)
std::vector<double> sigtot(const Matrix<double> &Kion,
                           const std::vector<double> &E_grid,
                           const std::vector<double> &q_grid);

} // namespace Scattering