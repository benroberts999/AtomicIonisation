#pragma once
#include <cmath>
#include <iostream>

//! For Standard Halo Model, and astrophysics constants
namespace Astro {

//! All Constants given in units of km/s (except cos_beta, which is
//! dimensionless, and rhoDM_GeVcm3)
namespace Constants {

//! DM energy density (in GeV/cm^3)
constexpr double rhoDM_GeVcm3 = 0.4; // GeV/cm^3

//! Speed of light, in km/s
constexpr double c_kms = 299792.458;

//! galactic escape velocity [1804.01231]
constexpr double v_esc = 550.0;

//! Most probable speed (& disk rotation speed). RevModPhys.85.1561
constexpr double v0 = 235.0;

//! Speed of sun, in galactic frame |{11, 12+v0, 7}|
constexpr double v_sun = 247.3;

//! [RevModPhys.85.1561]
constexpr double v_earth_orb = 29.8;

//! 1804.01231 (Earth inclination to sun dir)
constexpr double cos_beta_orb = 0.49;

//! earth rotation speed (approx) @ equator
constexpr double v_earth_rot_eq = 0.47;

constexpr double v_max = v_esc + v_sun + v_earth_orb + v_earth_rot_eq;

} // namespace Constants

//==============================================================================
// Approximate, non-normalised speed distribution
inline double f_speed_approx_nn(double v, double vobs, double v0, double vesc) {
  if (v <= 0.0 || v > vesc + vobs)
    return 0.0;
  const auto v2 = v * v;
  const auto v02 = v0 * v0;
  const auto vobs2 = vobs * vobs;
  const auto arg_minus = -(v2 + vobs2 - 2.0 * v * vobs) / v02;
  const auto arg_plus = -(v2 + vobs2 + 2.0 * v * vobs) / v02;
  const auto arg_vesc = -vesc * vesc / v02;
  return v < (vesc - vobs) ? v * (std::exp(arg_minus) - std::exp(arg_plus))
                           : v * (std::exp(arg_minus) - std::exp(arg_vesc));
}

//==============================================================================
inline double heaviside(double x) { return (x >= 0.0) ? 1.0 : 0.0; }

//==============================================================================
// Non-normalised speed distribution, maxwellian
inline double f_speed_maxwell_nn(double v, double vobs, double v0,
                                 double vesc) {
  if (v <= 0.0 || v > vesc + vobs)
    return 0.0;

  const auto v02 = v0 * v0;
  const auto vesc2 = vesc * vesc;
  const auto v_plus_vobs_2 = std::pow(v + vobs, 2);
  const auto arg1 = (vesc2 - v_plus_vobs_2) / v02;
  const auto arg2 = 4.0 * v * vobs / v02;

  if (v < vesc - vobs)
    return v * ((std::exp(arg1 + arg2) - std::exp(arg1)) - arg2);

  const auto hst_minus = heaviside(std::abs(v - vobs) - vesc);
  const auto hst_plus = heaviside(v + vobs - vesc);
  const auto a = (v - vesc + vobs) * (v + vesc + vobs) * hst_plus;
  const auto b = (v - vesc - vobs) * (v + vesc - vobs) * hst_minus;
  const auto c = vobs == 0.0 ? 0.0 : v02 * (hst_minus - hst_plus);
  const auto abc_v02 = vobs == 0.0 ? 0.0 : (a - b + c) / v02;
  const auto d = (hst_plus - 1.0) - std::exp(arg2) * (hst_minus - 1.0);
  return v * (std::exp(arg1) * d + abc_v02 - arg2);
}

//==============================================================================
//! Calculates approximate lab speed |vearth + vsun|. cos_phi_orbit is phase of
//! yearly orbit, [-1,1]. Does not account for daily modulation.
inline double abs_v_obs(double cos_phi_orbit, double vsun = Constants::v_sun,
                        double ve_orb = Constants::v_earth_orb,
                        double cos_beta = Constants::cos_beta_orb) {
  return vsun + cos_phi_orbit * ve_orb * cos_beta;
}

//==============================================================================
//! Simplistic implementation of standard halo model.
class StandardHaloModel {
private:
  double (*m_f_v_nn)(double, double, double, double) = f_speed_approx_nn;
  double m_vobs;
  double m_v0;
  double m_vesc;
  double m_norm;

  // calculates norm constant
  double normalise() const {
    const double vmax = m_vobs + m_vesc;
    constexpr int num_steps = 100 + 1;
    static_assert(num_steps % 2 != 0, "num_steps must be odd (Simpsons Rule)");
    double dv = vmax / (num_steps - 1);
    double integral = 0.0 + m_f_v_nn(vmax, m_vobs, m_v0, m_vesc);
    double v = dv;
    for (int i = 1; i < num_steps - 1; ++i) {
      const auto w = (i % 2 == 0) ? 2.0 : 4.0;
      integral += w * m_f_v_nn(v, m_vobs, m_v0, m_vesc);
      v += dv;
    }
    integral *= (dv / 3.0);
    return 1.0 / integral;
  }

public:
  StandardHaloModel(double vobs = Constants::v_sun, double v0 = Constants::v0,
                    double vesc = Constants::v_esc, bool use_approx = true)
      : m_f_v_nn(use_approx ? f_speed_approx_nn : f_speed_maxwell_nn),
        m_vobs(vobs), m_v0(v0), m_vesc(vesc), m_norm(normalise()) {}

  //! Lab-frame DM speed (|v|) distribution.
  //! Given v [in km/s], returns f(v) [in (km/s)^-1].
  //! f(v) is normalised as: Int[f(v) {v,0,infty}]=1.
  double f(double v) const {
    return m_norm * m_f_v_nn(v, m_vobs, m_v0, m_vesc);
  }

  //! maximum DM speed in local frame
  double vmax() const { return m_vesc + m_vobs; }
  double vobs() const { return m_vobs; }
  double v0() const { return m_v0; }
  double vesc() const { return m_vesc; }
};

} // namespace Astro