#pragma once

namespace Physics {

namespace Constants {

//! fine-structure constant
constexpr double alpha = 1.0 / 137.035999084;

//! speed of light in a.u. (=1/alpha)
constexpr double c_au = 1.0 / alpha;

//! speed of light, in m/s [CODATA 2018, exact]
constexpr double c_SI = 299792458.0;

//! Bohr radius, in m. CODATA 2018: a_B = 0.529177210903(80)e-10 m
constexpr double aB_m = 0.529177210903e-10;

//! Bohr radius, in cm.
constexpr double aB_cm = aB_m * (1.0e2);

//! barn = 1.0e-28 m^2
constexpr double barn_m2 = 1.0e-28;
//! barn = 1.0e-24 cm^2
constexpr double barn_cm2 = barn_m2 * (1.0e4);

//! Barn, in atomic units
constexpr double barn_au = barn_m2 / (aB_m * aB_m);

//! unified atomic mass unit; (nuclear mass unit, Dalton).
//! CODATA 2014: 1822.888 486 192(53)
constexpr double u_NMU = 1822.888486192;

//! Proton mass, in atomic units (mp/me). CODATA 2018: 1836.152 673 43(11)
constexpr double mp_on_me = 1836.15267343;

//! electron mass, in SI (kg). CODATA 2018: m_e = 9.109 383 7015(28) e-31
constexpr double m_e_kg = 9.1093837015e-31;

//! Electron mass (MeV/c^2); 2020 value 0.51099895000(15):
constexpr double m_e_MeV = 0.51099895000;

//! Hartree (atomic energy unit = 2Ry) in eV.
//! CODATA 2018: E_H = 27.211 386 245 988(53) eV
constexpr double Hartree_eV = 27.211386245988;

} // namespace Constants

namespace Conversions {

//! Mass: Converts atomic units to (GeV/c^2)
constexpr double M_to_GeV = Constants::m_e_MeV / 1000.0;

//! Energy: Converts atomic units to eV
constexpr double Energy_au_to_eV = Constants::Hartree_eV;

//! Energy: Converts atomic units to keV
constexpr double E_to_keV = Constants::Hartree_eV / 1000.0;

//! Momentum: Converts atomic units to MeV:
//! [hbar*q] = [hbar/a0] = (m_e*c*alpha) = E_H/c*alpha
constexpr double Momentum_to_MeV =
    Constants::Hartree_eV / Constants::alpha / 1.0e6;

//! Momentum: MeV -> au
constexpr double Momentum_MeV_to_au = 1.0 / Momentum_to_MeV;

//! Velocity: Converts atomic units to km/s
constexpr double V_to_kms = (Constants::c_SI / Constants::c_au) / 1000.0;

//! Velocity: Converts atomic units to cm/s
constexpr double V_to_cms = (Constants::c_SI * Constants::alpha) * 100.0;

//! Velocity: Converts atomic units to cm/day
constexpr double V_to_cmday = V_to_cms * (24.0 * 60.0 * 60.0);

//! Converts ds/dE from atomic units to sigma_0/keV
constexpr double dsdE_to_cm2keV = 1.0 / E_to_keV; // au -> cm^2/keV

//! Converts <ds.v>/dE from atomic units to sigma_0*cm/keV/day
constexpr double dsvdE_to_cm3keVday = dsdE_to_cm2keV * V_to_cmday;
} // namespace Conversions

} // namespace Physics