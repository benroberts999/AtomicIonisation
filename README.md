# xsection

[![arXiv][arXiv-badge]][arXiv-url]

This repository is the companion code for the above paper.
It provides tables of high-accuracy atomic ionisation factors (matrix elements),
which are required to calculate atomic ionisation rates, including from dark 
matter electron scattering.
It also provides example programs which read in the tables, and calculate 
example DM-electron-induced ionisation rates, accounting for the low-energy
detector response.

## Contents

1. [DM-electron scattering](#1-dm-electron-scattering-atomic-ionisation-factor)
2. [Tables of atomic ionisation factors](#tables-of-atomic-ionisation-factors)
3. [DM-induced ionisation rate: example](#3-dm-induced-ionisation-rate-example-calculation)
4. [Electron-impact ionisation](#4-electron-impact-ionisation-eimpact)

----

## 1. DM-electron scattering: atomic ionisation factor

We define the atomic ionisation factor:

$$ \tag{1}
K(E,q) = \sum_{n\kappa m}^{\rm bound}\sum_{\kappa' m'}^{\rm excited}
          E_H |\langle \varepsilon\kappa'm'|\Gamma e^{i\mathbf{q}\cdot\mathbf{r}}|n\kappa m\rangle|^2
$$

Where:
  - $E$ is the _energy deposition_, $E = \varepsilon - \varepsilon_{n\kappa}$
  - $\varepsilon>0$ is the final energy of the ionised (outgoing) electron
  - $\varepsilon_{n\kappa}<0$ is the initial energy of the bound atomic electron
    - (Though we use relativistic wavefunctions, the rest energy is subtracted from all energies)
  - $n$ and $\kappa$ are the principal and relativistic angular quantum numbers, with $m=j_z$
  - $q$ is the momentum transfer (technically $\hbar q$ is momentum, but we call $q$ momentum here)
  - $E_H \approx 27.211$ eV is the Hartree energy unit. It is included to make $K$ dimensionless
    - (All calculations are performed in atomic units, where $E_H=1$)

Note that the final (ionised) electron states, $|\epsilon\kappa'm'\rangle$,
are energy eigenstates, _not_ momentum eigenstates.
The destinction is important, and leads to a large Sommerfield-like enhancement, 
neglection of which can lead to results incorrect by many orders-of-magnitude.
Therefore, we normalise these states on the energy scale:

$$ \tag{2}
\int_{\varepsilon-\delta}^{\varepsilon+\delta}
\langle \varepsilon'\kappa'm'|\varepsilon\kappa m\rangle\ {\rm d}\varepsilon'
 = \delta_{\kappa \kappa'}\delta_{mm'}
$$

We consider only cases where 
  - $\Gamma=1$ (vector electron-coupling)
  - $\Gamma=\gamma^0$ (scalar)
  - $\Gamma=\gamma^5$ (pseudovector)
  - $\Gamma=\gamma^0\gamma^5$ (pseudoscalar)

----

## 2. Tables of atomic ionisation factors

In the directory [tables/](./tables/), we provide tables of atomic ionisation factors K.

- [K_Xe_v_6_hp_orth_mat.txt](./tables/K_Xe_v_6_hp_orth_mat.txt)
  - $K$ ionisation factor for Xe, vector couling, including up to $L=6$, accounting for hole-particle interaction, with explicitely orthogonalised continuum orbitals

These tables have been generated using [_ampsci_](https://ampsci.dev/), a C++ program for high-precision atomic structure calculations.
Example input/output files from the ampsci runs that generated these tables are also presented in the
[tables/ampsci](./tables/ampsci) subdirectory.

The tables are in the form of 3 blocks:
 - First block contains the energy deposition values, $E$, at which $K$ was calculated.
    - The $E$ values are distributed on a _logarithmic_ grid, and are in _atomic units_
    - The atomic unit of energy is the Hartree: $E_H\approx 27.211\ {\rm eV}$
 - Second block contains the momentum transfer values, $q$, at which $K$ was calculated.
    - The $q$ values are distributed on a _logarithmic_ grid, and are in _atomic units_
    - The atomic unit of momentum is $1/a_0\approx 3.7289\ {\rm keV}$
 - Third block contains the $K(E,q)$ factor, in matrix form.
    - $K$ is dimensionless
    - Each new row is for a new $E$ value, each column is for a new $q$ value

We have provided example programs that read in files of this form:

 * [dmex.cpp](./dmex.cpp) - a C++ version, which also calculates DM-induced ionisation event rates
 * [eimpact.cpp](./eimpact.cpp) - a C++ program, which calculates electron-impact ionisation cross-section
 * [InterpolateK.py](InterpolateK.py) - a short python script that reads in K, created an interpolating function for $K(E,q)$, and plots $K$ as an example

----

## 3. DM-induced ionisation rate: example calculation

In the example calculation ([dmex.cpp](./dmex.cpp)), we calculate atomic 
ionisation event rates due to dark-matter electron scattering.

For the case of a vector coupling (treating the DM particle non-relativistically),
the cross-section may be expressed:

$$ \tag{3}
\frac{d\langle\sigma v\rangle}{d E} = \frac{\bar\sigma_e c}{2m_e c^2}
  \int d v \frac{f(v)}{v/c}
  \int_{q_-}^{q_+} \, a_0^2 \, q {\rm d}q |F_\chi(q)|^2 K(E,q),
$$

where
 - $\bar\sigma_e$ is the free electron cross section at fixed $q=a_0^{-1}$
 - $a_0$ is the Bohr radius ($a_0$ and $c$ are introduced to make the dimensions explicit)
 - $q_{\pm}$ is the kinematically-allowed range for $q$ 
 (depends on the DM mass $m_\chi$, $E$, and $v$)
 - $f(v)$ is the normalised lab-frame DM speed distribution
 - $F_\chi$ is the DM form factor.
    - For a heavy mediator, $F_\chi=1$. For a light mediator, $F_\chi=q^{-2}$

The underlying event rate

$$
\frac{dR}{dE} = \frac{n_T \rho_\chi}{m_\chi c^2}\frac{d\langle\sigma v\rangle}{d E}
$$

where
  - $n_T$ is the number of target atoms
  - $m_\chi$ is the DM mass
  - $\rho_\chi\approx0.4\ {\rm GeV}\ {\rm cm}^{-3}$ is the local DM density

To model the detector response (energy resolution, efficiency, and acceptance), 
we follow the XENON Collaboration, ####.
The observable event rate, $S$, is written as a convolution of the underlying ionisation rate, $R$.

$$
\frac{dS}{dE} = \epsilon(E)\int \frac{dR}{dE}(E')g_{\sigma_E}(E-E')\ {\rm d}E'
$$

We use a simple model for the detector efficiency/acceptance,

$$
\epsilon(E) = \epsilon_0 \left[1 + \gamma e^{-\delta (E/{\rm keV} - 2.0)}\right]^{-1}.
$$

We take $\gamma=0.19$, $\delta=3.3$, $\epsilon_0=0.87$, which come from a fit to Fig. ## of ####

For the detector resolution, we also follow ####, and use a Gaussian, $g_{\sigma_E}$,
with energy-dependent standard deviation:

$$
\sigma_E / {\rm keV} = \alpha \sqrt{E/{\rm keV}} + \beta(E/{\rm keV}),
$$

where $\alpha=0.310$, $\beta=0.0037$ ####.

We stress, that for very low-energy electron recoil events, the observable event
rate is extremely sensitive to the low-energy detector resolution. This 
simplistic Gaussian model is almost certainly not appropriate for modelling of
such events. Care must be taken to interpret the results.

The values $\alpha$, $\beta$, $\gamma$, $\delta$, $\epsilon_0$ can all be modified in `dmex`, 
demonstrating their importance.

For the velocity distribution, the simple pseudo-Maxwellian Standard Halo Model is assumed.
Standard values for parameters of the model are chosen, these can be easily updated in the code.

----

### 3.1 example: running `dmex`

The program can be compiled using the provided Makefile.

Takes input directly from command line.
The first argument is name of the file containing Kion.
(File generated by ampsci in 'mat' format, assumed to be in atomic units).

The following input options are all optional, but must be given in order

  - Mediator mass, in MeV (may be 'inf')      [default: inf]
  - DM mass in GeV                            [default: 1.0]
  - Nuclear mass number (A)                   [default: 131]
  - cos(phi_orbit), for annual modulation     [default: 0.0]
  - alpha (Gaussian resolution parameter)     [default: 0.31]
  - beta (Gaussian resolution parameter)      [default: 0.0037]
  - gamma (detector efficiency parameter)     [default: 0.19]
  - delta (detector efficiency parameter)     [default: 3.3]
  - epsilon (detector efficiency parameter)   [default: 0.87]

e.g.,
  1. `./dmex K_Xe_v_6_hp_orth_mat.txt > rates.txt`
  2. `./dmex K_Xe_v_6_hp_orth_mat.txt inf 0.5 > rates.txt`
  3. `./dmex K_Xe_v_6_hp_orth_mat.txt inf 1 131 0 0.31 0.0037 0.019 3.3 0.87`

(1) Will calculate the rateswith all default inputs, using the 
    Kion form-factor stored in file: 'K_Xe_v_6_hp_orth_mat.txt'.
    The results are written to 'rates.txt', ready for plotting.

(2) Will calculate the rates for a heavy mediator, DM mass of 0.5 GeV.
    All other inputs will have their default values.
    The results are written to 'rates.txt', ready for plotting.

(3) Will calculate the rates, with all input options are set explicitly.
    Results written to screen.

----

## 4. Electron-impact ionisation: `eimpact`

The program [eimpact.cpp](./eimpact.cpp) calculated electron impact ionisation cross-sections, 
using the same atomic ionisation factors.

Electron-impact ionisation is a very similar problem to DM-induced ionisation.
There exist high-quality experimental data at high $sim\ {\rm keV}$ impact energies, 
which provide an ideal test case for testing the atomic ionisation factors.

$$
\sigma(E) = \frac{4\pi}{E}\int_0^E\int_{q_-}^{q_+}\ {\rm d}q \frac{K(E,q)}{q^3}\ {\rm d}E
$$

where the $K$ factor is exactly the same as in the (vector-coupled) DM case.

In particular, comparison to experiment for Xe shows our calculations are highly accurate, 
and demonstrates the importance of accurate modelling of the atomic wavefunctions.

Usage:

Takes input directly from command line.
First (and only) argument is file containing Kion (the file generated by ampsci
in 'mat' format, assumed to be in atomic units).

e.g.,
  `./eimpact K_Xe_v_6_hp_orth_mat.txt > sigma_e.txt`

sigma_e.txt will contain the cross-section, ready for plotting.

----

[arXiv-badge]: https://img.shields.io/badge/arXiv-xxxx.yyyyy-b31b1b.svg?style=flat
[arXiv-url]: https://arxiv.org/abs/xxxx.yyyyy