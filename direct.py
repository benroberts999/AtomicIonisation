#!/usr/bin/env python3
import numpy as np
from matplotlib import pyplot as plt
from io import StringIO  # StringIO behaves like a file object


def sigmaE(Kion, Egrid, qgrid, Ei):
    """Calculates electron-impact ionisation x-section by direct integration. 
    Assumes E and q are logarithmically-spaced grids. 
    Inputs must all be in atomic units; output is in units of 10^-15 cm^2"""
    qsteps = len(qgrid)
    esteps = len(Egrid)
    assert esteps == len(Kion)
    assert qsteps == len(Kion[0])

    EMAX, EMIN = max(Egrid), min(Egrid)
    QMAX, QMIN = max(qgrid), min(qgrid)

    # assume Logarithimc grid
    # dr = (dr/du)*du
    # (dr/du) = r
    # du = log(max / min) / (N - 1)
    duE = np.log(EMAX / EMIN) / (esteps - 1)
    duQ = np.log(QMAX / QMIN) / (qsteps - 1)

    sig = 0.0
    for ide in range(esteps):
        Et = Egrid[ide]
        if Et >= Ei:
            continue
        qminus = np.sqrt(2 * Ei) - np.sqrt(2 * (Ei - Et))
        qplus = np.sqrt(Ei) + np.sqrt(Ei - Et)
        for iq in range(qsteps):
            q = qgrid[iq]
            if q < qminus or q > qplus:
                continue
            dEdq = Et*q*duE*duQ
            FX2 = 1 / (q * q * q)
            sig += dEdq * FX2 * Kion[ide, iq]

    # Bohr radius, a0^2, in units of 10^-15 cm^2
    a02 = 0.0279841
    return (4 * np.pi / Ei) * a02 * sig


# filename = "testdata_mat.txt"
filename = "K_Xe_v_6_hp_orth_mat.txt"

text_file = open(filename).read()
out = text_file.split('\n\n')

# check the data makes sense:
assert len(out) == 4

E_array = np.loadtxt(StringIO(out[1]), comments='#')
q_array = np.loadtxt(StringIO(out[2]), comments='#')
K_array = np.loadtxt(StringIO(out[3]), comments='#')

# check the data makes sense:
assert len(K_array) == len(E_array)
assert len(K_array[0]) == len(q_array)

print("Energy deposition steps: ", len(E_array))
print("Momentum transfer steps: ", len(q_array))

es = np.logspace(np.log10(10/27.211), np.log10(1000/27.211), 150)
e_prev = 0.0
Esigma_prev = 0.0
# y = []
y2 = []
y3 = []
for e in es:
    sig2 = sigmaE(K_array, E_array, q_array, e)
    y3.append(sig2)

plt.xscale('log')
plt.plot(es*27.211, y3, label="direct")
leg = plt.legend(loc='best')
plt.show()
