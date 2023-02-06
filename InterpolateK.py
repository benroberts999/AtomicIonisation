#!/usr/bin/env python3
import numpy as np
from matplotlib import pyplot as plt
from io import StringIO  # StringIO behaves like a file object
from scipy.interpolate import LinearNDInterpolator, CloughTocher2DInterpolator
from scipy.integrate import dblquad, quad


def log_space_interp(X, Y, Z):
    """Linear interpolation in the log scale. X is 1D array of size SX, Y is 1D array of size SY, and Z is a 2D array of size SX*SY. Z_ij = Z(X_i, Y_j)"""
    xx, yy = np.meshgrid(X, Y)
    interp_log = LinearNDInterpolator(list(zip(np.log(xx.flatten('F')), np.log(
        yy.flatten('F')))), np.log(Z.flatten()))

    def interped_func(x, y):
        return np.exp(interp_log(np.log(x), np.log(y)))
    return interped_func


# filename = "testdata_mat.txt"
filename = "K_Xe_v_6_hp_orth_mat.txt"

text_file = open(filename).read()
# print(text_file)
out = text_file.split('\n\n')

# check the data makes sense:
assert len(out) == 4

E_array = np.loadtxt(StringIO(out[1]), comments='#')
q_array = np.loadtxt(StringIO(out[2]), comments='#')
K_array = np.loadtxt(StringIO(out[3]), comments='#')


Emin, Emax = min(E_array), max(E_array)
Qmin, Qmax = min(q_array), max(q_array)

# check the data makes sense:
assert len(K_array) == len(E_array)
assert len(K_array[0]) == len(q_array)

print("Energy deposition steps: ", len(E_array))
print("Momentum transfer steps: ", len(q_array))

# Create K(E,q) function by interpolation:
Kion = log_space_interp(E_array, q_array, K_array)

# Check to see if interpolation working:
# Just for testing: use a slightly different q-grid for interp'd function:
q0, qmax = min(q_array), max(q_array)
q_array_2 = np.logspace(np.log10(q0), np.log(qmax), 300)

plt.xscale('log')
plt.yscale('log')
for iE, E in enumerate(E_array):
    plt.plot(q_array, K_array[iE], label="E="+str(E))
    if (iE+1 < len(K_array)):
        E2 = np.exp(0.5*(np.log(E) + np.log(E_array[iE+1])))
        plt.plot(q_array, Kion(E2, q_array), ':', label="E=" +
                 str(E))
leg = plt.legend(loc='best')
plt.show()


# Electron-impact cross-section:

# sigma= (4*pi/E) Int_0^E dE Int_q0^qmax K(E,q)/q^3

def integral_q(E_i, delta_E, Kion):
    factor = 4*np.pi / E_i
    if (delta_E > E_i):
        return 0.0
    qmin = np.sqrt(2*E_i) - np.sqrt(2*(E_i - delta_E))
    qmax = np.sqrt(2*E_i) + np.sqrt(2*(E_i - delta_E))

    def f(q):
        return Kion(delta_E, q)/(q**3)
    return quad(f, qmin, qmax)[0]


def sigma_impact2(E_i, Kion):
    """Units 10^-15 cm^2"""
    factor = 4*np.pi / E_i
    a02 = 2.79841

    def f(delta_E):
        return integral_q(E_i, delta_E, Kion)
    return factor * a02 * quad(f, 0.0, E_i)[0]


def Esigma_impact(E_0, E_i, Kion):
    """Units 10^-15 cm^2"""
    factor = 4*np.pi
    a02 = 0.0279841

    def h(q, dE):
        if dE < Emin or dE > Emax or q < Qmin or q > Qmax:
            return 0.0
        if dE > E_i:
            return 0.0
        return Kion(dE, q)/(q**3)

    def qmin(dE): return np.sqrt(2*E_i) - np.sqrt(2*(E_i - dE))
    def qmax(dE): return np.sqrt(2*E_i) + np.sqrt(2*(E_i - dE))
    res = dblquad(h, E_0, E_i, qmin, qmax, epsabs=1.0e-6, epsrel=1.0e-3)[0]
    return factor * a02 * res

# XXX this integral is too unstable!


# es = np.linspace(10/27.211, 1000/27.211, 100)
es = np.logspace(np.log10(10/27.211), np.log(1000/27.211), 25)
e_prev = 0.0
Esigma_prev = 0.0
y = []
for e in es:
    Esig = Esigma_prev + Esigma_impact(e_prev, e, Kion)
    e_prev = e
    Esigma_prev = Esig
    y.append(Esig/e)
    # print(e*27.211, s/e)

plt.xscale('log')
plt.plot(es*27.211, y)
plt.show()
