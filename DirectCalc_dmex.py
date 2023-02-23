import numpy as np 
import matplotlib.pyplot as plt 
import SHM
from io import StringIO

# only heavy mediator rn
def F_chi_2(mu, q):
    return 1.0

# needs more testing, doesn't look quite right
# units?
def dsdE_Evmvmx(Kion, E, ide, qgrid, v, mv, mx):
    qsteps = len(qgrid)
    # esteps = len(Egrid)
    # assert esteps == len(Kion)
    assert qsteps == len(Kion[0])

    # EMAX, EMIN = max(Egrid), min(Egrid)
    QMAX, QMIN = max(qgrid), min(qgrid)
    duQ = np.log(QMAX / QMIN) / (qsteps - 1)

    m_e_MeV = 0.51099895000;
    M_au_to_GeV = m_e_MeV / 1000.0
    mx /= M_au_to_GeV
    mv /= M_au_to_GeV

    mu = mv * SHM.c_au
    
    arg = (mx * v)**2 - 2.0 * mx * E
    if (arg < 0): 
        return 0
    qminus = mx * v - np.sqrt(arg)
    qplus = mx * v + np.sqrt(arg)
    if qminus > QMAX or qplus < QMIN:
        return 0
    dsdE = 0.0
    for iq in range(qsteps):
        q = qgrid[iq]
        if q < qminus or q > qplus:
            continue
        qdq = q * duQ
        FX2 = F_chi_2(mu, q)
        dsdE += qdq * FX2 * Kion[ide, iq]
    return (0.5 / (v**2)) * dsdE

def dsvdE_Evmvmx(Kion, Egrid, qgrid, mv, mx, arr_fv, dv):
    esteps = len(Egrid)
    assert esteps == len(Kion)

    vsteps = len(arr_fv)
    
    dsvdE_array = []
    for ide in range(esteps):
        E = Egrid[ide]

        temp = 0.0
        v = dv
        for iv in range(vsteps):
            temp += dsdE_Evmvmx(Kion, E, ide, qgrid, v, mv, mx) * v * arr_fv[iv]
            v += dv
        dsvdE_array.append(temp * dv)
    return dsvdE_array

# filename = "testdata_mat.txt"
# filename = "K_Xe_v_6_hp_orth_mat.txt"
filename = "/home/uqcashle/ampsci/K_Xe_v_6_hp_orth_testdmex_q5_mat.txt"

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

print(len(K_array), len(K_array[0]))


fv_array, v_array, dv = SHM.fv_array()
print(len(fv_array))
dsvde_array = dsvdE_Evmvmx(K_array, E_array, q_array, 0.0, 10.0, fv_array, dv)

with open('dmex_py_20230223.txt','w') as f:
    for i in range(len(E_array)):
        print(str(E_array[i]*27.211) + " " + str(dsvde_array[i]), file=f)

plt.xscale('log')
plt.yscale('log')
plt.plot(E_array, dsvde_array)
plt.show()