import numpy as np
import scipy
from matplotlib import pyplot as plt
from numpy.polynomial import Polynomial

import fantastic.turbomachinery as tm

deg = 2 * np.pi / 360

f_0 = 0.5 * (0.7 + 1)
nu = 100 * 1e-6
Cf = 0.0025
d_m = 0.025
rho = 1.205
r_2 = 0.1
r_1 = 0

eff_c = 0.6
eff_t = 0.65

# dp0_comp = 6180
# dp0_turb = -17740
# dp0_vac = 11560
mdot = 0.0623
density = 1.204675848

W_VAC_MAX = 597.6087961307821


def main(rpm):
    w = rpm * 2 * np.pi / 60
    # DATA INPUT
    # inlet and outlet nodes
    inlet = tm.Node(101325, 293, 287)
    outlet = tm.Node(101325, 293, 287)

    # compressor blade rows (doesn't matter for now)
    com_rotor = tm.BladeRow(bi=0, bo=0, ri=0.025, ro=0.100, N=10)
    com_stator = tm.BladeRow(bi=0, bo=0, ri=0.110, ro=0.200, N=10)

    # turbine blade rows (doesn't matter for now)
    tur_rotor = tm.BladeRow(bi=0, bo=0, ri=0.020, ro=0.100, N=10)
    tur_stator = tm.BladeRow(bi=0, bo=0, ri=0.110, ro=0.200, N=10)

    # turbomachinery elements
    com = tm.Compressor(com_rotor, com_stator, h=0.012)
    tur = tm.Turbine(tur_rotor, tur_stator, h=0.012)
    vac = tm.Vacuum()

    # rig setup
    rig = tm.Rig(inlet, [com, tur, vac], outlet)
    rig.add_to_shaft(idxs=[0, 1])

    # COMPRESSOR DESIGN [notation: V_{(r)adial or (t)angential}_{(r)otor or (s)tator}{(i)nlet or (o)utlet}]
    dh0 = tm.dp0(w) / inlet.density / eff_c
    V_r_ri = mdot / 2 / np.pi / inlet.density / com.rotor.ri / com.h
    V_r_ro = mdot / 2 / np.pi / inlet.density / com.rotor.ro / com.h

    # solve quadratic for bo from omega

    def f(bo):
        return Polynomial(
            [-dh0, V_r_ro * com.rotor.ro * np.tan(bo), com.rotor.sigma * com.rotor.ro ** 2]
        ).roots()[1] - w

    # set parameters
    com.rotor.bo = scipy.optimize.root_scalar(f, bracket=(0.001 - np.pi / 2, np.pi / 2 - 0.001)).root
    rig.shaft.omega = f(com.rotor.bo) + w

    # find rotor vels
    V_t_ro = dh0 / w / com.rotor.ro
    U_ri = w * com.rotor.ri
    U_ro = w * com.rotor.ro
    # W is relative velocity
    W_ri = np.sqrt(V_r_ri ** 2 + U_ri ** 2)
    W_ro = np.sqrt(V_r_ro ** 2 + (V_t_ro - U_ro) ** 2)

    # find bi: beta_in = (alpha_in-5deg) as per handout
    # com.rotor.bi = np.arctan(U_ri / V_r_ri) - 5 * deg
    com.rotor.bi = np.arctan(U_ri / V_r_ri) - 5 * deg

    # rule of thumb
    print(f'W_ro / W_ri > 1 / 3: {W_ro / W_ri > 1 / 3}')
    print(f'backsweep: {V_t_ro < U_ro}, {com.rotor.bo>0}')

    # results are ok to show!
    print(f'beta_1: {com.rotor.bi / deg:.2f}deg')
    print(f'beta_2: {com.rotor.bo / deg:.2f}deg')
    print(f'omega: {rig.shaft.omega * 30 / np.pi:.0f}rpm')

    # now for blade numbers
    # first get blade profile for constant deltaP across blade
    b = tm.rotor_shape(com.rotor)
    # vels
    V_r_r = mdot / 2 / np.pi / inlet.density / r / com.h
    V_t_r = w * r + V_r_r * np.tan(b)
    W_r = V_r_r / np.cos(b)

    drVt_dr = scipy.interpolate.InterpolatedUnivariateSpline(r, r * V_t_r).derivative(1)
    Nb_min = mdot * drVt_dr(r) / (2 * inlet.density * W_r ** 2 * r * com.h)

    Nb = int(np.ceil(1.25 * np.max(Nb_min)))
    com.rotor.N = Nb
    print(f'num. compressor blades: {Nb}')

    for i in range(Nb):
        plt.polar(np.tan(b) * np.log(r) + c1 + i * 2 * np.pi / Nb, r, 'black')
    plt.yticks([])
    plt.ylim(0, None)
    plt.show()


if __name__ == "__main__":
    main(8930)
