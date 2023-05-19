import numpy as np
import scipy
from matplotlib import pyplot as plt
from numpy.polynomial import Polynomial

import fantastic.turbomachinery as tm

deg = 2 * np.pi / 360

eff_c = 0.6
eff_t = 0.65

dp_comp = 6250
dp_turb = -1780
dp_vac = 11560
mdot = 0.0623


def main():
    # DATA INPUT
    # inlet and outlet nodes
    inlet = tm.Node(101325, 293, 287)
    outlet = tm.Node(101325, 293, 287)

    # compressor blade rows (doesn't matter for now)
    com_rotor = tm.BladeRow(bi=0 * deg, bo=20 * deg, ri=0.02, ro=0.1, N=10)
    com_stator = tm.BladeRow(bi=-40 * deg, bo=-20 * deg, ri=0.11, ro=0.2, N=10)

    # turbine blade rows (doesn't matter for now)
    tur_rotor = tm.BladeRow(bi=0 * deg, bo=-20 * deg, ri=0.02, ro=0.1, N=10)
    tur_stator = tm.BladeRow(bi=40 * deg, bo=20 * deg, ri=0.11, ro=0.2, N=10)

    # turbomachinery elements
    com = tm.Compressor(com_rotor, com_stator, h=0.012)
    tur = tm.Turbine(tur_rotor, tur_stator, h=0.012)
    vac = tm.Vacuum()

    # rig setup
    rig = tm.Rig(inlet, [com, tur, vac], outlet)
    rig.add_to_shaft(idxs=[0, 1])

    # COMPRESSOR DESIGN [notation: V_{(r)adial or (t)angential}_{(r)otor or (s)tator}{(i)nlet or (o)utlet}]
    dh0 = dp_comp / inlet.density / eff_c
    V_r_ri = mdot / 2 / np.pi / inlet.density / com.rotor.ri / com.h
    V_r_ro = mdot / 2 / np.pi / inlet.density / com.rotor.ro / com.h

    # solve quadratic for bo from omega
    w = 8000 * 2 * np.pi / 60

    def f(bo, w):
        return Polynomial(
            [-dh0, V_r_ro * com.rotor.ro * np.tan(bo), com.rotor.sigma * com.rotor.ro ** 2]
        ).roots()[1] - w

    # set parameters
    com.rotor.bo = scipy.optimize.root_scalar(f, w, bracket=(0.01 - np.pi / 2, np.pi / 2 - 0.01)).root
    rig.shaft.omega = f(com.rotor.bo, 0) 

    # find rotor vels
    V_t_ro = dh0 / w / com.rotor.ro
    U_ri = w * com.rotor.ri
    U_ro = w * com.rotor.ro
    # W is relative velocity
    W_ri = np.sqrt(V_r_ri ** 2 + U_ri ** 2)
    W_ro = np.sqrt(V_r_ro ** 2 + (V_t_ro - U_ro) ** 2)

    # find bi: beta_in = (alpha_in-5deg) as per handout
    com.rotor.bi = np.arctan(U_ri / V_r_ri) - 5 * deg
    #com.rotor.bi = np.arctan(-U_ri / V_r_ri) + 5 * deg

    # rule of thumb
    assert (W_ro / W_ri > 1 / 3)

    # results are ok to show!
    print(f'beta_1: {com.rotor.bi / deg:.2f}deg')
    print(f'beta_2: {com.rotor.bo / deg:.2f}deg')
    print(f'omega: {rig.shaft.omega * 30 / np.pi:.0f}rpm')

    # now for blade numbers
    # first get blade profile for constant deltaP across blade
    r = np.linspace(com.rotor.ri, com.rotor.ro, 100)
    c1 = 2*(np.tan(com.rotor.bo) - np.tan(com.rotor.bi)) / (com.rotor.ro ** 2 - com.rotor.ri ** 2)
    c2 = np.tan(com.rotor.bi) - c1 * com.rotor.ri ** 2
    b = np.arctan(c1/2 * r ** 2 + c2)
    # vels
    V_r_r = mdot / 2 / np.pi / inlet.density / r / com.h
    V_t_r = w * r + V_r_r * np.tan(b)
    W_r = V_r_r / np.cos(b)

    drVt_dr = scipy.interpolate.InterpolatedUnivariateSpline(r, r * V_t_r).derivative(1)
    Nb_min = mdot * drVt_dr(r) / (2 * inlet.density * W_r ** 2 * r * com.h)

    Nb = np.ceil(1.25 * np.max(Nb_min))
    print(f'num. compressor blades: {Nb}')
    print(c1,c2)

    for i in range(1):
        plt.polar(np.tan(b)*np.log(r) + c1 + i * 2 * np.pi / 12, r, 'black')
    plt.yticks([])
    plt.ylim(0, None)
    plt.show()


if __name__ == "__main__":
    main()
