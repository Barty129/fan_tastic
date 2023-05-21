import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import root_scalar

from fantastic import turbomachinery as tm

eff_c = 0.6
eff_t = 0.65
mdot = 0.0623

deg = 2 * np.pi / 360


def V_r(r, rho, h):
    return mdot / (2 * np.pi * rho * r * h)


def main(rpm):
    omega = rpm * 2 * np.pi / 60
    # DATA INPUT
    # inlet and outlet nodes
    inlet = tm.Node(101325, 293, 287)
    outlet = tm.Node(101325, 293, 287)

    # compressor blade rows (doesn't matter for now)
    com_rotor = tm.BladeRow(bi=0, bo=0, ri=0.025, ro=0.100, N=10)
    com_stator = tm.BladeRow(bi=0, bo=0, ri=0.110, ro=0.200, N=10)

    # turbine blade rows (doesn't matter for now)
    tur_rotor = tm.BladeRow(bi=0, bo=0, ri=0.045, ro=0.100, N=10)
    tur_stator = tm.BladeRow(bi=0, bo=0, ri=0.110, ro=0.200, N=10)

    # turbomachinery elements
    com = tm.Compressor(com_rotor, com_stator, h=0.012)
    tur = tm.Turbine(tur_rotor, tur_stator, h=0.012)

    # ROTOR
    dh0_c, dh0_t = tm.dh0(omega)

    V_t_ro = dh0_t / omega / tur.rotor.ro
    V_r_ro = mdot / (2 * np.pi * inlet.density * tur.rotor.ro * tur.h)
    V_r_ri = mdot / (2 * np.pi * inlet.density * tur.rotor.ri * tur.h)

    def f(x):
        return omega * tur.rotor.ro * (1 - np.cos(x) ** 0.5 / tur.rotor.N ** 0.7) + V_r_ro * np.tan(x) - V_t_ro

    tur.rotor.bo = root_scalar(f, bracket=[0.1 - np.pi / 2, np.pi / 2 - 0.1]).root
    tur.rotor.bi = np.arctan(-omega * tur.rotor.ri / V_r_ri) - 5 * deg
    tur.rotor.sigma = 1 - np.cos(tur.rotor.bo) ** 0.5 / tur.rotor.N ** 0.7

    r_mid = (tur.rotor.ro + tur.rotor.ri) / 2
    b_mid = np.arctan((np.tan(tur.rotor.bi) + np.tan(tur.rotor.bo)) / 2)
    V_r_rm = V_r(r_mid, inlet.density, tur.h)
    W_rm = V_r_rm / np.cos(b_mid)

    drVt_dr = (tur.rotor.ro * V_t_ro) / (tur.rotor.ro - tur.rotor.ri)
    Nb_min = mdot * drVt_dr / (2 * inlet.density * W_rm ** 2 * r_mid * com.h)

    Nb = int(np.ceil(1.25 * np.max(Nb_min)))
    tur.rotor.N = Nb

    throat_width = 2 * np.pi * tur.rotor.ri * np.cos(tur.rotor.bi + 5 * deg) / Nb

    print(f'beta_3: {tur.rotor.bo / deg:.2f}deg\n'
          f'beta_4: {tur.rotor.bi / deg:.2f}deg\n'
          f'Nb: {Nb}\n'
          f'throat area: {throat_width * tur.h*1e6:.2f}mm2\n'
          f'resulting exit diameter: {np.sqrt(4 * throat_width * tur.h / np.pi) * 1000:.0f}mm')

    r, b, theta = tm.rotor_shape(tur.rotor)
    for i in range(tur.rotor.N):
        plt.polar(-theta + i * 2 * np.pi / tur.rotor.N, r, 'black')

    # STATOR

    plt.yticks([])
    plt.ylim(0, None)
    plt.show()


if __name__ == "__main__":
    main(9600)
