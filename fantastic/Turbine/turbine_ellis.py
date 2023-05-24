import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import root_scalar

from fantastic import turbomachinery as tm

eff_c = 0.6
eff_t = 0.65
mdot = 0.0623

PITCH_TO_CHORD = 0.5

deg = 2 * np.pi / 360


def V_r(r, rho, h):
    return mdot / (2 * np.pi * rho * r * h)



def main(rpm):
    omega = rpm * 2 * np.pi / 60
    # DATA INPUT
    # inlet and outlet nodes
    air = tm.Node(101325, 293, 287)

    # turbine blade rows (doesn't matter for now)
    tur_rotor = tm.BladeRow(bi=0, bo=0, ri=0.050, ro=0.100, N=16, rotor=True)
    tur_stator = tm.BladeRow(bi=0, bo=0, ri=0.105, ro=0.155, N=12, rotor=False)

    # turbine
    tur = tm.Turbine(tur_rotor, tur_stator, h=0.012)

    # ROTOR
    # get enthalpy change across turbine
    dh0_c, dh0_t = tm.dh0(omega)

    # get velocities from geometry where possible
    V_t_ro = dh0_t / omega / tur.rotor.ro
    V_r_ro = V_r(tur.rotor.ro, air.density, tur.h)
    V_r_ri = V_r(tur.rotor.ri, air.density, tur.h)

    # beta4 comes from the exit condition having zero swirl, and beta = alpha-5deg
    tur.rotor.bi = np.arctan(-omega * tur.rotor.ri / V_r_ri) - 5 * deg

    # beta outer
    # function to be solved when beta_3, N and the velocities are compatible
    def f(x, N):
        return omega * tur.rotor.ro * (1 - np.cos(x) ** 0.5 / N ** 0.7) + V_r_ro * np.tan(x) - V_t_ro

    # beta3 comes from solving the above function
    tur.rotor.bo = root_scalar(f, args=(tur.rotor.N,), bracket=[0.1 - np.pi / 2, np.pi / 2 - 0.1]).root

    # midpoint analysis for checking minimum number of blades
    r_mid = (tur.rotor.ro + tur.rotor.ri) / 2
    b_mid = np.arctan((np.tan(tur.rotor.bi) + np.tan(tur.rotor.bo)) / 2)
    V_r_rm = V_r(r_mid, air.density, tur.h)
    W_rm = V_r_rm / np.cos(b_mid)

    # estimate of derivative using values at stations 3 and 4
    drVt_dr = (tur.rotor.ro * V_t_ro) / (tur.rotor.ro - tur.rotor.ri)
    Nb_min = mdot * drVt_dr / (2 * air.density * W_rm ** 2 * r_mid * tur.h)

    # number of blades comes from taking the above minimum, adding 25% then clipping to a sensible range
    tur.rotor.N = int(np.clip(1.25 * np.max(Nb_min), 8, 21))

    # recalculate beta 3 and sigma with new bladenumber
    tur.rotor.bo = root_scalar(f, args=(tur.rotor.N,), bracket=[0.1 - np.pi / 2, np.pi / 2 - 0.1]).root
    tur.rotor.sigma = 1 - np.cos(tur.rotor.bo) ** 0.5 / tur.rotor.N ** 0.7
    # generate shape
    tur.rotor.generate_rotor()

    # calculate throat area for later checks
    # throat_area = 2 * np.pi * tur.rotor.ri * np.cos(tur.rotor.bi + 5 * deg) * tur.h

    r = np.linspace(tur.rotor.ri, tur.rotor.ro, 100)
    circle = np.linspace(-np.pi, np.pi, 100)
    plt.polar(circle, tur.rotor.ri * np.ones_like(circle), 'k--')
    plt.polar(circle, tur.rotor.ro * np.ones_like(circle), 'k--')
    for i in range(tur.rotor.N):
        plt.polar(tur.rotor.theta(r) + i * 2 * np.pi / tur.rotor.N, r, 'green')

    # STATOR
    tur.stator.bo = 0
    V_t_si = V_t_ro * tur.rotor.ro / tur.stator.ri
    V_r_si = V_r(tur.stator.ri, air.density, tur.h)
    tur.stator.bi = np.arctan(V_t_si / V_r_si)

    pitch = PITCH_TO_CHORD * tur.stator.chord
    tur.stator.N = int(np.ceil(np.clip(2 * np.pi * tur.stator.ri / pitch, 8, 21)))

    r = np.linspace(tur.stator.ri, tur.stator.ro, 100)
    plt.polar(circle, tur.stator.ri * np.ones_like(circle), 'k--')
    plt.polar(circle, tur.stator.ro * np.ones_like(circle), 'k--')
    for i in range(tur.stator.N):
        plt.polar(tur.stator.theta(r) + i * 2 * np.pi / tur.stator.N, r, 'red')

    tur.rotor.plot_xy("turbine-rotor.txt")
    tur.stator.plot_xy("turbine-stator.txt")

    print(f'ROTOR:\n'
          f'    beta-outer: {tur.rotor.bo / deg:.2f}°\n'
          f'    beta-inner: {tur.rotor.bi / deg:.2f}°\n'
          f'    Nb: {tur.rotor.N}\n'
          f'STATOR\n'
          f'    beta-outer: {tur.stator.bo / deg:.2f}°\n'
          f'    beta-inner: {tur.stator.bi / deg:.2f}°\n'
          f'    Nb: {tur.stator.N}\n'
          f'GENERAL:\n'
          f'    rpm: {rpm}')

    plt.yticks([])
    plt.xticks([])
    plt.ylim(0, None)
    plt.show()


if __name__ == "__main__":
    main(9290)
