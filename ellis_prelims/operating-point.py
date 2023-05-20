import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

EFF_C = 0.6
EFF_T = 0.65
# P_LOSS = 60.06721792343545
INLET_DIAM = 48e-3
DENSITY = 1.204675848

f_0 = 0.5 * (0.7 + 1)
nu = 100 * 1e-6
Cf = 0.0025
d_m = 0.025
rho = 1.205
r_2 = 0.1
r_1 = 0


def p_loss(omega):
    n = omega * 60 / 2 / np.pi
    M_0 = f_0 * 1e3 * (nu * n) ** (2 / 3) * d_m ** 3
    T_w = 0.2 * Cf * rho * omega ** 2 * np.pi * (r_2 ** 5 - r_1 ** 5)
    return (2 * M_0 + 4 * T_w) * omega


def main(rpm):
    omega = rpm*2*np.pi/60
    suction_df = pd.read_csv("../../prelims/suction-chics.csv")
    inlet_dp = -suction_df["Water Gauge /kPa"] * 1000
    vac_dp = -suction_df["Mercury Gauge /kPa"] * 1000
    vac_dp0 = vac_dp - inlet_dp
    A_in = np.pi * INLET_DIAM ** 2 / 4

    mdot = np.sqrt(2 * inlet_dp * DENSITY) * A_in

    W_vac = mdot * vac_dp0 / DENSITY
    W_vac_fit = np.poly1d(np.polyfit(mdot, W_vac, 3))
    mdot_W_vac_max = np.max(np.polyder(W_vac_fit).roots)
    W_vac_max = W_vac_fit(mdot_W_vac_max)

    eff_mech_max = (W_vac_max - p_loss(omega) / EFF_T) / (W_vac_max - p_loss(omega) * EFF_C)
    W_t_max = (EFF_T * W_vac_max) / (1 - eff_mech_max * EFF_C * EFF_T)
    W_c_max = eff_mech_max * W_t_max

    F = W_c_max / W_vac_max

    print(W_c_max, W_vac_max, W_t_max, eff_mech_max)

    print('AT MAX PERFORMANCE:')
    print(f'    Mass flow: {mdot_W_vac_max:.4f}kg/s')
    print(f'    Suction: {W_vac_max * DENSITY / mdot_W_vac_max / 1000:.2f}kPa ({W_vac_max:.2f}W)')
    print(f'    Compresor pressure rise: '
          f'{W_c_max * EFF_C * DENSITY / mdot_W_vac_max / 1000:.2f}kPa')
    print(f'    Turbine pressure drop: {W_t_max / EFF_T * DENSITY / mdot_W_vac_max / 1000:.2f}kPa')
    print(f'    Overall efficiency: {eff_mech_max * EFF_C * EFF_T:.4f}')
    print(f'    F: {F:.4f}')

    plt.plot(mdot, W_vac, 'x', label='$\\dot{W}_v$ (data)')
    plt.plot(mdot, W_vac_fit(mdot), label='$\\dot{W}_v$ (fit)')
    plt.vlines(mdot_W_vac_max, 0, 700, 'r')
    plt.hlines(W_vac_max, 0, 0.10, 'r')
    plt.ylabel('Suction Power, W')
    plt.xlabel('Mass flow, kg/s')
    plt.ylim(0, None)
    plt.xlim(0, None)
    plt.show()


if __name__ == "__main__":
    main(8000)
