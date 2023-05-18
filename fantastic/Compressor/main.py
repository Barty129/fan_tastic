import numpy as np
import globals as gb
import functions as fn
import scipy
from matplotlib import pyplot as plt

def main():
    dh0_comp = fn.dh0(gb.delta_p_comp,gb.rho,gb.eff_comp)
    V1_r = fn.v_radial(gb.m_dot, gb.r_1, gb.rho, gb.blade_height)
    V2_r = fn.v_radial(gb.m_dot, gb.r_2, gb.rho, gb.blade_height)

    # set parameters
    beta_2 = 1
    omega = fn.omega_opt(beta_2, V2_r, gb.r_2, gb.sigma_0, dh0_comp)

    V_t_ro = fn.euler
    Ub_1 = omega * gb.r_1
    Ub_2 = omega * gb.r_2
    # W is relative velocity
    W_ri = np.sqrt(V1_r ** 2 + Ub_1 ** 2)
    W_ro = np.sqrt(V2_r ** 2 + (V_t_ro - Ub_2) ** 2)

if __name__ == "__main__":
    main()
    
