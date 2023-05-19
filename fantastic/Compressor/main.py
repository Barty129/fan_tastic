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
    deg = np.pi/180

    omega_0 = 8000 * 2 * np.pi / 60

    # set parameters
    beta_2 = scipy.optimize.root_scalar(fn.omega_opt, args=(V1_r, gb.r_2, gb.sigma_0, dh0_comp, omega_0) , bracket=(0.01 - np.pi / 2, np.pi / 2 - 0.01)).root
    omega = fn.omega_opt(beta_2, V1_r, gb.r_2, gb.sigma_0, dh0_comp, 0) 

    v2_theta = fn.euler(dh0_comp, omega, gb.r_1, gb.r_2, gb.v1_theta)
    Ub_1 = omega * gb.r_1
    Ub_2 = omega * gb.r_2
    # Calculate relative velocities
    v1_t_rel = np.sqrt(V1_r ** 2 + Ub_1 ** 2)
    v2_t_rel = np.sqrt(V2_r ** 2 + (v2_theta - Ub_2) ** 2)

    # find bi: beta_in = (alpha_in-5deg) as per handout
    beta_1 = np.arctan(-Ub_1 / V1_r) + 5 * deg

    #De Haller test
    if v2_t_rel/v1_t_rel < 1/3:
        text = 'This is not valid'
        print(text)
    else:
        print('beta_1 = {} deg'.format(beta_1/deg))
        print('beta_2 = {} deg'.format(beta_2/deg))
        print('omega = {} rpm'.format(omega * 60/(2*np.pi)))
    
    # Calculate blade numbers, using midpoint analysis
    r_mid = 0.5*(gb.r_1 + gb.r_2)
    r_diff = 0.5 * (gb.r_2 - gb.r_1)
    r_mid_vals = np.linspace((r_mid - 0.3*r_diff), (r_mid + 0.3*r_diff), 100)
    V_r_vals = fn.v_radial(gb.m_dot, r_mid_vals, gb.rho, gb.blade_height)

    c1_mid = 2*(np.tan(beta_2) - np.tan(beta_1)) / (gb.r_2 ** 2 - gb.r_1 ** 2)
    c2_mid = np.tan(beta_1) - c1_mid * gb.r_1 ** 2
    beta_mid_vals = np.arctan(c1_mid/2 * r_mid_vals ** 2 + c2_mid)

    V_theta_vals = (omega * r_mid_vals) + (V_r_vals * np.tan(beta_mid_vals))
    W_av_vals = V_r_vals/np.cos(beta_mid_vals)
    dvthetdr = scipy.interpolate.InterpolatedUnivariateSpline(r_mid_vals, r_mid_vals * V_theta_vals).derivative(1)

    Nb_min = gb.m_dot * dvthetdr(r_mid_vals) / (2 * gb.rho * r_mid_vals * gb.blade_height * (W_av_vals**2))
    Nb = np.ceil(1.25 * np.max(Nb_min))
    if Nb > 21:
        Nb = 21
    elif Nb < 7:
        Nb = 7
    print('No. Blades = {}'.format(Nb))

    r_vals = np.linspace(gb.r_1, gb.r_2, 100)
    c1 = 2*(np.tan(beta_2) - np.tan(beta_1)) / (gb.r_2 ** 2 - gb.r_1 ** 2)
    c2 = np.tan(beta_1) - c1/2 * gb.r_1 ** 2
    beta_vals_tan = c1/2 * r_vals ** 2 + c2

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    for i in range(1):
        ax.plot(beta_vals_tan*np.log(r_vals) + c2 +  i * 2 * np.pi / int(Nb), r_vals, 'blue')
    plt.yticks([])
    plt.ylim(0, None)
    plt.show()

    r_3 = gb.r_2 * gb.G_val

if __name__ == "__main__":
    main()