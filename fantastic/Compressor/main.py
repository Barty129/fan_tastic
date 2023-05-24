import numpy as np
import globals as gb
import functions as fn
import scipy
from matplotlib import pyplot as plt

step = 1E-4
def main():
    #Find radial velocities
    V1_r = fn.v_radial(gb.m_dot, gb.r_1, gb.rho, gb.blade_height)
    V2_r = fn.v_radial(gb.m_dot, gb.r_2, gb.rho, gb.blade_height)

    # set parameters
    deg = np.pi/180
    beta_2 = gb.beta2
    #Calculate omega by iteration
    omega_first = fn.omega_opt(gb.T_loss, gb.eff_comp, gb.eff_turb, gb.suction_power_start, \
                            gb.suction, gb.delta_p_comp_start, V2_r, beta_2, gb.r_2, gb.sigma_0, gb.rho)

    #Velocity triangles first pass (incorrect slip)
    beta_1, v2_theta, v1_rel, v2_rel = fn.velocity2_triangles(gb.r_1, gb.r_2, V1_r, V2_r, beta_2, omega_first, 0)
    
    # Calculate blade numbers, using midpoint analysis
    r_mid = 0.5*(gb.r_1 + gb.r_2)
    V_r_mid = fn.v_radial(gb.m_dot, r_mid, gb.rho, gb.blade_height)
    beta_mid = np.arctan(0.5*(np.tan(beta_1) + np.tan(beta_2)))

    V_theta_vals = (omega_first * r_mid) + (V_r_mid * np.tan(beta_mid))
    W_av = V_r_mid/np.cos(beta_mid)
    dvthetdr = ((v2_theta * gb.r_2) - (gb.v1_theta * gb.r_1))/(gb.r_2 - gb.r_1)

    Nb_min = gb.m_dot * dvthetdr / (2 * gb.rho * r_mid * gb.blade_height * (W_av**2))
    Nb = round(1.25 * np.max(Nb_min))
    if Nb > 21:
        Nb = 21
    Nb = int(Nb)
    print('No. Blades = {}'.format(Nb))

    #Velocity triangle second pass (corrected slip)
    sigma2 = 1 - (np.sqrt(np.cos(beta_2))/(Nb**(0.7)))
    omega = fn.omega_opt(gb.T_loss, gb.eff_comp, gb.eff_turb, gb.suction_power_start, \
                    gb.suction, gb.delta_p_comp_start, V2_r, beta_2, gb.r_2, sigma2, gb.rho)
    beta_1, v2_theta, v1_rel, v2_rel = fn.velocity2_triangles(gb.r_1, gb.r_2, V1_r, V2_r, beta_2, omega, Nb)
    print('beta_1 = {} deg'.format(beta_1/deg))
    print('beta_2 = {} deg'.format(beta_2/deg))
    print('omega = {} rpm'.format(omega * 60/(2*np.pi)))
    print('sigma = {}'.format(sigma2))

    r_rot, theta_rot = fn.blades_plot(gb.r_1, gb.r_2, beta_1, beta_2, step)

    #Stator design
    r_3 = gb.r_2 * gb.G_val
    v3_theta = (gb.r_2 * v2_theta)/r_3
    v3_r = fn.v_radial(gb.m_dot, r_3, gb.rho, gb.blade_height)
    beta_3 = np.arctan(v3_theta/v3_r) - gb.inlet_angle_diff * deg

    rdiff, thetadiff = fn.stator_blades(r_3, gb.r_4, beta_3, gb.beta_4, step)
    (throat_min, index_1) = fn.throat_dist(gb.N_diff, thetadiff, rdiff, gb.blade_height, gb.r_4, 0)
    (throat_max, index_2) = fn.throat_dist(gb.N_diff, thetadiff, rdiff, gb.blade_height, gb.r_4, 1)
    throat_min = (throat_min - (gb.thickness*gb.blade_height)) #*0.95 For boundary layers
    area_ratio = throat_max/throat_min

    #Diffuser
    blade_length = fn.blade_length(rdiff, thetadiff, index_1[1])
    diffuser_angle = np.arctan((throat_max - throat_min)/(gb.blade_height*blade_length))

    print('beta_3 = {} deg'.format(beta_3/deg))
    print('beta_4 = {} deg'.format(gb.beta_4/deg))
    print('Area ratio = {}'.format(area_ratio))
    print('Diffuser Angle = {}'.format(diffuser_angle/deg))

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    angle = np.linspace(0,2*np.pi,1000)
    r1_array = gb.r_1*np.ones(len(angle))
    r2_array = gb.r_2*np.ones(len(angle))
    plt.plot(angle,r1_array,"--", c = 'black')
    plt.plot(angle,r2_array,"--", c = 'black')
    for i in range(gb.N_diff):
        theta_plot = thetadiff + (i * 2 *np.pi/gb.N_diff)
        ax.plot(theta_plot, rdiff, 'blue')
    for i in range(Nb):
        theta_plot = theta_rot + (i * 2 *np.pi/Nb)
        ax.plot(theta_plot, r_rot, 'red')
    r3_array = r_3*np.ones(len(angle))
    r4_array = gb.r_4*np.ones(len(angle))
    plt.plot(angle,r3_array,"--", c = 'black')
    plt.plot(angle,r4_array,"--", c = 'black')
    plt.yticks([])
    plt.ylim(0, None)
    plt.show()



if __name__ == "__main__":
    main()