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
    omega = fn.omega_opt(gb.T_loss, gb.eff_comp, gb.eff_turb, gb.suction_power_start, \
                            gb.suction, gb.delta_p_comp_start, V2_r, beta_2, gb.r_2, gb.sigma_0, gb.rho)

    #Velocity triangles first pass (incorrect slip)
    beta_1, v2_theta, v1_rel, v2_rel, sigma1 = fn.velocity2_triangles(gb.r_1, gb.r_2, V1_r, V2_r, beta_2, omega, 0)
    
    # Calculate blade numbers, using midpoint analysis
    r_mid = 0.5*(gb.r_1 + gb.r_2)
    V_r_mid = fn.v_radial(gb.m_dot, r_mid, gb.rho, gb.blade_height)
    beta_mid = np.arctan(0.5*(np.tan(beta_1) + np.tan(beta_2)))

    V_theta_vals = (omega * r_mid) + (V_r_mid * np.tan(beta_mid))
    W_av = V_r_mid/np.cos(beta_mid)
    dvthetdr = ((v2_theta * gb.r_2) - (gb.v1_theta * gb.r_1))/(gb.r_2 - gb.r_1)

    Nb_min = gb.m_dot * dvthetdr / (2 * gb.rho * r_mid * gb.blade_height * (W_av**2))
    Nb = np.ceil(1.25 * np.max(Nb_min))
    Nb = int(Nb)
    print('No. Blades = {}'.format(Nb))

    #Velocity triangle second pass (corrected slip)
    beta_1, v2_theta, v1_rel, v2_rel, sigma2 = fn.velocity2_triangles(gb.r_1, gb.r_2, V1_r, V2_r, beta_2, omega, Nb)
    print('beta_1 = {} deg'.format(beta_1/deg))
    print('beta_2 = {} deg'.format(beta_2/deg))
    print('omega = {} rpm'.format(omega * 60/(2*np.pi)))
    print('sigma = {}'.format(sigma2))
    
    fn.rotor_plot(gb.r_1, gb.r_2, beta_1, beta_2, step, Nb)

    #Stator design
    r_3 = gb.r_2 * gb.G_val
    v3_theta = (gb.r_2 * v2_theta)/r_3
    v3_r = fn.v_radial(gb.m_dot, r_3, gb.rho, gb.blade_height)
    beta_3 = np.arctan(v3_theta/v3_r) - gb.inlet_angle * deg
    beta_4 = 1

    rdiff, thetadiff = fn.stator_blades(r_3, gb.r_4, step, beta_3, beta_4, gb.N_diff)
    (throat_min, index_1, theta_off) = fn.throat_dist(gb.N_diff, thetadiff, rdiff, gb.blade_height, gb.r_4, 0)
    (throat_max, index_2, theta_off) = fn.throat_dist(gb.N_diff, thetadiff, rdiff, gb.blade_height, gb.r_4, 1)
    throat_min = (throat_min - (gb.thickness*gb.blade_height)) #*0.95 For boundary layers
    area_ratio = throat_max/throat_min
    r_exit = np.array([(rdiff[index_2[0]] * np.cos(thetadiff[index_2[0]])), (rdiff[index_2[0]] * np.sin(thetadiff[index_2[0]]))])
    r_throat = np.array([(rdiff[0] * np.cos(thetadiff[0])), (rdiff[0] * np.sin(thetadiff[0]))]) #End point
    delta_blade = r_exit - r_throat
    blade_length = np.sqrt((delta_blade[0]**2) + (delta_blade[1]**2))
    print(blade_length)
    diffuser_angle = np.arctan((throat_max - throat_min)/(gb.blade_height*blade_length))

    # length_ratio = bladelen*gb.blade_height/throat_min
    #div_angle = 360/gb.N_diff

    print('beta_3 = {} deg'.format(beta_3/deg))
    print('beta_4 = {} deg'.format(beta_4/deg))
    print('Area ratio = {}'.format(area_ratio))
    print('Diffuser Angle = {}'.format(diffuser_angle/deg))
    # print('Length ratio = {}'.format(length_ratio))

if __name__ == "__main__":
    main()