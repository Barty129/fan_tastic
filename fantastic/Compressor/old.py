import numpy as np
from numpy.polynomial import Polynomial
import globals as gb
from matplotlib import pyplot as plt
import numpy as np
import globals as gb
import functions as fn
import cmath

def stator_blades(r_3, r_4, beta_3, beta_4, step):
    # plots diffuser blade shapes
    r_val = np.arange(r_3, r_4, step)
    theta = np.zeros(len(r_val))
    b = beta_3 + (beta_4 - beta_3) * (r_val - r_3)/(r_4 - r_3)
    for j in range(1, len(r_val)):
        theta[j] = theta[j-1] + np.tan(b[j])/r_val[j] * step
    return(r_val, theta) 

''''
    file = open("Rot9.txt","w+")
    theta_plot = theta_rot + (8*2*np.pi)/int(Nb)
    for j in range(len(r_rot)):
        x = r_rot[j] * np.cos(theta_plot[j])
        y = r_rot[j] * np.sin(theta_plot[j])
        z = 0
        file.write(str(x))
        file.write(',')
        file.write(str(y))
        file.write(',')
        file.write(str(z))
        file.write('\n')
    file.close()
'''

'''
deg = np.pi/180
for N in range(7,21):
    beta_test = np.linspace(-85, 70, 100) *deg
    limit = 1/3 * np.ones(len(beta_test))
    omega = np.ones(len(beta_test))
    V_ratio = np.ones(len(beta_test))
    sigma2 = np.ones(len(beta_test))
    V1_r = fn.v_radial(gb.m_dot, gb.r_1, gb.rho, gb.blade_height)
    V2_r = fn.v_radial(gb.m_dot, gb.r_2, gb.rho, gb.blade_height)
    for i in range(len(beta_test)):
        sigma2[i] = 1 - (np.sqrt(np.cos(beta_test[i]))/(N**(0.7)))
        omega[i] = fn.omega_opt(gb.T_loss, gb.eff_comp, gb.eff_turb, gb.suction_power_start, \
                            gb.suction, gb.delta_p_comp_start, V2_r, beta_test[i], gb.r_2, sigma2[i], gb.rho)
        beta_1, v2_theta, v1_rel, v2_rel = fn.velocity2_triangles(gb.r_1, gb.r_2, V1_r, V2_r, beta_test[i], omega[i], N)
        V_ratio[i] = v2_rel/v1_rel
    plt.plot(beta_test/deg , V_ratio, label = 'N = {}'.format(N))
plt.plot(beta_test/deg , limit, 'black')
plt.legend()
plt.ylim(0,1)
plt.xlabel(r'$\/beta_2 (degrees)$')
plt.ylabel('Velocity Ratio')
plt.show()
'''


'''
    N_range = np.linspace(8,16,20)
    aratio = np.zeros(len(N_range))
    limit = 3 * np.ones(len(N_range))


    for N_b in range(len(N_range)):
        r_3 = gb.r_2 * gb.G_val
        v3_theta = (gb.r_2 * v2_theta)/r_3
        v3_r = fn.v_radial(gb.m_dot, r_3, gb.rho, gb.blade_height)
        beta_3 = np.arctan(v3_theta/v3_r) - gb.inlet_angle_diff * deg

        rdiff, thetadiff = fn.blades_plot(r_3, gb.r_4, beta_3, gb.beta_4, step)
        (throat_min, index_1) = fn.throat_dist(N_range[N_b], thetadiff, rdiff, gb.blade_height, gb.r_4, 0)
        (throat_max, index_2) = fn.throat_dist(N_range[N_b], thetadiff, rdiff, gb.blade_height, gb.r_4, 1)
        throat_min = (throat_min - (gb.thickness*gb.blade_height)) #*0.95 For boundary layers
        area_ratio = throat_max/throat_min
        aratio[N_b] = area_ratio
        print(area_ratio)
    plt.plot(N_range, aratio)
    plt.plot(N_range , limit, 'black')
    plt.xlabel('No. of Diffuser Blades')
    plt.ylabel('Area Ratio')
    plt.show()
'''

'''
N_range = np.linspace(6,16,20)
aratio = np.zeros(len(N_range))
limit = 3 * np.ones(len(N_range))
divangle = np.zeros(len(N_range))
limit2 = 12 * np.ones(len(N_range))


for N_b in range(len(N_range)):
    r_3 = gb.r_2 * gb.G_val
    v3_theta = (gb.r_2 * v2_theta)/r_3
    v3_r = fn.v_radial(gb.m_dot, r_3, gb.rho, gb.blade_height)
    beta_3 = np.arctan(v3_theta/v3_r) - gb.inlet_angle_diff * deg

    rdiff, thetadiff = fn.blades_plot(r_3, gb.r_4, beta_3, gb.beta_4, step)
    (throat_min, index_1) = fn.throat_dist(N_range[N_b], thetadiff, rdiff, gb.blade_height, gb.r_4, 0)
    (throat_max, index_2) = fn.throat_dist(N_range[N_b], thetadiff, rdiff, gb.blade_height, gb.r_4, 1)
    throat_min = (throat_min - (gb.thickness*gb.blade_height)) #*0.95 For boundary layers
    area_ratio = throat_max/throat_min
    blade_length = fn.blade_length(rdiff, thetadiff, index_1[1])
    diffuser_angle = np.arctan((throat_max - throat_min)/(gb.blade_height*blade_length))
    divangle[N_b] = diffuser_angle/deg
    print(diffuser_angle/deg)
plt.plot(N_range, divangle)
plt.plot(N_range , limit2, 'black')
plt.xlabel('No. of Diffuser Blades')
plt.ylabel('Diffuser Angle (degrees)')
plt.show()
'''