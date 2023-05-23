import numpy as np  

#Operating point
m_dot = 0.061 #kg/s
rho = 1.2046 #kg/m3
suction = 11.85 * 10**3 #Pa
suction_power_start = 598.99 #W
delta_p_comp_start = 6.41 * 10**3 #Pa
T_loss = 0.07168

#Geometry
r_1 = 0.05
r_2 = 0.10 #m
r_4 = 0.155
blade_height = 12 * 10**(-3)#m
thickness = 0.5E-3 #m
deg = np.pi/180

#For iteration
eff_comp = 0.6 #initial compressr efficiency estimate
eff_turb = 0.65
sigma_0 = 0.85

#Design choice
beta2 = -46 * deg #-0.7
beta_4 = 59 * deg
G_val = 1.075
v1_theta = 0
N_diff = 12 #14
inlet_angle_rotor = 5
inlet_angle_diff = 3

