import matplotlib.pyplot as plt 
import numpy as np

# Define input variables
rho = 1.2046 #kg/m^3
h = 12 *10**(-3) #blade heigh in m [12 mm]
#omega = 8000 * 2*np.pi/60 #Taken from compresor calculations, same shaft
etha_t = 0.65 
etha_mech = 0.899894 
del_p_t = 18260 #[kPa]
del_p_c = 6410
m_dot = 0.061 # [kg/m^3]
r3 = 10*10**(-2) #radius at turbine rotor exit [10 cm], measured
r4 = 45 * 10**(-3) #minimum radius is 50 mm as the blades will protude into the exhaust chamber

# Get vtheta and omega, matching with compressor calculations
Vtheta2 = 85.2651560166644
v03_actual = Vtheta2/etha_mech
omega = 1053.211292599318

# Define radial velocity function
def vradial(r):
    return - m_dot/(2*np.pi*rho*r*h)

# Define function for sigma and beta iteration given n
def getsig(n, v03_actual, vr3, r3, omega):
    sigma = 0.85 #initial guess
    sigma_1 = 0

    while round(sigma_1,4) != round(sigma,4):
        v03_ideal = (1-sigma)*(omega*r3) + v03_actual
        beta3 = np.arctan((v03_ideal-omega*r3)/vr3)
        sigma_1 = sigma
        sigma = 1 - ((np.cos(beta3))**(0.5))/(n**0.7)
    return sigma, beta3, v03_ideal

# Define blade number function
def bladen(beta3, beta4, r3, r4, v03_actual):
    beta = np.arctan((np.tan(beta3) + np.tan(beta4))*0.5)
    r_mid = (r3 + r4)/2
    vrm = vradial(r_mid)
    wavg = vrm/np.cos(beta)
    drvdr = (r3*v03_actual-r4*v04)/(r3-r4)
    n = m_dot*drvdr/(2*rho*(wavg**2)*r_mid*h)
    n = round(1.25 * n) # 25% extra safety factor
    print("Rotor blades number = {}".format(n))
    return n

# Define Blade Shape Equation function
def bshape(beta3, beta4, r3, r4):
    # Find metal angles variation with radius
    c1 = 2*(np.tan(beta4)-np.tan(beta3))/(r4**2-r3**2)
    c2 = np.tan(beta3)-c1*r3**2/2
    print('tan(\u03B2) = {} r^2 + {}'.format(c1/2,c2))

    # Define radius and beta vectors
    r = np.linspace(r4,r3,1000)
    beta = np.arctan(c1*r**2/2 + c2)
    return c1, c2, r

# Define plotting function
def plot_turbine(n, beta3, beta4, r3, r4, ns, beta33, beta2, r33, r2):
    # Plot rotor
    print('Rotor')
    c1r, c2r, rr = bshape(beta3, beta4, r3, r4)

    thetar = c1r*rr**2/4 + c2r*np.log(rr)
    thetar = (thetar - thetar[0])/2

    plt.subplots(subplot_kw={'projection': 'polar'})
    
    for i in range(n):
        theta_i = thetar + i*2*np.pi/n
        plt.plot(theta_i, rr, c = "c")
    plt.grid(True)
    
    # Plot stator
    print('Stator')
    c1s, c2s, rs = bshape(beta2, beta33, r2, r33)
    thetas = c1s*rs**2/4 + c2s*np.log(rs)
    thetas = (thetas - thetas[0])/4
    
    for i in range(ns):
        theta_i = thetas + i*2*np.pi/ns
        plt.plot(theta_i, rs, c = "r")
    plt.grid(True)

    # Show rotor/stator exit an inlet
    circle = np.linspace(0,2*np.pi,1000)
    radius = r4*np.ones(len(circle))
    radius_out = r3*np.ones(len(circle))
    stator_in = r33*np.ones(len(circle))
    stator_out = r2*np.ones(len(circle))
    plt.plot(circle,radius,"--", color = "black")
    plt.plot(circle,radius_out,"--", c = 'black')
    plt.plot(circle,stator_in,"--", c = 'black')
    plt.plot(circle,stator_out,"--", c = 'black')
    plt.show()
    return rr, thetar, rs, thetas


# Define function for throath area
def throat(r, thet):
    min = 100 #set min
    #adjacent blade
    thet_1 = thet + 2*np.pi/n

    for i in range(len(r)):
        for j in range(len(r)):
            x1 = r[i] * np.cos(thet[i])
            y1 = r[i] * np.sin(thet[i])
            x2 = r[j] * np.cos(thet_1[j])
            y2 = r[j] * np.sin(thet_1[j])
            d = ((x2 - x1) ** 2 + (y2 - y1) ** 2) **  0.5
            if d < min:
                min = d
    min = min - 0.0005 # accounts for blade thickness
    return min

# Define angle calculation function using vel triangles approach
def angle_calc(n=0, beta3=0):
    v04 = 0
    vr3 = vradial(r3)
    vr4 = vradial(r4)
    
    if n == 0:
        n = 15 #initial guess
        sigma, beta3, v03_ideal = getsig(n,v03_actual, vr3, r3, omega)
    else:
        sigma, beta3, v03_ideal = getsig(n,v03_actual, vr3, r3, omega)

    alpha4 = np.arctan(omega*r4/vr4)
    beta4 = alpha4 - 5*np.pi/180 #3D flow calculation suggests beta4 is 5 degrees more negative than alpha4,rel
    #this angle deviation slip will be lower if more blades are used
     
    print('\u03C3\u2083 = {}'.format(sigma))
    print('\u03B2\u2083 = {}'.format(beta3))
    print('\u03B2\u2084 = {}'.format(beta4))

    return beta3, beta4, v03_ideal, v03_actual, v04, alpha4

# Define function to get number of stator baldes


beta3, beta4, v03_ideal, v03_actual, v04, alpha4 = angle_calc()
n = 15
n = bladen(beta3, beta4, r3, r4, v03_actual)
beta3, beta4, v03_actual, v03_ideal, v04, alpha4 = angle_calc(n,beta3)

# Stator Velocity Analysis
r33 = r3*1.05 # 5% of local radius for stator-rotor gap, radius at stator exit
vr33 = vradial(r33)
v033 = v03_actual*r3/r33
beta33 = np.arctan(v033/vr33)
v02 = 0 # assume flow enters radially
beta2 = 0
r2 = 15 *10**(-2) #radius at turbine stator entry

# Find number of stator blades
theta33 = np.tan(beta33)/r33 * 0.001
theta2 = np.tan(beta2)/r2 * 0.001
x1 = r33 * np.cos(theta33)
y1 = r33 * np.sin(theta33)
x2 = r2 * np.cos(theta2)
y2 = r2 * np.sin(theta2)
chord = ((x2 - x1) ** 2 + (y2 - y1) ** 2) **  0.5
pitch_to_chord = 0.75 #as per handout
pitch = pitch_to_chord*chord
ns = round(2*np.pi*r33/pitch)
print("Stator blades number {}".format(ns))

rr, thetar, rs, thetas = plot_turbine(n, beta3, beta4, r3, r4, ns, beta33, beta2, r33, r2)

# Throat area
print(n*throat(rr,thetar))
print(2*np.pi*r4*np.cos(alpha4))