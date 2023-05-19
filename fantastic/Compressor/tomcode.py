import numpy as np
import matplotlib.pyplot as plt

# design parameters to be altered
beta2 = -1 # trailing edge blade angle in radians
r1 = 0.03 # leading edge radius
inletangle = 5 # angle of incidence of leading edge in degrees
vaneless = 1.1 # ratio of diffuser leading edge to rotor trailing edge radii
arearatio = 3
aspectratio = 10
# numerical parameters
step = 0.001

# design point values
nc = 0.6
h = 0.012
dp0 = 6410
density = 1.2046
mdot = 0.061
r2 = 0.1
dh0 = 1/nc * dp0/density

def bladeplot(N, r1, r2, step, beta1, beta2):
    # plots blade shapes
    r = np.arange(r1, r2, 0.001)
    theta = [[0 for j in range(len(r))] for i in range(N)]
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    beta = [0 for i in range(len(r))]
    c = np.tan(beta2)*2/r2**2
    for i in range(len(r)):
        beta[i] = -1*np.arctan(c*r[i]**2/2)
    for i in range(0, N):   
        for j in range(0, len(r)):
            theta[i][j] = i * 2 * np.pi/N + np.tan(beta[j])*np.log(r[j]) - np.tan(beta[0])*np.log(r[0])
        ax.plot(theta[i], r)
    ax.grid(True)
    plt.show()

def velocitycalc(N = 0):
    # calculates velocity triangles
    if N != 0:
        sigma = 1 - (np.cos(beta2)**0.5)/N**0.7
    else:
        sigma = 0.85
    Vr2 = mdot/(2 * np.pi * density * r2 * h)
    a = sigma * r2 ** 2
    b = r2 * Vr2 * np.tan(beta2)
    c = -dh0
    omega = (-b + (b**2 - 4*a*c)**0.5)/(2*a)
    Vtheta2 = sigma * omega * r2 + Vr2 * np.tan(beta2)
    alphaRel2 = np.arctan(Vtheta2 - omega * r2)/Vr2
    alpha2 = Vtheta2/Vr2
    Vtheta1 = omega * r1
    Vr1 = mdot/(2 * np.pi * density * r1 * h)
    Vtheta1 = Vtheta2 * r1/r2
    alphaRel1 = np.arctan(-Vtheta1/Vr1)
    if alphaRel1 >= 0:
        beta1 = alphaRel1 - inletangle * np.pi/180
    else:
        beta1 = alphaRel1 + inletangle * np.pi/180
    Vrel1 = (Vr1**2 + Vtheta1**2)**0.5
    Vrel2 = (Vr2**2 + Vtheta2**2)**0.5
    return(Vr2, omega, Vtheta1, Vtheta2, Vr1, alphaRel1, beta1, Vrel1, Vrel2)

def bladenumber():
    # calculates minimum number of blades
    rmid = (r2 + r1) * 0.5
    betamid = np.arctan(0.5*(np.tan(beta1) + np.tan(beta2)))
    drVthetadr = (r2 * Vtheta2 - r1 * Vtheta1)/(r2 - r1)
    Vr = mdot/(2 * np.pi * density * rmid * h)
    Wavg = Vr/np.cos(betamid)
    Nb = mdot * drVthetadr / (2 * density * Wavg ** 2 * rmid * h)
    return(Nb)

Vr2, omega, Vtheta1, Vtheta2, Vr1, alphaRel1, beta1, Vrel1, Vrel2 = velocitycalc()
N = bladenumber()
N = round(N * 1.25)
Vr2, omega, Vtheta1, Vtheta2, Vr1, alphaRel1, beta1, Vrel1, Vrel2 = velocitycalc(N) # recalculates values using a refined blade number

bladeplot(N, r1, r2, step, beta1, beta2)

r3 = r2 * vaneless
Vr3 = Vr2 * r2/r3
Vtheta3 = Vtheta2 * r2/r3
beta1 = -1 * np.arctan(Vtheta3/Vr3)
