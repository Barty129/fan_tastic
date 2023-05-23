from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad


@dataclass
class Node:
    p: float
    T: float
    R: float

    @property
    def density(self):
        return self.p / self.R / self.T


@dataclass
class BladeRow:
    """
    Turbomachinery blade row

    :parameter bi: inside blade angle
    :parameter bo: outside blade angle
    :parameter ri: inside radius
    :parameter ro: outside radius
    :parameter N: number of blades
    :parameter sigma: slip default 0.85
    """
    bi: float
    bo: float
    ri: float
    ro: float
    N: int

    sigma: float = 0.85

    c1: float = 0
    c2: float = 0
    c3: float = 0

    @property
    def chord(self):
        def f(r):
            return np.sqrt(1 + (self.c1 * r ** 2 / 2 + self.c3 * r + self.c2) ** 2)

        return quad(f, self.ri, self.ro)[0]

    def b(self, r):
        return np.arctan(self.c1 * r ** 2 / 2 + self.c3 * r + self.c2)

    def theta(self, r, rotor=True):
        if rotor:
            return self.c1 * r ** 2 / 4 + self.c3 * r + self.c2 * np.log(r)
        else:
            return np.array([quad(lambda x: np.tan(self.bi + (x - self.ri) * (self.bo - self.bi) / (self.ro - self.ri)) / x,
                                  self.ri,
                                  r_i
                                  )[0] for r_i in r])

    def generate_rotor(self):
        self.c1 = 2 * (np.tan(self.bi) - np.tan(self.bo)) / (self.ri ** 2 - self.ro ** 2)
        self.c2 = np.tan(self.bi) - self.c1 * self.ri ** 2 / 2


class Turbomachine:
    def __init__(self):
        self.source = None
        self.sink = None
        self.shaft_connection = None
        self.omega = 0

    @property
    def torque(self):
        return 0


class Compressor(Turbomachine):
    def __init__(self, rotor: BladeRow, stator: BladeRow, h: float):
        super().__init__()
        self.h = h
        if rotor.ro < stator.ri:
            self.rotor = rotor
            self.stator = stator
        else:
            raise ValueError("rotor and stator incompatible, negative gap size")


class Turbine(Turbomachine):
    def __init__(self, rotor: BladeRow, stator: BladeRow, h: float):
        super().__init__()
        self.h = h
        if rotor.ro < stator.ri:
            self.rotor = rotor
            self.stator = stator
        else:
            raise ValueError("rotor and stator incompatible, negative gap size")


class Vacuum(Turbomachine):
    def __init__(self):
        super().__init__()


class Shaft:
    def __init__(self, omega: float, turbomachines: list[Turbomachine]):
        self._omega = omega
        self.turbomachines = turbomachines

        self.omega = omega

    @property
    def torque(self):
        return sum(tm.torque for tm in self.turbomachines)

    @property
    def omega(self):
        return self._omega

    @omega.setter
    def omega(self, value):
        self._omega = value
        for tm in self.turbomachines:
            tm.omega = value


class Rig:
    def __init__(self, source: Node, turbomachines: list[Turbomachine], sink: Node):
        self.source = source
        self.turbomachines = turbomachines
        self.sink = sink

        self.connect_fluids()

        self.shaft = Shaft(0, [])

    def connect_fluids(self):
        elements = [self.source, *self.turbomachines, self.sink]
        for i in range(1, len(elements) - 1):
            elements[i].source = elements[i - 1]
            elements[i].sink = elements[i + 1]

    def add_to_shaft(self, idxs: list[int]):
        for i in idxs:
            self.shaft.turbomachines.append(self.turbomachines[i])


def p_loss(omega, f_0=0.85, nu=100e-6, d_m=0.025, Cf=0.0025, rho=1.205, r_1=0, r_2=0.1):
    n = omega * 60 / 2 / np.pi
    M_0 = f_0 * 1e3 * (nu * n) ** (2 / 3) * d_m ** 3
    T_w = 0.2 * Cf * rho * omega ** 2 * np.pi * (r_2 ** 5 - r_1 ** 5)
    return (2 * M_0 + 4 * T_w) * omega


def eff_mech_max(omega, W_vac_max=597.6087961307821, eff_c=0.6, eff_t=0.65):
    return (W_vac_max - p_loss(omega) / eff_t) / (W_vac_max - p_loss(omega) * eff_c)


def dh0(omega, eff_t=0.65, eff_c=0.6, W_vac_max=597.6087961307821, mdot=0.0623):
    eff_mech = eff_mech_max(omega, W_vac_max, eff_c, eff_t)
    W_t_max = (eff_t * W_vac_max) / (1 - eff_mech * eff_c * eff_t)
    W_c_max = eff_mech * W_t_max
    return W_c_max / mdot, W_t_max / mdot
