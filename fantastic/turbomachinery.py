from dataclasses import dataclass


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
    def __init__(self):
        super().__init__()


class Turbine(Turbomachine):
    def __init__(self):
        super().__init__()


class Vacuum(Turbomachine):
    def __init__(self):
        super().__init__()


@dataclass
class Node:
    p: float
    T: float
    R: float

    @property
    def density(self):
        return self.p / self.R / self.T


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

        self.shafts = []

    def connect_fluids(self):
        elements = [self.source, *self.turbomachines, self.sink]
        for i in range(1, len(elements) - 1):
            elements[i].source = elements[i - 1]
            elements[i].sink = elements[i + 1]

    def add_shaft(self, idxs: list[int]):
        self.shafts.append(Shaft(0, [self.turbomachines[i] for i in idxs]))
