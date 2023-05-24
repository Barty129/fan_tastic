import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fantastic.turbomachinery import p_loss, p_loss_analytic

MOM_INERTIA = 0.002456


def main():
    freewheel_df = pd.read_csv("../../prelims/freewheel-test.csv")
    w = freewheel_df["n (rpm)"].to_numpy() * 2 * np.pi / 60
    t_3000 = freewheel_df["t(n=3000) (s)"].to_numpy()

    w_fit = np.poly1d(np.polyfit(t_3000, w, 3))
    dw_dt_fit = np.polyder(w_fit)

    t_8000rpm = np.max((w_fit - 8000 * 2 * np.pi / 60).roots).real

    dw_dt_8000rpm = dw_dt_fit(t_8000rpm)
    torque = MOM_INERTIA * dw_dt_8000rpm
    P_loss = torque * 8000 * 2 * np.pi / 60

    print(f'dw/dt at 8000rpm = {dw_dt_8000rpm:.2f}rad/s^2')
    print(f'torque at 8000rpm: {torque:.2f}Nm')
    print(f'power loss at 8000rpm: {P_loss:.2f}W')

    plt.plot(t_3000, w_fit(t_3000), label='$\omega$ (fit)')
    plt.plot(t_3000, w, 'x', label='$\omega$ (data)')
    plt.legend()
    plt.ylabel('angular velocity, rad/s')
    plt.xlabel('time to reach 3000rpm, s')
    plt.show()


if __name__ == "__main__":
    w = np.pi / 30 * np.linspace(7000, 10000)
    plt.plot(w, p_loss(w), label='empirical')
    plt.plot(w, p_loss_analytic(w), label='analytic')
    plt.legend()
    plt.show()
