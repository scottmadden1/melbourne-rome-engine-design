"""
Q4.3 Fan design.

For each fan pressure ratio, sizes the fan diameter so that eq. (1) with
ACF = 1, less nacelle and engine weight drag, meets the thrust per engine.
"""

import matplotlib.pyplot as plt
import numpy as np

from engine import mission as ms
from engine.fan import size_fan
from engine.output import figure_path
from engine.thrust import propulsive_efficiency

FPR_RANGE = np.round(np.arange(1.45, 1.605, 0.01), 2)
ACF = 1.0


def main():
    M = ms.max_cruise_mach()
    F_req = ms.required_thrust(M)
    fans = [size_fan(FPR, M, F_req, acf=ACF) for FPR in FPR_RANGE]
    eta_p = [propulsive_efficiency(f.F_bare_eq1, f.m_air, f.exit.V_jb, f.V)
             for f in fans]

    ex0 = fans[0].exit
    print(f"M = {M:.2f}, V = {fans[0].V:.2f} m/s, required thrust = {F_req/1e3:.2f} kN per engine")
    print(f"Ambient pressure = 21.66 kPa, ACF = {ACF}\n")
    print(f"{'FPR':>5} {'P_c kPa':>8} {'nozzle':>9} {'V_jb m/s':>9} {'d_f m':>7} "
          f"{'m_air kg/s':>11} {'eta_p':>7}")
    for fan, ep in zip(fans, eta_p):
        ex = fan.exit
        print(f"{fan.FPR:5.2f} {ex.P_c/1e3:8.2f} {'choked' if ex.choked else 'unchoked':>9} "
              f"{ex.V_jb:9.2f} {fan.d_f:7.4f} {fan.m_air:11.2f} {ep:7.4f}")

    d_f = [f.d_f for f in fans]
    m_air = [f.m_air for f in fans]
    V_jb = [f.exit.V_jb for f in fans]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    ax1.plot(FPR_RANGE, d_f, marker="o")
    ax1.set(xlabel="Fan pressure ratio, FPR", ylabel="Required fan diameter, $d_f$ (m)",
            title="Fan diameter required to meet $F_{req}$")
    ax1.grid(True)
    ax2.plot(FPR_RANGE, m_air, marker="o", color="tab:orange")
    ax2.set(xlabel="Fan pressure ratio, FPR", ylabel="Mass flow rate of air, $\\dot{m}_{air}$ (kg/s)",
            title="Mass flow rate of air at the solved $d_f$")
    ax2.grid(True)
    fig.tight_layout()
    fig.savefig(figure_path("fan_diameter_mass_flow_rate.png"), dpi=150)

    fig2, (ax3, ax4) = plt.subplots(1, 2, figsize=(12, 5))
    ax3.plot(FPR_RANGE, V_jb, marker="o", color="tab:green")
    ax3.set(xlabel="Fan pressure ratio, FPR", ylabel="Bypass jet velocity, $V_{jb}$ (m/s)",
            title="Bypass jet velocity at the solved $d_f$")
    ax3.grid(True)
    ax4.plot(FPR_RANGE, eta_p, marker="o", color="tab:red")
    ax4.set(xlabel="Fan pressure ratio, FPR", ylabel="Propulsive efficiency, $\\eta_P$",
            title="Propulsive efficiency (assuming $f=0$, $V_{jc} = V_{jb}$)")
    ax4.grid(True)
    fig2.tight_layout()
    fig2.savefig(figure_path("fan_vjb_propulsive_efficiency.png"), dpi=150)
    plt.show()


if __name__ == "__main__":
    main()
