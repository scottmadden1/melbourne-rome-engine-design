"""
Q4.1 Choosing the flight Mach number.

Sweeps cruise Mach number and finds the fastest cruise whose take-off weight,
from the exact Breguet range equation with a 30 kn headwind, still fits
within the maximum take-off weight.
"""

import numpy as np

from engine import constants as c
from engine import mission as ms

MACH_RANGE = np.round(np.arange(0.75, 0.875, 0.01), 2)


def main():
    m_takeoff = c.MTOW * c.F_TAKEOFF
    m_landing = c.MTOW * c.F_LANDING
    fuel_capacity = c.MTOW - c.M_EMPTY - ms.payload() - m_takeoff - m_landing
    print(f"Mass of fuel available for cruise: {fuel_capacity/1e3:.2f} t")
    print(f"Mass of fuel burnt in take-off: {m_takeoff/1e3:.2f} t")
    print(f"Mass of fuel burnt in landing: {m_landing/1e3:.2f} t")
    print(f"Distance in cruise: {ms.cruise_distance_nm():.2f} nm\n")

    for M in MACH_RANGE:
        tow, m_fb = ms.actual_tow(M)
        status = "OK" if tow <= c.MTOW else "EXCEEDS CAP"
        print(f"M={M:.2f}: m_fb = {m_fb/1e3:.2f} t, ACTUAL_TOW = {tow/1e3:.2f} t  [{status}]")

    M = ms.max_cruise_mach(MACH_RANGE)
    tow, _ = ms.actual_tow(M)
    V = ms.flight_speed(M)
    print(f"\nMax allowable cruise Mach: {M:.2f}")
    print(f"Actual take-off weight at this Mach (with headwind): {tow/1e3:.2f} t "
          f"(cap: {c.MTOW/1e3:.0f} t)")
    print(f"the speed of sound at 37000 ft: {c.A_A:.2f} m/s")
    print(f"the ideal cruise speed: {V:.2f} m/s")
    print(f"the ideal cruise groundspeed: {V - c.HEADWIND_KN * c.KNOT_TO_MPS:.2f} m/s")
    print(f"Flight time (cruise) at M={M:.2f}, with headwind: "
          f"{ms.cruise_time(M, headwind=True)/3600:.2f} hr")
    print(f"Flight time (cruise) at M={M:.2f}, with no headwind: "
          f"{ms.cruise_time(M, headwind=False)/3600:.2f} hr")

    fuel_with = ms.cruise_fuel_burn(M, headwind=True)
    fuel_without = ms.cruise_fuel_burn(M, headwind=False)
    print(f"Fuel required with headwind at M={M:.2f}: {fuel_with/1e3:.2f} t")
    print(f"Fuel required without headwind at M={M:.2f}: {fuel_without/1e3:.2f} t")
    print(f"Extra fuel required due to headwind at M={M:.2f}: "
          f"{(fuel_with - fuel_without)/1e3:.2f} t")


if __name__ == "__main__":
    main()
