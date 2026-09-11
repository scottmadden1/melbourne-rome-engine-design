"""
Q4.2 Fuel burn metrics.

Cruise drag at the start of cruise sets the thrust per engine; the minimum
overall efficiency then sets the fuel flow and specific fuel consumption.
These are the thrust and fuel-flow targets every later question uses.
"""

from engine import mission as ms


def main():
    M = ms.max_cruise_mach()
    F_req = ms.required_thrust(M)           # per engine, N
    m_f = ms.cruise_fuel_flow(M, F_req)     # per engine, kg/s
    sfc = m_f * 1e3 / (F_req / 1e3)         # g/(kN s)

    print(f"Cruise Mach number (from Q4.1): {M:.2f}")
    print(f"Fuel burn rate in cruise: {m_f:.4f} kg/s per engine")
    print(f"Drag in cruise: {F_req/1e3:.2f} kN per engine")
    print(f"Specific fuel consumption in cruise: {sfc:.6f} g/kN/s")


if __name__ == "__main__":
    main()
