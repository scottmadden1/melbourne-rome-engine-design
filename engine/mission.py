"""Mission model: L/D, cruise distance, Breguet fuel burn and required thrust."""

import math

import numpy as np

from . import constants as c


def lift_to_drag(M):
    """Handout L/D correlation."""
    return 22.65 * (1.8 - M) ** 2


def flight_speed(M):
    """True airspeed at 37,000 ft, m/s."""
    return M * c.A_A


def _unit_vector(lat_deg, lon_deg):
    lat, lon = math.radians(lat_deg), math.radians(lon_deg)
    return np.array([math.cos(lat) * math.cos(lon),
                     math.cos(lat) * math.sin(lon),
                     math.sin(lat)])


def great_circle_nm(a=c.MELBOURNE, b=c.ROME):
    """Great-circle distance between two (lat, lon) points, nautical miles."""
    angle = math.acos(np.dot(_unit_vector(*a), _unit_vector(*b)))
    return c.KM_TO_NM * c.R_EARTH_KM * angle


def cruise_distance_nm():
    """Great-circle distance less the take-off and landing legs, nautical miles."""
    return great_circle_nm() - c.DIST_TAKEOFF_NM - c.DIST_LANDING_NM


def payload(cargo=c.MAX_CARGO):
    """Passengers plus cargo, kg."""
    return c.MAX_PASSENGERS * c.PAX_MASS + cargo


def cruise_time(M, headwind=True):
    """Time in cruise, s."""
    V = flight_speed(M)
    V_ground = V - c.HEADWIND_KN * c.KNOT_TO_MPS if headwind else V
    return cruise_distance_nm() * c.NM_TO_M / V_ground


def cruise_fuel_burn(M, eta_o=c.ETA_O_MISSION, cargo=c.MAX_CARGO, headwind=True):
    """
    Exact Breguet cruise fuel burn, kg:

        m_fb = (m_e + m_pl + m_landing) * (exp(R/H) - 1),   H = eta_o*LCV*(L/D)/g

    where R is the distance flown through the air mass, so a headwind
    lengthens it. Returns infinity when the mission is unflyable (R/H > 50).
    """
    R_air = flight_speed(M) * cruise_time(M, headwind)
    H = eta_o * c.LCV * lift_to_drag(M) / c.G
    ratio = R_air / H if H > 0 else float("inf")
    if ratio > 50:
        return float("inf")
    m_landing = c.MTOW * c.F_LANDING
    return (c.M_EMPTY + payload(cargo) + m_landing) * (math.exp(ratio) - 1)


def actual_tow(M, eta_o=c.ETA_O_MISSION, cargo=c.MAX_CARGO, headwind=True):
    """Take-off weight needed to fly the mission, and its cruise fuel, kg."""
    m_fb = cruise_fuel_burn(M, eta_o, cargo, headwind)
    tow = (c.M_EMPTY + payload(cargo) + m_fb
           + c.MTOW * c.F_LANDING + c.MTOW * c.F_TAKEOFF)
    return tow, m_fb


def required_thrust(M, tow=c.MTOW):
    """Cruise thrust per engine = drag at the start of cruise / N_engines, N."""
    m_start_cruise = tow * (1 - c.F_TAKEOFF)
    return m_start_cruise * c.G / (lift_to_drag(M) * c.N_ENGINES)


def cruise_fuel_flow(M, F_req, eta_o=c.ETA_O_MISSION):
    """Fuel flow per engine that delivers F_req at overall efficiency eta_o, kg/s."""
    return F_req * flight_speed(M) / (eta_o * c.LCV)


def max_cruise_mach(mach_range=None, eta_o=c.ETA_O_MISSION, cargo=c.MAX_CARGO):
    """Highest Mach in the range whose take-off weight (with headwind) fits MTOW."""
    if mach_range is None:
        mach_range = np.round(np.arange(0.75, 0.875, 0.01), 2)
    feasible = [M for M in mach_range if actual_tow(M, eta_o, cargo)[0] <= c.MTOW]
    return float(max(feasible)) if feasible else None
