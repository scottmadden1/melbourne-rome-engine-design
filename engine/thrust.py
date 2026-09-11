"""Thrust equations and installation penalties (pure formulas, no cycle logic)."""

from . import constants as c


def installation_drag(m_air, d_f, L_D):
    """Nacelle drag plus the drag of carrying the engine's weight, N."""
    D_nacelle = c.NACELLE_DRAG_PER_KGS * m_air
    D_weight = (c.ENGINE_WEIGHT_REF * c.G / L_D) * (d_f / 3) ** c.ENGINE_WEIGHT_EXP
    return D_nacelle + D_weight


def bare_thrust_eq1(m_air, fan_area, fan_exit, V, acf):
    """Eq. (1): single-stream bare thrust, pressure term over ACF * fan area, N."""
    return m_air * (fan_exit.V_jb - V) + acf * (fan_exit.P_jb - c.P_A) * fan_area


def bypass_exit_area(m_bypass, fan_exit):
    """Bypass nozzle exit area from continuity at the exit state, m^2."""
    return m_bypass / (fan_exit.rho_jb * fan_exit.V_jb)


def bare_thrust_eq2(m_core, V_jc, P_jc, A_c, BPR, V_jb, P_jb, A_b, V):
    """Eq. (2): core stream plus bypass stream bare thrust, N."""
    return (m_core * (V_jc - V) + (P_jc - c.P_A) * A_c
            + m_core * BPR * (V_jb - V) + (P_jb - c.P_A) * A_b)


def propulsive_efficiency(F_bare, m_air, V_j, V):
    """Single-stream propulsive efficiency, V*F / (V*F + 0.5*m*(V_j - V)^2)."""
    return V * F_bare / (V * F_bare + 0.5 * m_air * (V_j - V) ** 2)
