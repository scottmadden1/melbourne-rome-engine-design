"""Intake, fan and bypass nozzle, eq. (1) fan sizing and the area correction factor."""

import math
from dataclasses import dataclass

from scipy.optimize import brentq

from . import constants as c
from .mission import flight_speed, lift_to_drag
from .thrust import bare_thrust_eq1, installation_drag

_K = (c.GAMMA - 1) / c.GAMMA


@dataclass(frozen=True)
class Intake:
    M: float
    V: float            # m/s
    T_02: float         # K
    P_02: float         # Pa


def intake(M):
    """Stagnation conditions at the fan face (station 2)."""
    ram = 1 + c.ETA_INTAKE * (c.GAMMA - 1) / 2 * M ** 2
    return Intake(M=M, V=flight_speed(M), T_02=c.T_A * ram,
                  P_02=c.P_A * ram ** (c.GAMMA / (c.GAMMA - 1)))


def fan_exit_T(T_02, FPR):
    """Fan exit stagnation temperature with gamma = 1.4, K."""
    return T_02 * (1 + (FPR ** _K - 1) / c.ETA_FAN)


@dataclass(frozen=True)
class FanExit:
    FPR: float
    T_013: float        # K
    P_013: float        # Pa
    P_c: float          # critical pressure at the bypass throat, Pa
    choked: bool
    T_19: float         # bypass exit static temperature, K
    P_jb: float         # bypass exit static pressure, Pa
    V_jb: float         # bypass jet velocity, m/s
    rho_jb: float       # bypass exit density, kg/m^3


def fan_exit(FPR, M):
    """Fan exit (013) and bypass nozzle exit (19) state, choked or fully expanded."""
    inlet = intake(M)
    P_013 = inlet.P_02 * FPR
    T_013 = fan_exit_T(inlet.T_02, FPR)
    P_c = P_013 / ((c.GAMMA + 1) / 2) ** (c.GAMMA / (c.GAMMA - 1))
    choked = P_c > c.P_A
    if choked:
        T_19 = T_013 / ((c.GAMMA + 1) / 2)
        V_jb = math.sqrt(c.GAMMA * c.R * T_19)
        P_jb = P_c
    else:
        T_19 = T_013 * (c.P_A / P_013) ** _K
        V_jb = math.sqrt(2 * c.CP * T_013 * (1 - (c.P_A / P_013) ** _K))
        P_jb = c.P_A
    return FanExit(FPR=FPR, T_013=T_013, P_013=P_013, P_c=P_c, choked=choked,
                   T_19=T_19, P_jb=P_jb, V_jb=V_jb, rho_jb=P_jb / (c.R * T_19))


def matched_acf(FPR, BPR, M):
    """
    ACF from matching eqs. (1) and (2) with an unchoked core and V_jc = V_jb:

        ACF = A_b / A_fan = [BPR/(1+BPR)] * (rho_a*V) / (rho_jb*V_jb)

    It depends only on FPR, BPR and the flight condition, not on d_f.
    """
    ex = fan_exit(FPR, M)
    return (BPR / (1 + BPR)) * (c.RHO_A * flight_speed(M)) / (ex.rho_jb * ex.V_jb)


@dataclass(frozen=True)
class Fan:
    FPR: float
    M: float
    V: float
    L_D: float
    F_req: float        # installed thrust the fan was sized to, N
    acf: float
    d_f: float          # m
    area: float         # fan face area, m^2
    m_air: float        # kg/s
    exit: FanExit
    F_bare_eq1: float   # eq. (1) bare thrust at the solved size, N
    D_install: float    # nacelle plus weight drag, N


def _build(FPR, M, F_req, acf, ex, d_f):
    V, L_D = flight_speed(M), lift_to_drag(M)
    area = math.pi * d_f ** 2 / 4
    m_air = c.RHO_A * V * area
    return Fan(FPR=FPR, M=M, V=V, L_D=L_D, F_req=F_req, acf=acf, d_f=d_f,
               area=area, m_air=m_air, exit=ex,
               F_bare_eq1=bare_thrust_eq1(m_air, area, ex, V, acf),
               D_install=installation_drag(m_air, d_f, L_D))


def size_fan(FPR, M, F_req, acf=1.0, bracket=(0.5, 5.0)):
    """Solve eq. (1) less installation drag = F_req for the fan diameter."""
    ex = fan_exit(FPR, M)

    def residual(d_f):
        fan = _build(FPR, M, F_req, acf, ex, d_f)
        return fan.F_bare_eq1 - fan.D_install - F_req

    return _build(FPR, M, F_req, acf, ex, brentq(residual, *bracket))
