"""Core stream: booster, HPC, combustor, HPT, LPT and core nozzle."""

import math
from dataclasses import dataclass

from . import constants as c
from .fan import fan_exit_T, intake
from .gas import CONSTANT_CP, CycleError


def compressor_exit_T(T_in, pr, eta, gas):
    """Real compressor exit temperature: isentropic exit, then eta on enthalpy."""
    h_in = gas.h(T_in)
    h_out = h_in + (gas.h(gas.isentropic_T(T_in, pr)) - h_in) / eta
    return gas.T_from_h(h_out)


@dataclass(frozen=True)
class CoreState:
    FPR: float
    BPR: float
    f: float            # fuel-air ratio on the core flow
    m_core: float       # kg/s
    T_02: float
    P_02: float
    T_013: float
    P_013: float
    T_023: float
    P_023: float
    T_03: float
    P_03: float
    T_04: float
    P_04: float
    T_045: float
    P_045: float
    T_05: float         # LPT exit, also the nozzle inlet (no jet pipe loss)
    P_05: float
    P_c: float          # critical pressure at the core nozzle throat
    choked: bool
    T_9: float          # core nozzle exit static temperature
    P_jc: float         # core nozzle exit static pressure
    V_jc: float
    A_c: float          # core nozzle exit area, m^2


def core_cycle(FPR, BPR, M, m_air, f, gas=CONSTANT_CP, hpc_gas=CONSTANT_CP,
               tet_max=None):
    """
    Walk the core from the fan face to the core nozzle exit at fuel-air ratio f.

    gas      model for the combustor, turbine work balances and nozzle
    hpc_gas  model for the HPC's isentropic exit (the fan and booster always
             use gamma = 1.4)
    tet_max  optional turbine entry temperature ceiling on T_04, K

    Raises CycleError at the first station where the cycle stops being physical.
    """
    inlet = intake(M)
    T_02, P_02 = inlet.T_02, inlet.P_02
    pi_lpc = c.PI_FAN_LPC / FPR

    T_013, P_013 = fan_exit_T(T_02, FPR), P_02 * FPR
    T_023 = compressor_exit_T(T_013, pi_lpc, c.ETA_COMPRESSOR, CONSTANT_CP)
    P_023 = P_013 * pi_lpc
    T_03 = compressor_exit_T(T_023, c.PI_HPC, c.ETA_COMPRESSOR, hpc_gas) # hpc_gas can be both the constant_cp model or the enthalpy table model
    P_03 = P_023 * c.PI_HPC

    m_core = m_air / (1 + BPR)
    h_02, h_013, h_023, h_03 = (gas.h(T) for T in (T_02, T_013, T_023, T_03))

    # combustor: (1+f) h_04 = h_03 + f LCV, no pressure loss
    h_04 = (h_03 + f * c.LCV) / (1 + f)
    gas.check(h_04, "T_04")
    T_04, P_04 = gas.T_from_h(h_04), P_03
    if tet_max is not None and T_04 > tet_max:
        raise CycleError(f"T_04 = {T_04:.1f} K exceeds the {tet_max:.0f} K "
                         f"turbine entry temperature ceiling")

    # HPT drives the HPC
    h_045 = h_04 - (h_03 - h_023) / (1 + f)
    gas.check(h_045, "T_045")
    T_045 = gas.T_from_h(h_045)
    hpt = 1 - (1 - T_045 / T_04) / c.ETA_TURBINE
    if hpt <= 0:
        raise CycleError("HPT expansion term non-positive: HPT cannot drive the HPC")
    g = gas.gamma((T_04 + T_045) / 2)
    P_045 = P_04 * hpt ** (g / (g - 1))

    # LPT drives the fan on the full flow and the booster on the core flow
    h_05 = h_045 - ((1 + BPR) * (h_013 - h_02) + (h_023 - h_013)) / (1 + f)
    gas.check(h_05, "T_05")
    T_05 = gas.T_from_h(h_05)
    lpt = 1 - (T_045 - T_05) / (c.ETA_TURBINE * T_045)
    if lpt <= 0:
        raise CycleError("LPT expansion term non-positive: LPT cannot drive the fan")
    g = gas.gamma((T_045 + T_05) / 2)
    P_05 = P_045 * lpt ** (g / (g - 1))

    # core nozzle, fed directly from the LPT exit
    g, cp = gas.gamma(T_05), gas.cp(T_05)
    k = (g - 1) / g
    P_c = P_05 / ((g + 1) / 2) ** (g / (g - 1))
    choked = P_c > c.P_A
    if choked:
        T_9 = T_05 / ((g + 1) / 2)
        V_jc = math.sqrt(g * c.R * T_9)
        P_jc = P_c
    else:
        if P_05 < c.P_A:
            raise CycleError(f"P_05 = {P_05/1e3:.2f} kPa is below ambient "
                             f"({c.P_A/1e3:.2f} kPa): nozzle cannot expand")
        T_9 = T_05 * (c.P_A / P_05) ** k
        V_jc = math.sqrt(2 * cp * T_05 * (1 - (c.P_A / P_05) ** k))
        P_jc = c.P_A
    A_c = m_core / (P_jc / (c.R * T_9) * V_jc)

    return CoreState(FPR=FPR, BPR=BPR, f=f, m_core=m_core,
                     T_02=T_02, P_02=P_02, T_013=T_013, P_013=P_013,
                     T_023=T_023, P_023=P_023, T_03=T_03, P_03=P_03,
                     T_04=T_04, P_04=P_04, T_045=T_045, P_045=P_045,
                     T_05=T_05, P_05=P_05, P_c=P_c, choked=choked, T_9=T_9,
                     P_jc=P_jc, V_jc=V_jc, A_c=A_c)
