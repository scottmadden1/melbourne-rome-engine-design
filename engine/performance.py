"""Whole-engine performance: a sized fan plus a core, combined through eq. (2)."""

from dataclasses import dataclass

from .core import CoreState, core_cycle
from .fan import Fan
from .gas import CONSTANT_CP
from .thrust import bare_thrust_eq2, bypass_exit_area


@dataclass(frozen=True)
class Performance:
    fan: Fan
    core: CoreState
    BPR: float
    A_b: float          # bypass nozzle exit area, m^2
    F_bare: float       # eq. (2), N
    F_net: float        # eq. (2) less nacelle and weight drag, N

    @property
    def velocity_ratio(self):
        return self.core.V_jc / self.fan.exit.V_jb


def fuel_air_ratio(fan, BPR, m_f):
    """Fuel-air ratio when fuel flow m_f is burnt in the core stream of this fan."""
    return m_f / (fan.m_air / (1 + BPR))


def evaluate(fan, BPR, f, gas=CONSTANT_CP, hpc_gas=CONSTANT_CP, tet_max=None):
    """Run the core behind a sized fan and sum both streams with eq. (2)."""
    core = core_cycle(fan.FPR, BPR, fan.M, fan.m_air, f, gas, hpc_gas, tet_max)
    ex = fan.exit
    A_b = bypass_exit_area(core.m_core * BPR, ex)
    F_bare = bare_thrust_eq2(core.m_core, core.V_jc, core.P_jc, core.A_c,
                             BPR, ex.V_jb, ex.P_jb, A_b, fan.V)
    return Performance(fan=fan, core=core, BPR=BPR, A_b=A_b, F_bare=F_bare,
                       F_net=F_bare - fan.D_install)
