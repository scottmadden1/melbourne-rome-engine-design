"""
Shared propulsion library for the Assignment 1 report.

Each module holds one piece of the physics, written once:

    constants  handout values and unit conversions, all SI
    mission    L/D, cruise distance, Breguet fuel burn, required thrust
    gas        ConstantCp and EnthalpyTable gas models (same interface)
    fan        intake, fan exit and bypass nozzle, eq. (1) fan sizing, ACF
    core       booster, HPC, combustor, HPT, LPT and core nozzle
    thrust     eq. (1) and (2) thrust formulas and installation penalties
    performance  sized fan plus core, combined through eq. (2)
    output     figure and table output paths

Nothing in this package prints or plots; the Q*.py scripts do that.
"""

from .gas import CONSTANT_CP, ENTHALPY_TABLE, CycleError

__all__ = ["CONSTANT_CP", "ENTHALPY_TABLE", "CycleError"]

