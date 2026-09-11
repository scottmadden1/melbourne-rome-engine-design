"""
Gas property models.

Both models expose the same methods, so the core cycle is written once in
enthalpy form and either model can be passed in:

    h(T)                    specific enthalpy, J/kg
    T_from_h(h)             inverse of h(T), K
    cp(T)                   J/(kg K)
    gamma(T)
    isentropic_T(T1, pr)    exit temperature of an isentropic process with
                            pressure ratio pr = p2/p1
    check(h, label)         raise CycleError if h is outside the model's range

For ConstantCp, h = cp*T and gamma is fixed, which reduces every enthalpy
balance to the familiar constant-cp temperature form.
"""

import math

import numpy as np

from . import constants as c


class CycleError(ValueError):
    """The cycle does not close at this operating point."""


class ConstantCp:
    name = "Constant cp"

    def __init__(self, cp=c.CP, gamma=c.GAMMA):
        self._cp, self._gamma = cp, gamma

    def h(self, T):
        return self._cp * T

    def T_from_h(self, h):
        return h / self._cp

    def cp(self, T):
        return self._cp

    def gamma(self, T):
        return self._gamma

    def isentropic_T(self, T1, pr):
        return T1 * pr ** ((self._gamma - 1) / self._gamma)

    def check(self, h, label):
        pass


class EnthalpyTable:
    """
    Ideal-gas air, Cengel Table A-17. cp is the local slope of h(T) and
    gamma = cp/(cp - R). Isentropic processes use the relative pressure
    column, p2/p1 = Pr(T2)/Pr(T1), interpolated in ln(Pr).
    """
    name = "Enthalpy tables"

    T_TABLE = np.array([200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100,
                        1200, 1300, 1400, 1500, 1600, 1700, 1800], dtype=float)
    H_TABLE_KJ = np.array([199.97, 300.19, 400.98, 503.02, 607.32, 713.27,
                           822.20, 933.33, 1046.04, 1161.07, 1277.79, 1395.97,
                           1515.42, 1635.97, 1757.57, 1880.10, 2003.30])
    PR_TABLE = np.array([0.3363, 1.3860, 3.806, 8.411, 16.28, 28.80, 47.75,
                         75.29, 114.0, 167.1, 238.0, 330.9, 450.5, 601.9,
                         791.2, 1025.0, 1310.0])

    def __init__(self):
        self._h = self.H_TABLE_KJ * 1e3
        cp_kj = np.gradient(self.H_TABLE_KJ, self.T_TABLE)
        self._cp = cp_kj * 1e3
        self._gamma = cp_kj / (cp_kj - c.R / 1e3)
        self._ln_pr = np.log(self.PR_TABLE)

    def h(self, T):
        return float(np.interp(T, self.T_TABLE, self._h))

    def T_from_h(self, h):
        return float(np.interp(h, self._h, self.T_TABLE))

    def cp(self, T):
        return float(np.interp(T, self.T_TABLE, self._cp))

    def gamma(self, T):
        return float(np.interp(T, self.T_TABLE, self._gamma))

    def isentropic_T(self, T1, pr):
        ln_pr2 = np.interp(T1, self.T_TABLE, self._ln_pr) + math.log(pr)
        return float(np.interp(ln_pr2, self._ln_pr, self.T_TABLE))

    def check(self, h, label):
        """np.interp silently clamps outside the table, so catch it here."""
        if h > self._h[-1]:
            raise CycleError(f"{label} above the enthalpy table ceiling "
                             f"({self.T_TABLE[-1]:.0f} K)")
        if h < self._h[0]:
            raise CycleError(f"{label} below the enthalpy table floor "
                             f"({self.T_TABLE[0]:.0f} K)")


CONSTANT_CP = ConstantCp()
ENTHALPY_TABLE = EnthalpyTable()
