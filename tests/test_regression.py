"""
Regression check: the refactored library must reproduce the numbers the
original scripts printed before the refactor. The original scripts typed in
58.36 kN and 0.7488 kg/s by hand, so those are used here; the Q scripts now
compute them from the mission model instead.

Run from the "final code" folder:   python tests/test_regression.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from engine import CONSTANT_CP, ENTHALPY_TABLE              # noqa: E402
from engine import mission as ms                            # noqa: E402
from engine.fan import size_fan                             # noqa: E402
from engine.performance import evaluate, fuel_air_ratio     # noqa: E402
import Q5_sensitivity_studies as q5                         # noqa: E402

OLD_F_REQ, OLD_M_F, M = 58.36e3, 0.7488, 0.86
failures = []


def check(name, got, expected, decimals):
    ok = round(got, decimals) == round(expected, decimals)
    if not ok:
        failures.append(name)
    print(f"{'ok  ' if ok else 'FAIL'} {name:<34} got {got:.{decimals}f}  expected {expected:.{decimals}f}")


# Q4.1
check("Q1 cruise distance (nm)", ms.cruise_distance_nm(), 8293.37, 2)
check("Q1 max cruise Mach", ms.max_cruise_mach(), 0.86, 2)
tow, m_fb = ms.actual_tow(0.86)
check("Q1 TOW at M=0.86 (t)", tow / 1e3, 242.32, 2)
check("Q1 cruise fuel at M=0.86 (t)", m_fb / 1e3, 79.12, 2)
check("Q1 headwind extra fuel (t)",
      (ms.cruise_fuel_burn(0.86) - ms.cruise_fuel_burn(0.86, headwind=False)) / 1e3, 5.78, 2)
check("Q1 TOW at M=0.87 (t)", ms.actual_tow(0.87)[0] / 1e3, 244.33, 2)

# Q4.2
F_req = ms.required_thrust(M)
m_f = ms.cruise_fuel_flow(M, F_req)
check("Q2 thrust per engine (kN)", F_req / 1e3, 58.36, 2)
check("Q2 fuel flow per engine (kg/s)", m_f, 0.7488, 4)
check("Q2 SFC (g/kN/s)", m_f * 1e3 / (F_req / 1e3), 12.829044, 6)

# Q4.3 (fan at FPR 1.5 and 1.6, old hard-coded thrust)
fan15 = size_fan(1.5, M, OLD_F_REQ)
check("Q3 d_f at FPR 1.5 (m)", fan15.d_f, 2.8282, 4)
check("Q3 m_air at FPR 1.5 (kg/s)", fan15.m_air, 555.21, 2)
check("Q3 V_jb at FPR 1.5 (m/s)", fan15.exit.V_jb, 307.66, 2)
check("Q3 bypass P_c at FPR 1.45 (Pa)", size_fan(1.45, M, OLD_F_REQ).exit.P_c, 26889.50, 2)

# Q4.4 design point, both gas models (old hard-coded thrust and fuel flow)
fan16 = size_fan(1.6, M, OLD_F_REQ)
check("Q4 fan face area (m^2)", fan16.area, 5.0951, 4)
expected = {
    CONSTANT_CP: dict(T_03=875.67, T_04=1493.95, T_045=949.85, P_045=341.74, T_05=563.20,
                      P_05=41.55, V_jc=434.29, A_c=0.7071, A_b=2.9947, F_net=47.92),
    ENTHALPY_TABLE: dict(T_03=875.67, T_04=1409.57, T_045=917.50, P_045=284.57, T_05=559.83,
                         P_05=33.34, V_jc=361.60, A_c=0.9117, A_b=2.9947, F_net=44.08),
}
for gas, exp in expected.items():
    p = evaluate(fan16, 8, fuel_air_ratio(fan16, 8, OLD_M_F), gas)
    got = dict(T_03=p.core.T_03, T_04=p.core.T_04, T_045=p.core.T_045,
               P_045=p.core.P_045 / 1e3, T_05=p.core.T_05, P_05=p.core.P_05 / 1e3,
               V_jc=p.core.V_jc, A_c=p.core.A_c, A_b=p.A_b, F_net=p.F_net / 1e3)
    for k, v in exp.items():
        dec = 4 if k in ("A_c", "A_b") else 2
        check(f"Q4 {gas.name[:8]} {k}", got[k], v, dec)
    check(f"Q4 {gas.name[:8]} choked", float(p.core.choked), float(gas is CONSTANT_CP), 0)

# Q4.5 baseline operating point
r = q5.evaluate(0.75, 0.35, 21e3, 1.5, 8)
check("Q5 baseline TOW (t)", r["TOW"], 248.50, 2)
check("Q5 baseline F_req (kN)", r["F_req"] / 1e3, 47.84, 2)
check("Q5 baseline f_sup", r["f_sup"], 0.01176, 5)
check("Q5 baseline f_req", r["f_req"], 0.01165, 5)
check("Q5 baseline f_min", r["f_min"], 0.01087, 5)
check("Q5 baseline f_max", r["f_max"], 0.01914, 5)
check("Q5 baseline F_N(sup) (kN)", r["F_N_sup"] / 1e3, 44.26, 2)
check("Q5 baseline T_04 (K)", r["T_04_sup"], 1248, 0)
r = q5.evaluate(0.75, 0.35, 21e3, 1.5, 12)
check("Q5 BPR 12 gap (kN)", r["gap_sup"] / 1e3, 0.04, 2)

print(f"\n{'ALL PASSED' if not failures else f'{len(failures)} FAILED: ' + ', '.join(failures)}")
sys.exit(1 if failures else 0)
