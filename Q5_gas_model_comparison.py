"""
Final design: enthalpy tables vs constant cp.

The final design from Q5_final_design_sweep.py, found with its default gas
model (GAS there), is held fixed (eta_o, FPR, BPR, cargo) and evaluated with
both gas models, which apply from the HPC to the core nozzle:

  Enthalpy tables  Cengel A-17 air, cp and gamma at each station's temperature,
                   and the HPC exit from the relative pressure column
  Constant cp      cp = 1004.5 J/(kg K) and gamma = 1.4 throughout

The intake, fan, booster and bypass stream use gamma = 1.4 in both, so the fan
sizing, m_air, m_core and f are identical and the two columns differ only in
the gas model. Every design constraint is checked under each model. The
comparison is printed and written as the LaTeX source of Table 3.
"""

import numpy as np

import Q5_final_design_sweep as sweep
from engine import CONSTANT_CP, ENTHALPY_TABLE
from engine.output import table_path

MODELS = (ENTHALPY_TABLE, CONSTANT_CP)
DEFAULT = sweep.GAS    # the model the sweep reports its final design with


def use(gas):
    """Run the sweep's cycle with this gas model from the HPC onwards."""
    sweep.GAS = sweep.HPC_GAS = gas


def final_design():
    """The sweep's search: first cargo with a matched design, then eta_o refined."""
    for cargo in sweep.CARGOS:
        best = sweep.highest_matched(cargo, sweep.ETA_SEARCH)
        if best:
            break
    else:
        return None
    eta_hi = min(best["eta_o"] + 0.005, sweep.ETAS[0])
    for eta_f in np.round(np.arange(eta_hi - sweep.FINE_STEP, best["eta_o"],
                                    -sweep.FINE_STEP), 4):
        b = sweep.best_matched(sweep.grid(eta_f, best["cargo"]))
        if b:
            return b
    return best


def verdict(r):
    if r["matched"]:
        return "matched"
    return "valid, Vjc/Vjb off" if r["valid"] else f"fails: {r['why']}"


def constraint_checks(r):
    """Each design constraint as (name, value, passes) for one evaluated design."""
    p = r["perf"]
    return [
        ("T_04 <= TET ceiling", f"{p.core.T_04:.0f} K", p.core.T_04 <= sweep.TET_MAX),
        ("|F_N - F_req| <= 4.5 kN", f"{r['gap'] / 1e3:+.2f} kN",
         abs(r["gap"]) <= sweep.THRUST_TOL),
        ("TOW <= MTOW", f"{r['TOW']:.1f} t", r["mtow_ok"]),
        ("V_jc > V", f"{p.core.V_jc:.1f} m/s", p.core.V_jc > p.fan.V),
        ("|V_jc/V_jb - 1| <= 0.01", f"{r['ratio']:.3f}",
         abs(r["ratio"] - 1) <= sweep.RATIO_TOL),
    ]


def signed(d, fmt):
    """Difference with an explicit sign, and a plain zero if it rounds to zero."""
    return f"{0:{fmt}}" if float(f"{abs(d):{fmt}}") == 0 else f"{d:+{fmt}}"


def station_rows(r):
    """(label, LaTeX label, unit, LaTeX unit, value, format) for one evaluated design."""
    p = r["perf"]
    core = p.core
    return [
        ("T_03", "$T_{03}$", "K", "K", core.T_03, ".1f"),
        ("T_04", "$T_{04}$", "K", "K", core.T_04, ".1f"),
        ("T_045", "$T_{045}$", "K", "K", core.T_045, ".1f"),
        ("p_045", "$p_{045}$", "kPa", "kPa", core.P_045 / 1e3, ".2f"),
        ("T_05", "$T_{05}$", "K", "K", core.T_05, ".1f"),
        ("p_05", "$p_{05}$", "kPa", "kPa", core.P_05 / 1e3, ".2f"),
        ("Core nozzle", "Core nozzle", "-", "--",
         "Choked" if core.choked else "Unchoked", None),
        ("V_jc", "$V_{jc}$", "m/s", "m/s", core.V_jc, ".1f"),
        ("V_jc/V_jb", "$V_{jc}/V_{jb}$", "-", "--", r["ratio"], ".3f"),
        ("A_c", "$A_c$", "m^2", "m$^2$", core.A_c, ".3f"),
        ("F_N", "$F_N$", "kN", "kN", p.F_net / 1e3, ".2f"),
        ("F_N - F_req", "$F_N - F_{req}$", "kN", "kN", r["gap"] / 1e3, ".2f"),
        ("sfc", "sfc", "g/(kN s)", "g/(kN\\,s)", r["m_f"] * 1e6 / p.F_net, ".2f"),
        ("eta_p", "$\\eta_p$", "-", "--", sweep.two_stream_propulsive_efficiency(p), ".4f"),
        ("eta_o", "$\\eta_o$", "-", "--", r["eta_ach"], ".4f"),
    ]


def compare_at(design):
    """The final design evaluated under each gas model, printed and written as Table 3."""
    results = {}
    for gas in MODELS:
        use(gas)
        results[gas] = sweep.evaluate(design["eta_o"], design["FPR"], design["BPR"],
                                      design["cargo"])
    use(DEFAULT)
    a, b = results[ENTHALPY_TABLE], results[CONSTANT_CP]
    fan = a["fan"]

    print(f"\n{'=' * 78}\nThe {DEFAULT.name} final design under both gas models"
          f"\n{'=' * 78}")
    print(f"Design point: M = {sweep.M}, eta_o = {design['eta_o']:.3f}, "
          f"FPR = {design['FPR']:.2f}, BPR = {design['BPR']:g}, "
          f"cargo = {design['cargo'] / 1e3:.0f} t")
    for gas, r in results.items():
        if "perf" not in r:
            print(f"{gas.name}: the cycle does not run at this point ({r['why']}).")
            return
    print(f"Common to both: TOW = {a['TOW']:.2f} t, F_req = {a['F_req'] / 1e3:.2f} kN, "
          f"d_f = {fan.d_f:.3f} m, m_air = {fan.m_air:.1f} kg/s,")
    print(f"                m_f = {a['m_f']:.4f} kg/s, f = {a['f']:.5f}, "
          f"V_jb = {fan.exit.V_jb:.1f} m/s, A_b = {a['perf'].A_b:.3f} m^2\n")

    print(f"{'Quantity':<13} {'Unit':<9} {ENTHALPY_TABLE.name:>15} {CONSTANT_CP.name:>12} "
          f"{'Difference':>11}")
    print("-" * 64)
    tex = []
    for (label, tex_label, unit, tex_unit, va, fmt), vb in zip(
            station_rows(a), (row[4] for row in station_rows(b))):
        if fmt is None:
            va_s, vb_s, d_s, d_tex = va, vb, "-", "--"
        else:
            va_s, vb_s = f"{va:{fmt}}", f"{vb:{fmt}}"
            d_s = signed(vb - va, fmt)
            d_tex = f"${d_s}$"
        print(f"{label:<13} {unit:<9} {va_s:>15} {vb_s:>12} {d_s:>11}")
        tex.append(f"{tex_label:<18}& {tex_unit:<11}& {va_s:<9}& {vb_s:<9}& {d_tex:<10}\\\\")

    print(f"\n{'Constraint':<26} {ENTHALPY_TABLE.name:>22} {CONSTANT_CP.name:>22}")
    print("-" * 72)
    for (name, va, oka), (_, vb, okb) in zip(constraint_checks(a), constraint_checks(b)):
        print(f"{name:<26} {va + (' PASS' if oka else ' FAIL'):>22} "
              f"{vb + (' PASS' if okb else ' FAIL'):>22}")
    print(f"{'Verdict':<26} {verdict(a):>22} {verdict(b):>22}")

    latex = (
        "\\begin{table}[H]\n\\centering\n\\begin{tabular}{l c c c c}\n\\toprule\n"
        "Quantity & Unit & Enthalpy tables & Constant $c_p$, $\\gamma = 1.4$ & Difference \\\\\n"
        "\\midrule\n" + "\n".join(tex) + "\n\\bottomrule\n\\end{tabular}\n"
        f"\\caption{{Final design ($M = {sweep.M:.2f}$, $\\eta_o = {design['eta_o']:.3f}$, "
        f"FPR $= {design['FPR']:.2f}$, BPR $= {design['BPR']:g}$, "
        f"{design['cargo'] / 1e3:.0f} t cargo) evaluated\nwith the ideal-gas enthalpy tables "
        "and with a constant $c_p$ and $\\gamma$ from the HPC onwards.}\n"
        "\\label{tab:final_cp_comparison}\n\\end{table}\n")
    path = table_path("table3_final_design_cp_comparison.tex")
    path.write_text(latex)
    print(f"\nLaTeX written to {path.relative_to(path.parent.parent)}")


def main():
    print(__doc__)
    use(DEFAULT)
    design = final_design()
    if design is None:
        print(f"The {DEFAULT.name} sweep finds no final design, so there is nothing "
              "to compare.")
        return
    compare_at(design)


if __name__ == "__main__":
    main()
