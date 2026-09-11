"""
Final design sweep: overall efficiency stepped down at a fixed cruise Mach.

The Q5 sensitivity studies showed that lowering eta_o is the only input that
gives the engine enough fuel to run its cycle, and that the fuel deficit is
smallest around M = 0.75. So M is fixed at 0.75 and eta_o is stepped DOWN from
the handout's 0.46. At every eta_o the full mission is re-run (TOW -> F_req),
and the Q4 core design is repeated over an FPR/BPR grid, with the fan sized by
the matched ACF. The core uses a single constant cp and gamma = 1.4 (GAS
below); Q5_gas_model_comparison.py reruns the sweep with the enthalpy tables.

A design is valid when all five hold:
  1. the cycle runs            -- p_05 >= p_a (and inside the table, if used)
  2. T_04 <= TET_MAX           -- turbine entry temperature ceiling
  3. |F_N - F_req| <= 4.5 kN   -- the handout's thrust design tolerance
  4. TOW <= MTOW
  5. V_jc > V                  -- the core must produce forward thrust

and it is matched when, in addition, |V_jc/V_jb - 1| <= RATIO_TOL, the
handout's core design aim and the assumption behind the matched ACF. Among the
matched designs at one eta_o, the one with the highest achieved eta_o is taken.

Lowering eta_o burns more fuel, which pushes TOW over MTOW before a matched
design appears at the full payload. The compromise made here is cargo: it is
reduced 1 t at a time until a matched design exists, and the final design is
the matched design at the highest eta_o for that cargo, refined to 0.001.
A final check holds that cargo and asks whether any faster cruise, up to
Q1's M_max, still gives a matched design.
"""

import math

import numpy as np

from engine import CONSTANT_CP, CycleError
from engine import constants as c
from engine import mission as ms
from engine.fan import matched_acf, size_fan
from engine.output import table_path
from engine.performance import evaluate as engine_performance, fuel_air_ratio

M = 0.75
GAS = HPC_GAS = CONSTANT_CP    # engine.ENTHALPY_TABLE for variable cp
TET_MAX = 1500.0       # K, the same ceiling as Q5 (not in the handout)
THRUST_TOL = 4.5e3     # N, the handout's design tolerance
RATIO_TOL = 0.01       # |V_jc/V_jb - 1| for a matched core

FPR_SCAN = np.round(np.arange(1.45, 1.6001, 0.01), 2)
BPR_SCAN = np.round(np.arange(4.0, 14.01, 0.5), 1)
ETAS = np.round(np.arange(0.46, 0.1999, -0.01), 3)        # full-payload table
ETA_SEARCH = np.round(np.arange(0.46, 0.2999, -0.005), 3)  # per-cargo search
CARGOS = np.arange(c.MAX_CARGO, -1.0, -1e3)                 # kg, 21 t down to 0
FINE_STEP = 0.001


def payload_t(cargo):
    return (c.MAX_PASSENGERS * c.PAX_MASS + cargo) / 1e3


def evaluate(eta_o, FPR, BPR, cargo=c.MAX_CARGO, mach=None):
    """Mission, fan and core at one (eta_o, FPR, BPR, cargo) point, with the verdict."""
    mach = M if mach is None else mach
    TOW, m_fb = ms.actual_tow(mach, eta_o, cargo)
    out = dict(M=mach, eta_o=eta_o, FPR=FPR, BPR=BPR, cargo=cargo, TOW=TOW / 1e3,
               m_fb=m_fb / 1e3, mtow_ok=TOW <= c.MTOW, valid=False, matched=False)
    if not math.isfinite(TOW):
        return {**out, "why": "unflyable"}
    F_req = ms.required_thrust(mach, TOW)
    try:
        fan = size_fan(FPR, mach, F_req, acf=matched_acf(FPR, BPR, mach), bracket=(0.3, 9.0))
    except ValueError:
        return {**out, "why": "no fan d_f"}
    m_f = ms.cruise_fuel_flow(mach, F_req, eta_o)
    f = fuel_air_ratio(fan, BPR, m_f)
    out.update(F_req=F_req, fan=fan, m_f=m_f, f=f)
    try:
        p = engine_performance(fan, BPR, f, GAS, HPC_GAS)
    except CycleError:
        return {**out, "why": "cycle"}

    gap = p.F_net - F_req
    out.update(perf=p, gap=gap, T_04=p.core.T_04, ratio=p.velocity_ratio,
               eta_ach=p.F_net * fan.V / (m_f * c.LCV))
    checks = (("MTOW", out["mtow_ok"]), ("TET", p.core.T_04 <= TET_MAX),
              ("thrust", abs(gap) <= THRUST_TOL), ("Vjc<V", p.core.V_jc > fan.V))
    fails = [name for name, ok in checks if not ok]
    out.update(valid=not fails, why=fails[0] if fails else "")
    out["matched"] = out["valid"] and abs(out["ratio"] - 1) <= RATIO_TOL
    return out


def grid(eta_o, cargo=c.MAX_CARGO, mach=None):
    return [evaluate(eta_o, FPR, BPR, cargo, mach) for FPR in FPR_SCAN for BPR in BPR_SCAN]


def best_matched(rows):
    """Highest achieved eta_o among the matched designs, or None."""
    matched = [r for r in rows if r["matched"]]
    return max(matched, key=lambda r: r["eta_ach"]) if matched else None


def highest_matched(cargo, etas, mach=None):
    """Best matched design at the first (highest) eta_o that has one, or None."""
    for eta_o in etas:
        best = best_matched(grid(eta_o, cargo, mach))
        if best:
            return best
    return None


def two_stream_propulsive_efficiency(p):
    """V*F / (V*F + kinetic energy left in both jets), with F the eq. (2) bare thrust."""
    V, core = p.fan.V, p.core
    wasted = (0.5 * core.m_core * (1 + core.f) * (core.V_jc - V) ** 2
              + 0.5 * core.m_core * p.BPR * (p.fan.exit.V_jb - V) ** 2)
    return V * p.F_bare / (V * p.F_bare + wasted)


def spec_rows(r):
    """
    The handout's Table 1 for one design, as (label, LaTeX label, value).
    None marks the start of the per-engine rows.
    """
    p, fan = r["perf"], r["fan"]
    return [
        ("Max Mach number, M_max", "Max Mach number ($M_{max}$)",
         f"{ms.max_cruise_mach():.2f}"),
        ("Design M", "Design $M$ (if different to $M_{max}$)", f"{M:.2f}"),
        ("Flight speed (m/s)", "Flight speed (m/s)", f"{fan.V:.1f}"),
        None,
        ("Fuel flow rate, m_f (kg/s)", "Fuel flow rate, $\\dot{m}_f$ (kg/s)",
         f"{r['m_f']:.4f}"),
        ("Net thrust F_N (kN)", "Net thrust $F_N$ (kN)", f"{p.F_net / 1e3:.2f}"),
        ("sfc (g/(kN s))", "$sfc$ (g/(kN\\,s))", f"{r['m_f'] * 1e6 / p.F_net:.2f}"),
        ("FPR", "FPR", f"{r['FPR']:.2f}"),
        ("Bypass jet velocity, V_jb (m/s)", "Bypass jet velocity, $V_{jb}$ (m/s)",
         f"{fan.exit.V_jb:.1f}"),
        ("Fan diameter, d_f (m)", "Fan diameter, $d_f$ (m)", f"{fan.d_f:.3f}"),
        ("Mass flow rate, m_air (kg/s)", "Mass flow rate, $\\dot{m}_{air}$ (kg/s)",
         f"{fan.m_air:.1f}"),
        ("BPR", "BPR", f"{r['BPR']:g}"),
        ("Turbine inlet temperature, T_04 (K)", "Turbine inlet temperature, $T_{04}$ (K)",
         f"{p.core.T_04:.0f}"),
        ("Core jet velocity, V_jc (m/s)", "Core jet velocity, $V_{jc}$ (m/s)",
         f"{p.core.V_jc:.1f}"),
        ("Propulsive efficiency, eta_p", "Propulsive efficiency, $\\eta_p$",
         f"{two_stream_propulsive_efficiency(p):.4f}"),
        ("Overall efficiency, eta_o", "Overall efficiency, $\\eta_o$", f"{r['eta_ach']:.4f}"),
        ("Bypass nozzle area, A_b (m^2)", "Bypass nozzle area, $A_b$ (m$^2$)", f"{p.A_b:.3f}"),
        ("Core nozzle area, A_c (m^2)", "Core nozzle area, $A_c$ (m$^2$)", f"{p.core.A_c:.3f}"),
        ("ACF", "ACF", f"{fan.acf:.4f}"),
    ]


def print_design(r):
    for row in spec_rows(r):
        if row is None:
            print("  -- per engine --")
            continue
        label, _, value = row
        print(f"  {label:<38} {value}")
    extras = [
        ("F_req (kN) / thrust gap (kN)", f"{r['F_req'] / 1e3:.2f} / {r['gap'] / 1e3:+.2f}"),
        ("V_jc/V_jb", f"{r['ratio']:.3f}"),
        ("Mission eta_o", f"{r['eta_o']:.3f}"),
        ("Cargo (t) / payload (t)", f"{r['cargo'] / 1e3:.1f} / {payload_t(r['cargo']):.2f}"),
        ("TOW (t) / cruise fuel (t)", f"{r['TOW']:.2f} / {r['m_fb']:.2f}"),
    ]
    print("  -- design point --")
    for label, value in extras:
        print(f"  {label:<38} {value}")


def write_final_design_table(r):
    """The handout's Table 1 for the final design, as LaTeX."""
    lines = []
    for row in spec_rows(r):
        if row is None:
            lines.append("\\midrule\n & Per Engine \\\\\n\\midrule")
            continue
        _, tex_label, value = row
        lines.append(f"{tex_label:<46}& {value} \\\\")
    gas = ("a constant $c_p$ and $\\gamma = 1.4$" if GAS is CONSTANT_CP
           else "the ideal-gas enthalpy tables")
    latex = (
        "\\begin{table}[H]\n\\centering\n\\begin{tabular}{l c}\n\\toprule\n"
        "Specification & Your Design \\\\\n\\midrule\n" + "\n".join(lines)
        + "\n\\bottomrule\n\\end{tabular}\n"
        f"\\caption{{Specification sheet for the propulsion system. Final design with "
        f"{gas} from the HPC\nonwards, {r['cargo'] / 1e3:.0f} t cargo "
        f"(payload {payload_t(r['cargo']):.2f} t) and a take-off weight of "
        f"{r['TOW']:.1f} t.}}\n\\label{{tab:final_design}}\n\\end{{table}}\n")
    path = table_path("table1_final_design.tex")
    path.write_text(latex)
    print(f"\nLaTeX written to {path.relative_to(path.parent.parent)}")


def full_payload_table():
    print(f"STEP 1 -- eta_o sweep at the full payload ({payload_t(c.MAX_CARGO):.2f} t)")
    hdr = (f"{'eta_o':>6} {'TOW':>7} {'MTOW':>5} {'runs':>5} {'valid':>6} {'match':>6} | "
           f"{'FPR':>5} {'BPR':>5} {'Vjc/Vjb':>8} {'gap kN':>7} {'T_04':>6} {'status':>8}")
    print(hdr)
    print("-" * len(hdr))
    for eta_o in ETAS:
        rows = grid(eta_o)
        runs = [r for r in rows if "gap" in r]
        valid = [r for r in runs if r["valid"]]
        n_match = sum(r["matched"] for r in runs)
        head = (f"{eta_o:6.3f} {rows[0]['TOW']:7.2f} "
                f"{'ok' if rows[0]['mtow_ok'] else 'OVER':>5} {len(runs):5d} "
                f"{len(valid):6d} {n_match:6d} | ")
        if not runs:
            print(head + "cycle fails everywhere")
            continue
        if n_match:
            b, status = best_matched(runs), "MATCHED"
        elif valid:                     # valid, but V_jc/V_jb too far from 1
            b, status = min(valid, key=lambda r: abs(r["ratio"] - 1)), "Vjc/Vjb"
        else:                           # nearest thrust, and why it fails
            b = min(runs, key=lambda r: abs(r["gap"]))
            status = b["why"]
        print(head + f"{b['FPR']:5.2f} {b['BPR']:5.1f} {b['ratio']:8.3f} "
              f"{b['gap'] / 1e3:+7.2f} {b['T_04']:6.0f} {status:>8}")

def cargo_search():
    """Reduce cargo until a matched design exists; the best one at the first such cargo."""
    print(f"\nSTEP 2 -- cargo reduced until a matched design exists "
          f"(eta_o from {ETA_SEARCH[0]:.2f} in steps of 0.005)")
    hdr = (f"{'cargo t':>7} {'payload':>7} {'eta_o':>6} {'FPR':>5} {'BPR':>5} "
           f"{'Vjc/Vjb':>8} {'gap kN':>7} {'T_04':>6} {'TOW':>7}")
    print(hdr)
    print("-" * len(hdr))
    for cargo in CARGOS:
        b = highest_matched(cargo, ETA_SEARCH)
        if b is None:
            print(f"{cargo / 1e3:7.1f} {payload_t(cargo):7.2f}   none")
            continue
        print(f"{cargo / 1e3:7.1f} {payload_t(cargo):7.2f} {b['eta_o']:6.3f} {b['FPR']:5.2f} "
              f"{b['BPR']:5.1f} {b['ratio']:8.3f} {b['gap'] / 1e3:+7.2f} {b['T_04']:6.0f} "
              f"{b['TOW']:7.2f}")
        return b
    return None


def higher_mach_check(best):
    """Can a faster cruise, up to Q1's M_max, still find a matched design at this cargo?"""
    machs = np.round(np.arange(M + 0.01, ms.max_cruise_mach() + 1e-9, 0.01), 2)
    print(f"\nHIGHER MACH CHECK -- cargo held at {best['cargo'] / 1e3:.0f} t, FPR and BPR "
          f"re-searched, eta_o from {ETA_SEARCH[0]:.2f} in steps of 0.005")
    hdr = (f"{'M':>5} {'eta_o':>6} {'FPR':>5} {'BPR':>5} {'Vjc/Vjb':>8} {'gap kN':>7} "
           f"{'T_04':>6} {'TOW':>7}")
    print(hdr)
    print("-" * len(hdr))
    faster = []
    for mach in machs:
        b = highest_matched(best["cargo"], ETA_SEARCH, mach)
        if b is None:
            # eta_o cannot rise above the final design's without breaking the TET
            # ceiling, so the TOW at that eta_o shows whether MTOW is what blocks it
            tow, _ = ms.actual_tow(mach, best["eta_o"], best["cargo"])
            print(f"{mach:5.2f}   none -- TOW at eta_o {best['eta_o']:.3f} would be "
                  f"{tow / 1e3:.1f} t (MTOW {c.MTOW / 1e3:.0f} t)")
            continue
        faster.append(b)
        print(f"{mach:5.2f} {b['eta_o']:6.3f} {b['FPR']:5.2f} {b['BPR']:5.1f} "
              f"{b['ratio']:8.3f} {b['gap'] / 1e3:+7.2f} {b['T_04']:6.0f} {b['TOW']:7.2f}")
    if faster:
        print(f"Highest Mach with a matched design: {faster[-1]['M']:.2f}")
    else:
        print(f"No Mach above {M} gives a matched design at this cargo, "
              f"so M = {M} stands.")


def main():
    print(__doc__)
    print(f"Gas model: {GAS.name} (HPC to core nozzle); fan and bypass gamma = 1.4")
    print(f"M = {M}, FPR {FPR_SCAN[0]:.2f}-{FPR_SCAN[-1]:.2f}, "
          f"BPR {BPR_SCAN[0]:g}-{BPR_SCAN[-1]:g}, TET ceiling {TET_MAX:.0f} K, "
          f"tolerance +/-{THRUST_TOL / 1e3:.1f} kN, matched |Vjc/Vjb - 1| <= {RATIO_TOL}\n")

    full_payload_table()
    best = cargo_search()
    if best is None:
        print(f"\nNo matched design at M = {M} for any cargo from "
              f"{CARGOS[0] / 1e3:.0f} t down to 0.")
        return

    # refine eta_o between the coarse step that worked and the one above it
    eta_hi = min(best["eta_o"] + 0.005, ETAS[0])
    for eta_f in np.round(np.arange(eta_hi - FINE_STEP, best["eta_o"], -FINE_STEP), 4):
        b = best_matched(grid(eta_f, best["cargo"]))
        if b:
            best = b
            break

    print(f"\n{'=' * 78}\nFINAL DESIGN -- highest eta_o with a matched design "
          f"(refined to {FINE_STEP})\n{'=' * 78}")
    print_design(best)
    write_final_design_table(best)
    higher_mach_check(best)


if __name__ == "__main__":
    main()
