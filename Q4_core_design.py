"""
Q4.4 Core design, with constant cp and with ideal-gas enthalpy tables.

The fan is sized once per FPR from eq. (1) with ACF = 1, which involves no cp.
The same core cycle is then run with each gas model, so the only difference
between the two sets of results is the gas property model. The design point
comparison is printed and written as the LaTeX source of Table 2.
"""

import math

import matplotlib.pyplot as plt
import numpy as np

from engine import CONSTANT_CP, ENTHALPY_TABLE, CycleError
from engine import constants as c
from engine import mission as ms
from engine.fan import size_fan
from engine.output import figure_path, table_path
from engine.performance import evaluate, fuel_air_ratio

GAS_MODELS = (CONSTANT_CP, ENTHALPY_TABLE)
FPR_GRID = np.round(np.arange(1.45, 1.605, 0.01), 2)
BPR_GRID = range(1, 9)
DESIGN_FPR, DESIGN_BPR = 1.60, 8

# velocity ratio figure, constant cp only, as used in the report
PLOT_BPR = np.arange(1.0, 14.01, 0.5)
PLOT_SPEC = {
    CONSTANT_CP: dict(fprs=(1.45, 1.50, 1.55, 1.60), reference_lines=False,
                      file="vjc_vjb_vs_bpr_single_cp.png"),
}


def run(fan, BPR, m_f, gas):
    """Engine performance with the mission fuel flow burnt in the core."""
    return evaluate(fan, BPR, fuel_air_ratio(fan, BPR, m_f), gas)


def grid_search(fans, m_f, gas):
    """Closest thrust match over the FPR/BPR grid, overall and with V_jc > V."""
    best = best_valid = None
    for FPR in FPR_GRID:
        for BPR in BPR_GRID:
            try:
                p = run(fans[FPR], BPR, m_f, gas)
            except CycleError:
                continue
            gap = abs(p.F_net - fans[FPR].F_req)
            if best is None or gap < best[0]:
                best = (gap, p)
            if p.core.V_jc > p.fan.V and (best_valid is None or gap < best_valid[0]):
                best_valid = (gap, p)
    return best, best_valid


def velocity_ratio_curve(fan, m_f, gas):
    ratios = []
    for BPR in PLOT_BPR:
        try:
            ratios.append(run(fan, BPR, m_f, gas).velocity_ratio)
        except CycleError:
            ratios.append(np.nan)
    return np.array(ratios)


def plot_velocity_ratio(fans, m_f, gas):
    spec = PLOT_SPEC[gas]
    fig, ax = plt.subplots(figsize=(8, 6))
    for k, FPR in enumerate(spec["fprs"]):
        fan = fans[FPR]
        ratio = velocity_ratio_curve(fan, m_f, gas)
        ax.plot(PLOT_BPR, ratio, marker="o", markersize=4, color=f"C{k}",
                label=f"FPR = {FPR:.2f}")

        ok = ~np.isnan(ratio)
        cross = np.where(ok[:-1] & ok[1:] & (ratio[:-1] < 1) & (ratio[1:] >= 1))[0]
        if cross.size:
            i = cross[0]
            bpr_match = np.interp(1.0, ratio[i:i + 2], PLOT_BPR[i:i + 2])
            print(f"  FPR = {FPR:.2f}: V_jc = V_jb at BPR = {bpr_match:.2f}")
            if spec["reference_lines"]:
                ax.axvline(bpr_match, color="gray", linestyle=":", linewidth=1)
                ax.annotate(f"$V_{{jc}} = V_{{jb}}$ at BPR $\\approx$ {bpr_match:.1f}",
                            xy=(bpr_match, 1.0), xytext=(8, -25), textcoords="offset points")
        else:
            print(f"  FPR = {FPR:.2f}: V_jc never reaches V_jb for BPR "
                  f"{PLOT_BPR[0]:.1f}-{PLOT_BPR[-1]:.1f}")
        if np.isnan(ratio).any():
            print(f"  FPR = {FPR:.2f}: no valid core for BPR <= "
                  f"{PLOT_BPR[np.isnan(ratio)].max():.1f} (p_05 < p_a)")

        if spec["reference_lines"]:
            ax.axhline(1.0, color="gray", linestyle="--", linewidth=1,
                       label="$V_{jc} = V_{jb}$")
            ax.axhline(fan.V / fan.exit.V_jb, color="tab:red", linestyle=":",
                       linewidth=1, label="$V_{jc} = V$ (core thrust = 0)")

    ax.set_xlabel("Bypass ratio, BPR")
    ax.set_ylabel("$V_{jc}/V_{jb}$")
    ax.set_title(f"Core-to-bypass jet velocity ratio vs BPR ({gas.name})")
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(figure_path(spec["file"]), dpi=150)


def comparison_rows(a, b):
    """Table 2 rows: (label, LaTeX label, unit, LaTeX unit, value A, value B, format)."""
    def row(label, tex, unit, tex_unit, get, fmt, scale=1.0):
        return label, tex, unit, tex_unit, get(a) * scale, get(b) * scale, fmt

    def nozzle(p):
        return "Choked" if p.core.choked else "Unchoked"

    return [
        row("T_03", "$T_{03}$", "K", "K", lambda p: p.core.T_03, ".2f"),
        row("T_04", "$T_{04}$", "K", "K", lambda p: p.core.T_04, ".2f"),
        row("T_045", "$T_{045}$", "K", "K", lambda p: p.core.T_045, ".2f"),
        row("p_045", "$p_{045}$", "kPa", "kPa", lambda p: p.core.P_045, ".2f", 1e-3),
        row("T_05", "$T_{05}$", "K", "K", lambda p: p.core.T_05, ".2f"),
        row("p_05", "$p_{05}$", "kPa", "kPa", lambda p: p.core.P_05, ".2f", 1e-3),
        ("Core nozzle", "Core nozzle", "-", "--", nozzle(a), nozzle(b), None),
        row("V_jc", "$V_{jc}$", "m/s", "m/s", lambda p: p.core.V_jc, ".2f"),
        row("V_jc/V_jb", "$V_{jc}/V_{jb}$", "-", "--", lambda p: p.velocity_ratio, ".3f"),
        row("A_c", "$A_c$", "m^2", "m$^2$", lambda p: p.core.A_c, ".4f"),
        row("F_N", "$F_N$", "kN", "kN", lambda p: p.F_net, ".2f", 1e-3),
    ]


def signed(d, fmt):
    """Difference with an explicit sign, and a plain zero if it rounds to zero."""
    return f"{0:{fmt}}" if float(f"{abs(d):{fmt}}") == 0 else f"{d:+{fmt}}"


def print_and_write_table2(a, b):
    rows = comparison_rows(a, b)
    print(f"\nTable 2: design point FPR = {DESIGN_FPR}, BPR = {DESIGN_BPR} "
          f"(F_N is eq. (2) less nacelle and weight drag)")
    print(f"{'Quantity':<12} {'Unit':<5} {CONSTANT_CP.name:>15} "
          f"{ENTHALPY_TABLE.name:>16} {'Difference':>11}")
    tex = []
    for label, tex_label, unit, tex_unit, va, vb, fmt in rows:
        if fmt is None:
            va_s, vb_s, d_s = va, vb, "-"
            d_tex = "--"
        else:
            va_s, vb_s, d_s = f"{va:{fmt}}", f"{vb:{fmt}}", signed(vb - va, fmt)
            d_tex = f"${d_s}$"
        print(f"{label:<12} {unit:<5} {va_s:>15} {vb_s:>16} {d_s:>11}")
        tex.append(f"{tex_label:<22}& {tex_unit:<8}& {va_s:<8}& {vb_s:<9}& {d_tex:<10}\\\\")

    latex = (
        "\\begin{table}[H]\n\\centering\n\\begin{tabular}{l c c c c}\n\\toprule\n"
        "Quantity & Unit & Constant $c_p$, $\\gamma = 1.4$ & Enthalpy tables & Difference \\\\\n"
        "\\midrule\n" + "\n".join(tex) + "\n\\bottomrule\n\\end{tabular}\n"
        f"\\caption{{Engine conditions at FPR $= {DESIGN_FPR:g}$, BPR $= {DESIGN_BPR}$ and "
        f"$M = {a.fan.M:.2f}$ using a constant\n$c_p$ and $\\gamma$ compared with the "
        "ideal-gas enthalpy tables.}\n\\label{tab:cp_comparison}\n\\end{table}\n")
    path = table_path("table2_cp_comparison.tex")
    path.write_text(latex)
    print(f"\nLaTeX written to {path.relative_to(path.parent.parent)}")


def main():
    M = ms.max_cruise_mach()
    F_req = ms.required_thrust(M)
    m_f = ms.cruise_fuel_flow(M, F_req)
    fans = {FPR: size_fan(FPR, M, F_req, acf=1.0) for FPR in FPR_GRID}
    fan = fans[DESIGN_FPR]
    print(f"M = {M:.2f}, V = {fan.V:.2f} m/s, F_req = {F_req/1e3:.2f} kN, "
          f"m_f = {m_f:.4f} kg/s per engine, ACF = 1")
    print(f"Design fan at FPR = {DESIGN_FPR}: d_f = {fan.d_f:.4f} m, "
          f"m_air = {fan.m_air:.2f} kg/s, V_jb = {fan.exit.V_jb:.2f} m/s, "
          f"fan face area = {fan.area:.4f} m^2")

    design = {}
    for gas in GAS_MODELS:
        print(f"\n{'=' * 72}\n{gas.name}\n{'=' * 72}")
        p = run(fan, DESIGN_BPR, m_f, gas)
        design[gas] = p
        print(f"At FPR = {DESIGN_FPR}, BPR = {DESIGN_BPR}:")
        print(f"  V_jc = {p.core.V_jc:.2f} m/s, V_jc/V_jb = {p.velocity_ratio:.4f}, "
              f"core nozzle {'choked' if p.core.choked else 'unchoked'}")
        print(f"  A_c = {p.core.A_c:.4f} m^2, A_b = {p.A_b:.4f} m^2")
        print(f"  eq. (2) bare thrust = {p.F_bare/1e3:.2f} kN, installation drag = "
              f"{fan.D_install/1e3:.2f} kN")
        print(f"  net thrust = {p.F_net/1e3:.2f} kN vs required {F_req/1e3:.2f} kN "
              f"(difference {(p.F_net - F_req)/1e3:+.2f} kN)")

        if gas in PLOT_SPEC:
            print("Velocity ratio crossings:")
            plot_velocity_ratio(fans, m_f, gas)

        best, best_valid = grid_search(fans, m_f, gas)
        b = best[1]
        print(f"Minimum thrust difference over the grid: {best[0]/1e3:.2f} kN "
              f"at FPR = {b.fan.FPR:.2f}, BPR = {b.BPR}")
        if best_valid is None:
            print("No grid point gives a physically valid core (V_jc > V).")
        else:
            v = best_valid[1]
            eta_o = v.fan.V * v.F_net / (m_f * c.LCV)
            print(f"Best match with V_jc > V: {best_valid[0]/1e3:.2f} kN at "
                  f"FPR = {v.fan.FPR:.2f}, BPR = {v.BPR}, F_N = {v.F_net/1e3:.2f} kN, "
                  f"eta_o = {eta_o:.4f}")

    print_and_write_table2(design[CONSTANT_CP], design[ENTHALPY_TABLE])
    plt.show()


if __name__ == "__main__":
    main()
