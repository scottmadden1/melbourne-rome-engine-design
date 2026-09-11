"""
Q4.5 Sensitivity studies at a fixed FPR and BPR.

Each study varies ONE input and holds the others at their baseline, and
reports two fuel-air ratios that must be reconciled for a workable design:

  f_supplied  the fuel the mission budget allows,
                  f = m_f/m_core   with   m_f = F_req*V/(eta_o*LCV)
  f_required  the smallest f at which the engine makes F_req to within the
              handout's 4.5 kN tolerance

A design closes when the mission's own fuel makes the required thrust AND the
point is physically valid: the cycle closes, T_04 stays under the TET
ceiling, and TOW stays under MTOW.

  Study 1  Mach number        -- M varied, eta_o and cargo held
  Study 2  overall efficiency -- eta_o varied, M and cargo held
  Study 3  cargo weight       -- cargo varied, M and eta_o held

FPR and BPR are fixed at 1.50 and 8 in all three, so cargo and efficiency act
through the weight chain, not just the fuel budget. Everything downstream of the mission (TOW -> F_req -> d_f -> m_air -> m_core)
is recomputed in every case. The fan is sized with the matched ACF and the
core uses the enthalpy tables, including for the HPC exit.
"""

import math

import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt

from engine import ENTHALPY_TABLE, CycleError
from engine import constants as c
from engine import mission as ms
from engine.fan import matched_acf, size_fan
from engine.output import figure_path
from engine.performance import evaluate as engine_performance

# ------------------------------------------------------- fixed design point

FPR_FIXED, BPR_FIXED = 1.50, 8
TET_MAX = 1500.0       # K, turbine entry temperature ceiling (not in the handout)
GAS = HPC_GAS = ENTHALPY_TABLE

BASE_M, BASE_ETA, BASE_CARGO = 0.75, c.ETA_O_MISSION, c.MAX_CARGO
BASE_FPR, BASE_BPR = FPR_FIXED, BPR_FIXED

PAX_TOTAL = c.MAX_PASSENGERS * c.PAX_MASS
THRUST_TOL = 4.5e3     # N, the handout's design tolerance
F_SCAN = np.linspace(0.002, 0.045, 4000)


def mission(M, eta_o, cargo):
    """Take-off weight, cruise fuel and required thrust per engine for one mission."""
    TOW, m_fb = ms.actual_tow(M, eta_o, cargo)
    F_req = ms.required_thrust(M, TOW) if math.isfinite(TOW) else float("inf")
    return TOW, m_fb, F_req, ms.lift_to_drag(M), ms.flight_speed(M)


def size_engine(M, F_req, FPR, BPR):
    """Fan sized to F_req with the matched ACF for this FPR and BPR."""
    return size_fan(FPR, M, F_req, acf=matched_acf(FPR, BPR, M), bracket=(0.3, 9.0))


def net_thrust(fan, BPR, f):
    p = engine_performance(fan, BPR, f, GAS, HPC_GAS, TET_MAX)
    return p.F_net, p.core


def evaluate(M, eta_o, cargo, FPR=FPR_FIXED, BPR=BPR_FIXED):
    """Both fuel-air ratios and the validity verdict at one operating point."""
    TOW, m_fb, F_req, L_D, V = mission(M, eta_o, cargo)
    out = dict(M=M, eta_o=eta_o, cargo=cargo, TOW=TOW / 1e3, m_fb=m_fb,
               F_req=F_req, L_D=L_D, V=V, mtow_ok=TOW <= c.MTOW)
    if not np.isfinite(F_req):
        return {**out, "f_sup": float("nan"), "f_req": float("nan"),
                "note": "mission unflyable"}

    try:
        fan = size_engine(M, F_req, FPR, BPR)
    except ValueError:
        return {**out, "f_sup": float("nan"), "f_req": float("nan"),
                "note": "no fan diameter"}
    m_core = fan.m_air / (1 + BPR)

    m_f = ms.cruise_fuel_flow(M, F_req, eta_o)
    f_sup = m_f / m_core
    out.update(d_f=fan.d_f, m_air=fan.m_air, m_core=m_core, m_f=m_f, f_sup=f_sup)

    # Map the valid f window. It is bounded BELOW by the LPT (p_05 must clear
    # ambient) and ABOVE by the TET ceiling, and it can be very narrow -- so
    # F_N need never reach F_req at all. f_required is therefore the smallest
    # f meeting the 4.5 kN tolerance, not an exact root.
    valid = []
    for f in F_SCAN:
        try:
            F_N, st = net_thrust(fan, BPR, f)
        except CycleError:
            continue
        valid.append((f, F_N, st.T_04))

    if not valid:
        return {**out, "f_sup": f_sup, "f_req": float("nan"), "note": "no valid f"}

    f_lo, f_hi = valid[0][0], valid[-1][0]
    f_best, F_N_max, T04_best = max(valid, key=lambda t: t[1])
    out.update(f_min=f_lo, f_max=f_hi, F_N_max=F_N_max, f_at_max=f_best,
               T_04=T04_best, shortfall=F_req - F_N_max)

    # Smallest f meeting the tolerance, refined with a root solve so it is not
    # quantised to the scan grid.
    target = F_req - THRUST_TOL
    out["f_req"] = float("nan")
    for (fa, ga, _), (fb, gb, _) in zip(valid, valid[1:]):
        if (ga - target) * (gb - target) <= 0:
            out["f_req"] = brentq(lambda x: net_thrust(fan, BPR, x)[0] - target,
                                  fa, fb, xtol=1e-12)
            out["T_04"] = net_thrust(fan, BPR, out["f_req"])[1].T_04
            break
    else:
        within = [t for t in valid if t[1] >= target]
        if within:                      # already satisfied at the window floor
            out["f_req"] = within[0][0]
            out["T_04"] = within[0][2]

    # The decisive test: the thrust the engine makes on the mission's own fuel.
    out["sup_valid"] = f_lo <= f_sup <= f_hi
    if out["sup_valid"]:
        F_N_sup, st_s = net_thrust(fan, BPR, f_sup)
        out.update(F_N_sup=F_N_sup, gap_sup=F_N_sup - F_req, T_04_sup=st_s.T_04)
    else:
        out.update(F_N_sup=float("nan"), gap_sup=float("nan"), T_04_sup=float("nan"))
    return out


def table(title, var, values, build):
    print(f"\n{'='*104}\n{title}\n{'='*104}")
    hdr = (f"{var:>8} {'TOW':>7} {'MTOW':>5} {'F_req':>7} {'f_sup':>8} "
           f"{'f_req':>8} {'sup/req':>8} {'f_min':>8} {'f_max':>8} "
           f"{'F_N(sup)':>9} {'gap':>7} {'T_04':>6} {'verdict':>8}")
    print(hdr); print('-'*len(hdr))
    rows = []
    for v in values:
        r = evaluate(*build(v))
        rows.append((v, r))
        if "f_min" not in r:
            print(f"{v:8.3f} {r['TOW']:7.2f} {'--':>5} {'--':>7} {'--':>8} {'--':>8} "
                  f"{'--':>8} {'--':>8} {'--':>8} {'--':>9} {'--':>7} {'--':>6} "
                  f"{r.get('note','--'):>8}")
            continue
        fr = f"{r['f_req']:8.5f}" if np.isfinite(r["f_req"]) else f"{'--':>8}"
        rat = (f"{r['f_sup']/r['f_req']:8.3f}"
               if np.isfinite(r["f_req"]) else f"{'--':>8}")
        fn = (f"{r['F_N_sup']/1e3:9.2f}" if np.isfinite(r["F_N_sup"])
              else f"{'--':>9}")
        gp = (f"{r['gap_sup']/1e3:+7.2f}" if np.isfinite(r["gap_sup"])
              else f"{'--':>7}")
        t4 = (f"{r['T_04_sup']:6.0f}" if np.isfinite(r["T_04_sup"])
              else f"{'--':>6}")
        # closes only if the mission's own fuel makes the required thrust
        ok = (r["mtow_ok"] and r["sup_valid"]
              and np.isfinite(r["gap_sup"]) and abs(r["gap_sup"]) <= THRUST_TOL)
        why = "CLOSES" if ok else ("MTOW" if not r["mtow_ok"] else
              ("lean" if r["f_sup"] < r["f_min"] else
               "hot" if r["f_sup"] > r["f_max"] else "thrust"))
        print(f"{v:8.3f} {r['TOW']:7.2f} {'ok' if r['mtow_ok'] else 'OVER':>5} "
              f"{r['F_req']/1e3:7.2f} {r['f_sup']:8.5f} {fr} {rat} "
              f"{r['f_min']:8.5f} {r['f_max']:8.5f} {fn} {gp} {t4} {why:>8}")
    return rows


def main():
    print(__doc__)
    print(f"Fixed: FPR = {FPR_FIXED}, BPR = {BPR_FIXED}, variable ACF, "
          f"TET ceiling {TET_MAX:.0f} K")
    print(f"Baseline: M = {BASE_M}, eta_o = {BASE_ETA}, "
          f"payload = {(PAX_TOTAL + BASE_CARGO)/1e3:.2f} t, "
          f"FPR = {BASE_FPR}, BPR = {BASE_BPR}")

    machs = np.arange(0.60, 0.86, 0.01)
    etas = np.arange(0.24, 0.461, 0.02)
    # cargo swept 0 to the 21 t maximum, reported as the total payload in tonnes
    payloads = (PAX_TOTAL + np.linspace(0.0, c.MAX_CARGO, 11)) / 1e3

    r1 = table("STUDY 1 -- Mach number varied (all else at baseline)",
               "M", machs, lambda v: (v, BASE_ETA, BASE_CARGO, BASE_FPR, BASE_BPR))
    r2 = table("STUDY 2 -- overall efficiency varied (all else at baseline)",
               "eta_o", etas, lambda v: (BASE_M, v, BASE_CARGO, BASE_FPR, BASE_BPR))
    r3 = table("STUDY 3 -- cargo varied (all else at baseline)",
               "payload", payloads,
               lambda v: (BASE_M, BASE_ETA, max(0.0, v * 1e3 - PAX_TOTAL),
                          BASE_FPR, BASE_BPR))

    # ------------------------------------------------------------------ plots
    def series(rows):
        """
        x, f_supplied, f_required for the points where both are defined, plus
        the FULL requested x range. Points where the cycle has no valid f are
        dropped from the curves -- but the axis is still spanned over the whole
        sweep and the dead region shaded, so a range truncated by the TET
        ceiling reads as truncated rather than as never having been swept.
        """
        xs, sup, req, all_x, dead = [], [], [], [], []
        for v, r in rows:
            all_x.append(v)
            if "f_min" not in r or not np.isfinite(r.get("f_req", np.nan)):
                dead.append(v)
                continue
            xs.append(v); sup.append(r["f_sup"]); req.append(r["f_req"])
        return (np.array(xs), np.array(sup), np.array(req),
                np.array(all_x), np.array(dead))

    def span_full(ax, all_x, dead):
        """Span the axis over the whole sweep and grey out the unsolvable part."""
        if len(all_x) < 2:
            return
        pad = 0.02 * (all_x.max() - all_x.min())
        ax.set_xlim(all_x.min() - pad, all_x.max() + pad)
        if len(dead) == 0:
            return
        step = np.min(np.diff(np.sort(all_x))) if len(all_x) > 1 else 0.0
        d = np.sort(dead)
        runs, start, prev = [], d[0], d[0]
        for x in d[1:]:
            if x - prev > step * 1.5:
                runs.append((start, prev)); start = x
            prev = x
        runs.append((start, prev))
        for k, (lo, hi) in enumerate(runs):
            ax.axvspan(lo - step/2, hi + step/2, color="0.6", alpha=0.28,
                       hatch="//", lw=0, zorder=0,
                       label="no valid $f$" if k == 0 else None)

    panels = ((r1, "Cruise Mach number, $M$", "Mach"),
              (r2, "Overall efficiency, $\\eta_o$", "Overall Efficiency"),
              (r3, "Cargo + passenger mass (t)", "Cargo/Passenger weight"))

    # ---- figure 1: the two fuel-air ratios
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    for ax, (rows, xlabel, name) in zip(axes, panels):
        xs, sup, req, all_x, dead = series(rows)
        if len(xs) < 2:
            ax.text(0.5, 0.5, "insufficient valid points", ha="center",
                    va="center", transform=ax.transAxes)
            ax.set_title(name); continue
        span_full(ax, all_x, dead)
        ax.plot(xs, sup, marker="o", ms=4, color="tab:orange", lw=1.8,
                label="$f_{supplied}$")
        ax.plot(xs, req, marker="s", ms=4, color="tab:red", lw=1.8,
                label="$f_{required}$")
        ax.fill_between(xs, req, sup, where=sup >= req, alpha=0.18,
                        color="tab:green", interpolate=True, label="surplus")
        ax.fill_between(xs, req, sup, where=sup < req, alpha=0.18,
                        color="tab:red", interpolate=True, label="deficit")
        ax.set_xlabel(xlabel); ax.set_ylabel("fuel-air ratio, $f$")
        ax.set_title(name)
        ax.grid(True, alpha=0.3); ax.legend(fontsize=9, loc="best")
    fig.tight_layout()
    fig.savefig(figure_path("sensitivity_studies.png"), dpi=150, bbox_inches="tight")

    # ---- figure 2: the margin between them
    fig2, ax2 = plt.subplots(1, 3, figsize=(16, 4.8))
    print(f"\n{'='*72}\nFuel margin  f_sup - f_req\n{'='*72}")
    print(f"{'study':>22} {'margin range':>26} {'swept':>16} {'solvable':>16} "
          f"{'dead':>7}")
    print('-' * 92)

    for ax, (rows, xlabel, name) in zip(ax2, panels):
        xs, sup, req, all_x, dead = series(rows)
        if len(xs) < 2:
            ax.text(0.5, 0.5, "insufficient valid points", ha="center",
                    va="center", transform=ax.transAxes)
            ax.set_title(name)
            print(f"{name:>22}   insufficient valid points")
            continue
        diff = sup - req
        span_full(ax, all_x, dead)
        ax.axhline(0, color="k", lw=1.2, ls="--", zorder=1)
        ax.plot(xs, diff, marker="o", ms=4, color="tab:purple", lw=1.9,
                label="$f_{sup} - f_{req}$")
        ax.fill_between(xs, 0, diff, where=diff >= 0, alpha=0.20,
                        color="tab:green", interpolate=True, label="surplus")
        ax.fill_between(xs, 0, diff, where=diff < 0, alpha=0.20,
                        color="tab:red", interpolate=True, label="deficit")

        ax.set_xlabel(xlabel); ax.set_ylabel("$f_{sup} - f_{req}$")
        ax.set_title(name)
        ax.grid(True, alpha=0.3); ax.legend(fontsize=9, loc="best")

        rng = f"{all_x.min():g} .. {all_x.max():g}"
        sol = f"{xs.min():g} .. {xs.max():g}"
        print(f"{name:>22} {f'{diff.min():+.5f} .. {diff.max():+.5f}':>26} "
              f"{rng:>16} {sol:>16} {len(dead):>7}")

    fig2.tight_layout()
    fig2.savefig(figure_path("sensitivity_difference.png"), dpi=150, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
