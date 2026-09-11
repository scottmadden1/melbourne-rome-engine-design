"""
Q5.1 Reduction of fuel burn: Melbourne to Rome with 0, 1 or 2 refuelling stops.

The final engine from Q5_final_design_sweep.py (its M, achieved eta_o and cargo)
flies the route split into n_stops + 1 legs. Following the handout:

  * every leg covers 180 nm from take-off to cruise and 150 nm from cruise to
    landing, and cruises for the rest of its length
  * take-off/climb and descent/landing fuel are 2% and 0.2% of the 0-stop
    take-off weight (not MTOW), on every leg
  * no headwind on any leg
  * cruise fuel from the exact Breguet range equation at the updated eta_o

Stops can be anywhere along the route. Breguet fuel grows exponentially with
cruise distance, so for a fixed number of stops equal legs burn the least.
"""

import math

from engine import constants as c
from engine import mission as ms
from engine.output import table_path
from Q5_gas_model_comparison import final_design

STOPS = (0, 1, 2)


def get_final_design():
    """M, achieved eta_o, cargo and the 0-stop take-off weight of the final design."""
    r = final_design()
    if r is None:
        raise RuntimeError("the final design sweep finds no design")
    return dict(M=r["M"], eta_o=r["eta_ach"], cargo=r["cargo"], TOW0=r["TOW"] * 1e3)


def leg_distances_nm(n_stops):
    """
    Great-circle length of each leg for n_stops refuelling stops, nm.

    """
    return ms.great_circle_nm() / (n_stops + 1)


def leg_cruise_distance_nm(leg_nm):
    """
    Cruise distance of one leg, nm: the leg less the 180 nm climb and 150 nm
    descent (c.DIST_TAKEOFF_NM, c.DIST_LANDING_NM).

    """
    return leg_nm - c.DIST_TAKEOFF_NM - c.DIST_LANDING_NM


def leg_cruise_fuel(M, eta_o, cargo, cruise_nm, m_landing):
    """
    Exact Breguet cruise fuel for one leg with no headwind, kg:

        m_fb = (m_e + m_pl + m_landing) * (exp(R/H) - 1),   H = eta_o*LCV*(L/D)/g

    with R = cruise_nm * c.NM_TO_M, m_pl = ms.payload(cargo), and m_landing the
    descent/landing fuel still aboard at the end of cruise.
    Compare ms.cruise_fuel_burn, which does the same for the whole route.

    """
    R = cruise_nm * c.NM_TO_M
    H = eta_o * c.LCV * ms.lift_to_drag(M) / c.G
    return (c.M_EMPTY + ms.payload(cargo) + m_landing) * (math.exp(R / H) - 1)


def journey_fuel(design, n_stops):
    """
    Fuel for the whole journey with n_stops stops, kg.

    Returns dict(n_stops, legs, takeoff, cruise, landing, total), where legs is a
    list of per-leg dicts with keys leg_nm, cruise_nm, takeoff, cruise, landing.

    """

    leg_nm = leg_distances_nm(n_stops)
    cruise_nm = leg_cruise_distance_nm(leg_nm)
    leg_cruise = leg_cruise_fuel(design["M"], design["eta_o"], design["cargo"], cruise_nm, design["TOW0"] * c.F_LANDING)
    leg_takeoff = design["TOW0"] * c.F_TAKEOFF
    leg_landing = design["TOW0"] * c.F_LANDING
    total_takeoff = leg_takeoff * (n_stops + 1)
    total_cruise = leg_cruise * (n_stops + 1)
    total_landing = leg_landing * (n_stops + 1)
    total_fuel = total_takeoff + total_cruise + total_landing
    leg = dict(leg_nm=leg_nm, cruise_nm=cruise_nm, takeoff=leg_takeoff,
               cruise=leg_cruise, landing=leg_landing)
    return dict(n_stops=n_stops, legs=[leg] * (n_stops + 1), takeoff=total_takeoff,
                cruise=total_cruise, landing=total_landing, total=total_fuel)



def print_results(design, results):
    print(f"Final design: M = {design['M']:.2f}, eta_o = {design['eta_o']:.4f}, "
          f"cargo = {design['cargo'] / 1e3:.0f} t, 0-stop TOW = {design['TOW0'] / 1e3:.2f} t")
    print(f"Route: {ms.great_circle_nm():.0f} nm great circle, no headwind\n")
    hdr = (f"{'stops':>5} {'leg nm':>8} {'take-off t':>11} {'cruise t':>9} "
           f"{'landing t':>10} {'total t':>8} {'saving':>8}")
    print(hdr)
    print("-" * len(hdr))
    base = results[0]["total"]
    for r in results:
        saving = 100 * (base - r["total"]) / base
        print(f"{r['n_stops']:5d} {r['legs'][0]['leg_nm']:8.0f} {r['takeoff'] / 1e3:11.2f} "
              f"{r['cruise'] / 1e3:9.2f} {r['landing'] / 1e3:10.2f} "
              f"{r['total'] / 1e3:8.2f} {saving:7.1f}%")


def write_table(design, results):
    """The fuel burn for each number of stops, as LaTeX."""
    base = results[0]["total"]
    rows = []
    for r in results:
        saving = 100 * (base - r["total"]) / base
        saving_s = "--" if r["n_stops"] == 0 else f"${saving:.1f}$"
        rows.append(f"{r['n_stops']} & {r['legs'][0]['leg_nm']:.0f} & "
                    f"{r['legs'][0]['cruise_nm']:.0f} & {r['takeoff'] / 1e3:.2f} & "
                    f"{r['cruise'] / 1e3:.2f} & {r['landing'] / 1e3:.2f} & "
                    f"{r['total'] / 1e3:.2f} & {saving_s} \\\\")
    latex = (
        "\\begin{table}[H]\n\\centering\n\\small\n\\setlength{\\tabcolsep}{4pt}\n"
        "\\begin{tabular}{c c c c c c c c}\n\\toprule\n"
        "Stops & Leg (nm) & Cruise (nm) & Take-off (t) & Cruise fuel (t) & Landing (t) "
        "& Total (t) & Saving (\\%) \\\\\n"
        "\\midrule\n" + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n"
        f"\\caption{{Total fuel burnt from Melbourne to Rome with 0, 1 and 2 refuelling "
        f"stops, using the final design\n($M = {design['M']:.2f}$, "
        f"$\\eta_o = {design['eta_o']:.4f}$, {design['cargo'] / 1e3:.0f} t cargo), "
        f"equal legs, no headwind, and take-off\nand landing fuel of 2\\% and 0.2\\% of "
        f"the non-stop take-off weight ({design['TOW0'] / 1e3:.1f} t) on every leg.}}\n"
        "\\label{tab:refuelling_stops}\n\\end{table}\n")
    path = table_path("table4_refuelling_stops.tex")
    path.write_text(latex)
    print(f"\nLaTeX written to {path.relative_to(path.parent.parent)}")


def main():
    print(__doc__)
    design = get_final_design()
    results = [journey_fuel(design, n) for n in STOPS]
    print_results(design, results)
    write_table(design, results)
    # TODO: compare the 0-stop, no-headwind total with the headwind cruise fuel
    # of the final design itself (ms.cruise_fuel_burn).


if __name__ == "__main__":
    main()
