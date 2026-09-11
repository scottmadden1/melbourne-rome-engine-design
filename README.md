# Final code: map to the report

The physics lives once, in the `engine/` package. Each `Q*.py` script is a
thin driver that runs a sweep, prints results, and saves figures.

```
final code/
├── engine/                     shared library (no printing or plotting)
│   ├── constants.py            handout values, all in SI units
│   ├── mission.py              L/D, cruise distance, Breguet fuel burn, required thrust
│   ├── gas.py                  ConstantCp and EnthalpyTable gas models
│   ├── fan.py                  intake, fan exit, bypass nozzle, eq. (1) sizing, matched ACF
│   ├── core.py                 booster, HPC, combustor, HPT, LPT, core nozzle
│   ├── thrust.py               eq. (1) and (2) formulas, installation drag
│   ├── performance.py          sized fan plus core, combined through eq. (2)
│   └── output.py               figure and table output paths
├── Q1_flight_mach_number.py    Q4.1
├── Q2_fuel_burn_metrics.py     Q4.2
├── Q3_fan_design.py            Q4.3
├── Q4_core_design.py           Q4.4, both gas models, writes Table 2
├── Q5_sensitivity_studies.py   Q4.5
├── Q5_final_design_sweep.py    Q4.5 final design: eta_o sweep at M = 0.75
├── Q5_gas_model_comparison.py  Q4.5 final design, enthalpy tables vs constant cp
├── Q6_reduction_of_fuel_burn.py Q5.1 fuel burn with 0, 1 or 2 refuelling stops, writes Table 4
├── tests/test_regression.py    checks the library against the pre-refactor numbers
├── figures/                    written by the scripts
└── tables/                     written by the scripts
```

## Running

From this folder, with the assignment's virtual environment:

```bash
../.venv/bin/python Q4_core_design.py
../.venv/bin/python tests/test_regression.py
```

Figures and tables are always written to `figures/` and `tables/` here,
whatever folder you run from.

## How the questions chain together

| Script | Uses | Produces |
|---|---|---|
| Q1 | mission model | cruise Mach number (0.86) |
| Q2 | Q1's Mach, via `mission.max_cruise_mach()` | thrust per engine and fuel flow |
| Q3 | Q2's thrust | fan diameter and mass flow over FPR, ACF = 1 |
| Q4 | Q2's thrust and fuel flow, Q3's fan sizing | core stations for both gas models, Table 2 |
| Q5 | mission and cycle, with its own baseline (M = 0.75, eta_o = 0.46, full payload) | Mach, efficiency and cargo sensitivity tables and figures |
| Q5 final design sweep | mission and cycle at M = 0.75, the Q4 FPR/BPR design over an eta_o sweep | highest eta_o with a valid V_jc ~ V_jb design, the cargo it costs, and the final design as Table 1 |
| Q5 gas model comparison | the final design sweep's final design | that design evaluated with both gas models, station by station, with every constraint checked (Table 3) |

Nothing downstream is typed in by hand, so a change to a handout value in
`engine/constants.py` flows through every question.

## The gas models

The core cycle in `engine/core.py` is written once in enthalpy form. Passing
`CONSTANT_CP` makes h = cp*T with gamma = 1.4, which reduces every balance to
the familiar constant-cp temperature form. Passing `ENTHALPY_TABLE` uses the
Cengel A-17 air table. The fan and booster always compress with gamma = 1.4.
The HPC does too, unless `hpc_gas=ENTHALPY_TABLE` is passed, which
`Q5_sensitivity_studies.py` does and Q4 does not. That is why T_03 is
identical in both columns of Table 2. `Q5_final_design_sweep.py` uses one
model from the HPC onwards, set by `GAS` (constant cp by default), and
`Q5_gas_model_comparison.py` runs it with both.

## Differences from the original scripts

- **Chained inputs.** Q3 and Q4 now use the exact 58.364 kN and 0.748797 kg/s
  from the mission model instead of the rounded 58.36 kN and 0.7488 kg/s typed
  in by hand. Table 2 shifts in the second decimal place as a result, for
  example net thrust with enthalpy tables goes from 44.08 to 44.07 kN.
  `tests/test_regression.py` confirms that the old inputs still give the old
  numbers exactly.
- **Q1 label.** The 79.80 t line is the fuel available for cruise after
  take-off and landing allowances, so it is no longer labelled as fuel burnt.
- **Q3 output** is a table per FPR instead of two lines of choke messages per FPR.
- **Q4** is one script running both gas models, and it writes
  `tables/table2_cp_comparison.tex` directly.
- **Q5 studies** are Mach number, overall efficiency and cargo, at FPR = 1.50,
  BPR = 8 and a baseline eta_o of 0.46. The original script's docstring
  described these, but its code ran Mach, FPR and BPR studies instead.
  f_required is the smallest f within the 4.5 kN tolerance, so it sits about
  0.0004 below the report's Figures 4 and 5, which came from an earlier
  version of the script.
- **Not carried over:** the Mach sweep in the old `f_vs_mach.py` (only its
  cycle functions were used by the report) and its unused "entropy" HPC
  option. The originals remain in the parent `Assignment 1` folder.
