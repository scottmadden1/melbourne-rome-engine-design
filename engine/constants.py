import math

# --- ideal gas (handout)
R = 287.05              # J/(kg K)
GAMMA = 1.4
CP = 1004.5             # J/(kg K)

# --- static air at 37,000 ft (handout)
T_A = 216.65            # K
P_A = 21.66e3           # Pa
A_A = math.sqrt(GAMMA * R * T_A)    # speed of sound, m/s
RHO_A = P_A / (R * T_A)             # kg/m^3, from the gas law (0.3483)

G = 9.81                # m/s^2
LCV = 43e6              # J/kg

# --- component efficiencies and pressure ratios (handout)
ETA_INTAKE = 1.0
ETA_FAN = 0.9
ETA_COMPRESSOR = 0.9
ETA_TURBINE = 0.9
PI_FAN_LPC = 2.3        # fan and booster combined
PI_HPC = 26.0

# --- installation penalties (handout)
NACELLE_DRAG_PER_KGS = 9.25         # N per kg/s of air
ENGINE_WEIGHT_REF = 12000.0         # kg, engine mass at d_f = 3 m
ENGINE_WEIGHT_EXP = 2.4

# --- aircraft and mission (handout)
N_ENGINES = 2
MTOW = 243e3            # kg, maximum take-off weight
M_EMPTY = 115e3         # kg
MAX_PASSENGERS = 230
PAX_MASS = 95.0         # kg per passenger
MAX_CARGO = 21e3        # kg
F_TAKEOFF = 0.02        # fraction of MTOW burnt in take-off and climb
F_LANDING = 0.002       # fraction of MTOW burnt in descent and landing
ETA_O_MISSION = 0.46    # minimum overall efficiency
HEADWIND_KN = 30.0

MELBOURNE = (-37.67, 144.85)        # (latitude, longitude), degrees
ROME = (41.90, 12.48)
DIST_TAKEOFF_NM = 180.0
DIST_LANDING_NM = 150.0

# --- conversions
NM_TO_M = 1852.0
KM_TO_NM = 0.539957
KNOT_TO_MPS = 0.514444
R_EARTH_KM = 6371.0
