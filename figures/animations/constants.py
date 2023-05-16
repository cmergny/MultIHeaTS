"""
File that stores all constants of the code.
Do not modify unless absolutely sure of what you are doing.
"""
# ASTRONOMIC
UA = 149.597e9  # m
SOLAR_CST = 1.361e3  # W/m^2
YEAR = 86400 * 365  # s (Earth)

# EARTH
EARTH_DAY = 86400  # s
EARTH_ALB = 0.015
EARTH_DIST = UA  # m

# PHYSICS
SIGMA = 5.670374419e-8  # Stephan Boltzmann W.m-2.K-4
RGAS = 8.31446261815324  # Gas constant J.K-1.mol-1

# WATER
MOL_WATER = 18e-3  # Molar mass of water kg.mol-1

# METAMORPHISM
P_WATER_TRIPLE = 611.73  # (Pa) for triple pt of water
T_WATER_TRIPLE = 273.16  # (K) for triple pt of water
D_VAPOR_AIR = 21.0e-6  # (m2/s) diff of water vapor in air Collbeck 93
D_VAPOR_SNOW = 85.0e-6  # (m2/s) diff of water vapor in snow (ok for ice?)
K_ICE = 2.2  # Conductivity of ice (W.m-1.K-1)
K_AIR = K_ICE / 100  # Conductivity of air (W.m-1.K-1)
RHO_ICE = 917  # density pure ice SI (valid for T=0 C only!!)
LSUB_ICE = 2.838e6  # (J/kg) latent heat sublimation of ice
R_ICE = 461.9  # (?) gas constant for ice
GAMMA_S = 0.06  # (J/m2) surfance tension /!\ may vary in litt
## SNOWPACK PARAMS
FGG = 3  # geom factor eq9 snowpackII
FGB = 0.35  # geom factor eq14 snowpackII
NSGZ = 0.15e-3  # (m) new snow grain size
PSAT_C2 = 21.88  # empirical coeff found in Murrar1966
PSAT_C3 = 7.66  # empirical coeff found in Murrar1966
