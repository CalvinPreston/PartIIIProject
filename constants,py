'''Script to define all constants to keep things clean'''
import numpy as np

h=0.6744
H0= h / ( 9.78*10**(9) * 365 * 24 * 60**2 )
G = 6.674*10**(-11)
c = 3*10**(8)
RHOB = 4.185*10**(-28)      # Baryon density in LCDM in kg/m^3
RHOC = 2.248*10**(-27)      # CDM density in LCDM in kg/m^3
RHOL = 5.870*10**(-27)      # Cosmological constant energy density in LCDM in kg/m^3
RHOR = 4.64511*10**(-31)    # Radiation density in LCDM in kg/m^3
OmegaM = 0.31               # Density parameter of ALL matter in LCDM from PLANCK
OmegaL = 0.69               # Density parameter of ALL dark energy in LCDM from PLANCK
omegab = 0.02226            # Density parameter for baryons multipled by h**2 from PLANCK
omegac = 0.1196             # Density parameter for LCDM multipled by h**2 from PLANCK
Mpc = 3.0857*10**(19)       # Parsec definition
zstart = 10                 # Redshift from which to consider start

PreFac = ((8*np.pi*G)/3)

TStart = (2 / (3 * H0 * (1-OmegaM)) ) * np.arcsinh(  ((1-OmegaM) / OmegaM)**0.5 * (1 + zstart)**(-1.5) )
TEnd = (15*10**9)*(365*24*60**2)
