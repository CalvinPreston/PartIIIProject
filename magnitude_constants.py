# Script to define all constants to keep things clean (for magnitude plot)
import numpy as np

h=0.6742
H0= h / ( 9.78*10**(9) * 365 * 24 * 60**2 )
G = 6.674*10**(-11)
c = 3*10**(8)
RHOB = 4.185*10**(-28)      # Baryon density in LCDM in kg/m^3
RHOC = 2.248*10**(-27)      # CDM density in LCDM in kg/m^3
RHOL = 5.870*10**(-27)      # Cosmological constant energy density in LCDM in kg/m^3
RHOR = 4.64511*10**(-31)    # Radiation density in LCDM in kg/m^3
OmegaM = 0.31584412154822133               # Density parameter of ALL matter in LCDM from PLANCK
OmegaL = 1-OmegaM              # Density parameter of ALL dark energy in LCDM from PLANCK
OmegaC = 0.26502750760903215
OmegaB = 0.049388982622703914

omegab = 0.022383          # Density parameter for baryons multipled by h**2 from PLANCK
omegac = 0.12011             # Density parameter for LCDM multipled by h**2 from PLANCK
Mpc = 3.0857*10**(22)       # Megaparsec definition
Kpc = 3.0857*10**(19)       # Kilaparsec definition
Pc = 3.0857*10**(16)  # Parsec definition

zstart = 10               # Redshift from which to consider start
MB = -19.387
PreFac = ((8*np.pi*G)/3)

TStart = (2 / (3 * H0 * (1-OmegaM)) ) * np.arcsinh(  ((1-OmegaM) / OmegaM)**0.5 * (1 + zstart)**(-1.5) )
TEnd = (13.9*10**9)*(365*24*60**2)
