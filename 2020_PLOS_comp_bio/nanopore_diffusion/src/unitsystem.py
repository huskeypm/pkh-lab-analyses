"""Units and physical constants

**TODO**: list defined units and constants"""

import scipy.constants

pc=scipy.constants.physical_constants

#Boltzmann constant
kB=pc['Boltzmann constant in eV/K'][0]

#Avogadro's number
avogadro=pc['Avogadro constant'][0]

#Unit conversion factors
#These are designed to go from the designated unit to model units
#This way, inputs are specfied as X*unit to convert them from "unit" to the corresponding model uit
#and outputs are converted to other units as Y/unit.

#Length
meter = 1e9
micron = 1e3
nanometer = 1

#Number of particles
particle = 1
mole = avogadro
dozen = 12

#Time
nanosecond=1
microsecond=1e3
millisecond=1e6
second=1e9

#Energy
electronVolt = 1
joule = 1.0/pc['electron volt'][0]

#Temperature
kelvin = 1

#Electric charge
coulomb = 1.0/pc['elementary charge'][0]

#Concentration
molar = avogadro/1e24
millimolar = molar*1e-3
micromolar = molar*1e-6

#Electric potential
volt = 1
millivolt = 1e-3

#Vacuum permittivity
eps_0_SI = pc['electric constant'][0]
eps_0 = eps_0_SI * coulomb / volt / meter
