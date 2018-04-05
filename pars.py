import numpy as np

# gN = 6.67408e-08 #Newtons gravitational constant
# mSun= 1.9884754153381438e+33 #mass of the Sun (in grams)
# mEarth = 5.9722e5 #mass of Earth (in grams)
# au = 1.495978707e13 #astronomical unit (in cm)
#
# yr = 2*np.pi /np.sqrt(gN*mSun/au**3) #1 year in seconds

mSun = 1.                # Mass of Sun (in Solar mass)
gN   = 4. * np.pi**2      # Newtons gravitational constant [au]3[yr]-2[Solarmass]-1
au   = 1.                # Astronomical unit (in AU)
yr   = 1.                # Year (in yr)
mEarth = 1/1.9884754153381438e+33   # mass of Earth (in Solar mass)




#problem parameters
Np = 2 #2 particles