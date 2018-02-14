import numpy as np

gN = 6.67408e-08 #Newtons gravitational constant
mSun= 1.9884754153381438e+33 #mass of the Sun (in grams)
mEarth = 5.9722e27 #mass of Earth (in grams)
au = 1.495978707e13 #astronomical unit (in cm)
yr = 2*np.pi /np.sqrt(gN*mSun/au**3) #1 year in seconds

#problem parameters
Np = 2 #2 particles