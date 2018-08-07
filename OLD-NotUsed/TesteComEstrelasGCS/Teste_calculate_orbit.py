# Import libraries
import numpy as np
import galpy
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014 as MW14
from matplotlib import pyplot as plt
from druvw import uvw
# Galactic model parameters
ro = 8.0  # kpc
vo = 220.0  # km/s
ts = np.linspace(0, 100, 200000) # Integration time. Translates to approx 7 Gyr
pot = MW14
solarmotion = 'schoenrich' # or 'hogg' or 'dehnen', or value in [-U,V,W]

# Dados observacionais
# GCS
ra = 1.280833 # deg
dec = -52.151667 # deg
plx =  23.86/1000. # arcsec
d = (1./plx) / 1000. # kpc
pmra = -111 # mas/yr
pmde = -131 # mas/yr
rv = 34.3 # km/s

eGCS = 0.16
RminGCS = 6.39 # kpc
RmaxGCS = 8.88 # kpc
zmaxGCS = 0.14
radcon  
UGCS = 40 # km/s
VGCS = -22 # km/s
WGCS = -16 # km/s

U, V, W = uvw(ra, dec, d*1000, pmra, pmde, rv)

vxvv = [ra, dec, d, pmra, pmde, rv]
orbit = Orbit(vxvv=vxvv, radec = True, ro = ro, vo = vo, solarmotion = solarmotion)
orbit.integrate(ts, pot, method = 'leapfrog')

vxvv2 = [ra, dec, d, UGCS, VGCS, WGCS]
orbit2 = Orbit(vxvv=vxvv2, radec = True, uvw = True, ro = ro, vo = vo, solarmotion = solarmotion)
orbit2.integrate(ts, pot, method = 'leapfrog')

vxvv3 = [ra, dec, d, U, V, W]
orbit3 = Orbit(vxvv=vxvv3, radec = True, uvw = True, ro = ro, vo = vo, solarmotion = solarmotion)
orbit3.integrate(ts, pot, method = 'leapfrog')

print "eGCS: {0}".format(eGCS)
print "eobs: {0}".format(orbit.e())
print "eobs(UVWcalc): {0}".format(orbit2.e())
print "eobs(UVWGCS): {0}".format(orbit3.e())

x = orbit.x(ts) # kpc
y = orbit.y(ts) # kpc
z = orbit.z(ts) # kpc
r = np.sqrt(x*x+y*y) # kpc
    
