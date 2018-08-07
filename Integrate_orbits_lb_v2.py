# Import libraries
import numpy as np
import galpy
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014 as MW14
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib import pyplot as plt

# Galactic model parameters
ro = 8.0  # kpc
vo = 220.0  # km/s
ti = -0.15 # Gyr
tf = 0.15 # Gyr
pot = MW14
solarmotion = 'schoenrich' # or 'hogg' or 'dehnen', or value in [-U,V,W]

# Load data
star_id = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [0], dtype = "str")
ra = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [1]) # degrees
dec = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [2]) # degrees
d = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [3])/1000. # kpc
rv = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [4]) # km/s
pmra = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [5]) # mas/yr
pmde = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [6]) # mas/yr

Nstars = len(star_id)

# Convert radec to lb
l = []
b = []

for i in range(Nstars):
    c = SkyCoord(ra=ra[i]*u.degree, dec=dec[i]*u.degree, frame='icrs')
    l.append(c.galactic.l.value)
    b.append(c.galactic.b.value)

l = np.array(l)
b = np.array(b)

# Integrate orbits for all stars and plot
ts_fwd = np.linspace(0, tf, 10000)*u.Gyr
ts_bwd = np.linspace(0, ti, 10000)*u.Gyr

colors = ["#0000AA", "#AA0000", "#00AA00", "#AA00AA", "#AAAA00", "#00AAAA", "#AA6600", "#AA0066", "#00AA66", "#66AA00", "#0066AA", "#6600AA", "#666666", "#000000", "#000000", "#000000", "#000000", "#000000", "#000000"]
radconv = np.pi/180.

def cut_l(l, b):
    for i in range(len(l)-1):
        if np.sign(l[i+1]) != np.sign(l[i]):
            l[i] = np.nan
            b[i] = np.nan
    
    return l, b 

def plot_lb(i, orbit_fwd, orbit_bwd, alpha = 0.5):
    l_fwd = orbit_fwd.ll(ts_fwd) # deg
    b_fwd = orbit_fwd.bb(ts_fwd) # deg
    l_bwd = orbit_bwd.ll(ts_bwd) # deg
    b_bwd = orbit_bwd.bb(ts_bwd) # deg
    
    l_fwd[l_fwd > 180] = l_fwd[l_fwd > 180] - 360
    l_bwd[l_bwd > 180] = l_bwd[l_bwd > 180] - 360
    
    l_fwd, b_fwd = cut_l(l_fwd, b_fwd)
    l_bwd, b_bwd = cut_l(l_bwd, b_bwd)
    
    li = orbit_bwd.ll(ti*u.Gyr)
    bi = orbit_bwd.bb(ti*u.Gyr) 
    lf = orbit_fwd.ll(tf*u.Gyr)
    bf = orbit_fwd.bb(tf*u.Gyr)
    
    if li > 180: li = li -360
    if lf > 180: lf = lf -360
    
    plt.plot(l_fwd*radconv, b_fwd*radconv, '-', color = colors[i], linewidth = 1, alpha = alpha)
    plt.plot(l_bwd*radconv, b_bwd*radconv, '--', color = colors[i], linewidth = 1, alpha = alpha)
    plt.plot(l[i]*radconv, b[i]*radconv, marker = 'o', markeredgecolor = colors[i], color = colors[i], alpha = alpha)
    plt.plot(li*radconv, bi*radconv, marker = 'o', markersize = 4, markeredgecolor = colors[i], color = "#FFFFFF", alpha = alpha)
    plt.plot(lf*radconv, bf*radconv, marker = 'o', markersize = 4, markeredgecolor = colors[i], color = "#000000", alpha = alpha)

    
plt.figure()
plt.subplot(111, projection="aitoff")

for i in range(Nstars):  
    print "Integrating orbit for {0}".format(star_id[i])
    
    # prepare data
    vxvv = [ra[i], dec[i], d[i], pmra[i], pmde[i], rv[i]]
    orbit_fwd = Orbit(vxvv=vxvv, radec = True, ro = ro, vo = vo, solarmotion = solarmotion)
    orbit_fwd.integrate(ts_fwd, pot, method = 'leapfrog')
    orbit_bwd = Orbit(vxvv=vxvv, radec = True, ro = ro, vo = vo, solarmotion = solarmotion)
    orbit_bwd.integrate(ts_bwd, pot, method = 'leapfrog')
    
    plot_lb(i, orbit_fwd, orbit_bwd)

plt.show()

i = 4

plt.figure()
plt.subplot(111, projection="aitoff")

vxvv = [ra[i], dec[i], d[i], pmra[i], pmde[i], rv[i]]
orbit_fwd = Orbit(vxvv=vxvv, radec = True, ro = ro, vo = vo, solarmotion = solarmotion)
orbit_fwd.integrate(ts_fwd, pot, method = 'leapfrog')
orbit_bwd = Orbit(vxvv=vxvv, radec = True, ro = ro, vo = vo, solarmotion = solarmotion)
orbit_bwd.integrate(ts_bwd, pot, method = 'leapfrog')
    
plot_lb(i, orbit_fwd, orbit_bwd)

plt.show()

#
#

#i = 0
#vxvv = [ra[i], dec[i], d[i], pmra[i], pmde[i], rv[i]]
#orbit_fwd = Orbit(vxvv=vxvv, radec = True, ro = ro, vo = vo, solarmotion = solarmotion)
#orbit_fwd.integrate(ts_fwd, pot, method = 'leapfrog')
#orbit_bwd = Orbit(vxvv=vxvv, radec = True, ro = ro, vo = vo, solarmotion = solarmotion)
#orbit_bwd.integrate(ts_bwd, pot, method = 'leapfrog')

#l_fwd = orbit_fwd.ll(ts_fwd) # deg
#b_fwd = orbit_fwd.bb(ts_fwd) # deg
#l_bwd = orbit_bwd.ll(ts_bwd) # deg
#b_bwd = orbit_bwd.bb(ts_bwd) # deg

#li = orbit_bwd.ll(ti*u.Gyr)
#bi = orbit_bwd.bb(ti*u.Gyr) 
#lf = orbit_fwd.ll(tf*u.Gyr)
#bf = orbit_fwd.bb(tf*u.Gyr)
#radconv = np.pi/180.

#plt.figure()
#plt.subplot(111, projection="aitoff")
#color="#0000AA"
#plt.plot(ll_fwd*radconv, bb_fwd*radconv, '-', color = color, linewidth = 1)
#plt.plot(ll_bwd*radconv, bb_bwd*radconv, '--', color = color, linewidth = 1)
#plt.plot(l[i]*radconv, b[i]*radconv, marker = 'o', markeredgecolor = color, color = color)
#plt.plot(li*radconv, bi*radconv, marker = 'o', markersize = 4, markeredgecolor = color, color = "#FFFFFF")
#plt.plot(lf*radconv, bf*radconv, marker = 'o', markersize = 4, markeredgecolor = color, color = "#000000")
#plt.show()

#        
