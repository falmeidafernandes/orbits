# Import libraries
import numpy as np
from matplotlib import pyplot as plt

import galpy
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014 as MW14
import galpy.util.bovy_coords as coords


# Galactic model parameters
ro = 8.0  # kpc
vo = 220.0  # km/s
ts = np.linspace(0, 200, 400000) # Integration time. Translates to approx 7 Gyr
pot = MW14


# Solar Peculiar Velocity
Usun = 9.8 # km/s
Vsun = 12.5 # km/s
Wsun = 7.2 # km/s


# Calculate galpy coordinates vxvv
def get_vxvv(U, V, W, X, Y, Z): 
    # Calculate uvw
    u = U + Usun
    v = V + Vsun
    w = W + Wsun
    
    # Calculate vxvv
    Rgal = np.sqrt(ro**2 + X**2 + Y**2 + 2*ro*X)
    R = Rgal/ro
    phi = np.arctan2(Y, ro+X)*(180/np.pi) # degrees
    z = Z/ro
    vR = u/vo
    vT = (vo + v)/vo
    vz = w/vo
        
    vxvv = [R, vR, vT, z, vz, phi]
    return vxvv
    
# Load data
star_id = np.loadtxt("porbitas_UVWXYZ.dat", delimiter = ',', usecols = [0], dtype = "str")
ra = np.loadtxt("porbitas_UVWXYZ.dat", delimiter = ',', usecols = [1]) # degrees
dec = np.loadtxt("porbitas_UVWXYZ.dat", delimiter = ',', usecols = [2]) # degrees
d = np.loadtxt("porbitas_UVWXYZ.dat", delimiter = ',', usecols = [3]) # pc
rv = np.loadtxt("porbitas_UVWXYZ.dat", delimiter = ',', usecols = [4]) # km/s
pmra = np.loadtxt("porbitas_UVWXYZ.dat", delimiter = ',', usecols = [5]) # mas/yr
pmde = np.loadtxt("porbitas_UVWXYZ.dat", delimiter = ',', usecols = [6]) # mas/yr
U = np.loadtxt("porbitas_UVWXYZ.dat", delimiter = ',', usecols = [7]) # degrees
V = np.loadtxt("porbitas_UVWXYZ.dat", delimiter = ',', usecols = [8]) # degrees
W = np.loadtxt("porbitas_UVWXYZ.dat", delimiter = ',', usecols = [9]) # pc
X = np.loadtxt("porbitas_UVWXYZ.dat", delimiter = ',', usecols = [10]) # km/s
Y = np.loadtxt("porbitas_UVWXYZ.dat", delimiter = ',', usecols = [11]) # mas/yr
Z = np.loadtxt("porbitas_UVWXYZ.dat", delimiter = ',', usecols = [12]) # mas/yr

Nstars = len(star_id)

#Iterate orbit integration for every star
for i in range(Nstars):  
    print "Integrating orbit for {0}".format(star_id[i])

    # prepare data
    vxvv_calcUVW = get_vxvv(U=U[i], V=V[i], W=W[i], X=X[i], Y=Y[i], Z=Z[i])
    orbit_calcUVW = Orbit(vxvv=vxvv_calcUVW)
    
    # vxvv -> [ra[deg], dec[deg], d[kpc], pmra[mas/yr], pmde[mas/yr], rv[km/s]]
    orbit_obs = Orbit(vxvv=[ra[i], dec[i], d[i]/1000., pmra[i], pmde[i], rv[i]], radec = True, ro = ro, vo = vo)
    
    print orbit_calcUVW._orb.vxvv
    print orbit_obs._orb.vxvv
