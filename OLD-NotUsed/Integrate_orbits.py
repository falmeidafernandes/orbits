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


# Prepare objects to save calculated data
Nstars = len(star_id)
e = np.empty(Nstars); e[:] = np.nan
Rmax = np.empty(Nstars); Rmax[:] = np.nan
Rmin = np.empty(Nstars); Rmin[:] = np.nan
zmax = np.empty(Nstars); zmax[:] = np.nan

#Iterate orbit integration for every star
for i in range(Nstars):  
    print "Integrating orbit for {0}".format(star_id[i])

    # prepare data
    vxvv = get_vxvv(U=U[i], V=V[i], W=W[i], X=X[i], Y=Y[i], Z=Z[i])
    orbit = Orbit(vxvv=vxvv)
    orbit.integrate(ts, pot, method = 'leapfrog')
    
    # saving orbit data for plot
    x = orbit.x(ts)*ro # kpc
    y = orbit.y(ts)*ro # kpc
    z = orbit.z(ts)*ro # kpc
    r = np.sqrt(x*x+y*y) # kpc
    
    np.save("./Orbits/{}/x.npy".format(star_id[i]), x)
    np.save("./Orbits/{}/y.npy".format(star_id[i]), y)
    np.save("./Orbits/{}/z.npy".format(star_id[i]), z)
    np.save("./Orbits/{}/r.npy".format(star_id[i]), r)
    
    e[i] = orbit.e()
    Rmax[i] = orbit.rap()*ro # kpc
    Rmin[i] = orbit.rperi()*ro # kpc
    zmax[i] = orbit.zmax()*ro # kpc


# Save data
with open("porbitas_full.dat", "w") as f:
    # Save header
    f.write("#Estrela,ra,dec,dist,RV,pmra,pmdec,U,V,W,X,Y,Z,e,zmax,Rmax,Rmin\n")
    
    # Iterate to save data for each star
    for i in range(Nstars):
        f.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(star_id[i],
                                                                              ra[i],
                                                                              dec[i],
                                                                              d[i],
                                                                              rv[i],
                                                                              pmra[i],
                                                                              pmde[i],
                                                                              U[i],
                                                                              V[i],
                                                                              W[i],
                                                                              X[i],
                                                                              Y[i],
                                                                              Z[i],
                                                                              e[i],
                                                                              zmax[i],
                                                                              Rmax[i],
                                                                              Rmin[i]))

        
