# Import libraries
import numpy as np
import galpy
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014 as MW14

# Galactic model parameters
ro = 8.0  # kpc
vo = 220.0  # km/s
ts = np.linspace(0, 200, 400000) # Integration time. Translates to approx 7 Gyr
pot = MW14
solarmotion = 'schoenrich' # or 'hogg' or 'dehnen', or value in [-U,V,W]

# Load data
ra = np.loadtxt("GCSporbitas.dat", delimiter = '|', usecols = [0]) # degrees
dec = np.loadtxt("GCSporbitas.dat", delimiter = '|', usecols = [1]) # degrees
plx = np.loadtxt("GCSporbitas.dat", delimiter = '|', usecols = [2]) # mas
plx = plx/1000. # arcsec
d = 1./plx #pc
d = d/1000. #kpc

pmra = np.loadtxt("GCSporbitas.dat", delimiter = '|', usecols = [3]) # mas/yr
pmde = np.loadtxt("GCSporbitas.dat", delimiter = '|', usecols = [4]) # mas/yr
rv = np.loadtxt("GCSporbitas.dat", delimiter = '|', usecols = [5]) # km/s

RminGCS = np.loadtxt("GCSporbitas.dat", delimiter = '|', usecols = [6]) # kpc
RmaxGCS = np.loadtxt("GCSporbitas.dat", delimiter = '|', usecols = [7]) # kpc
eGCS = np.loadtxt("GCSporbitas.dat", delimiter = '|', usecols = [8])
zmaxGCS = np.loadtxt("GCSporbitas.dat", delimiter = '|', usecols = [9]) # kpc

Nstars = 10

# Prepare objects to save calculated data
e = np.empty(Nstars); e[:] = np.nan
Rmax = np.empty(Nstars); Rmax[:] = np.nan
Rmin = np.empty(Nstars); Rmin[:] = np.nan
zmax = np.empty(Nstars); zmax[:] = np.nan

# Integrate orbits for all stars
for i in range(Nstars):  
    print "Integrating orbit for Star_{0}".format(i)

    # prepare data
    vxvv = [ra[i], dec[i], d[i], pmra[i], pmde[i], rv[i]]
    orbit = Orbit(vxvv=vxvv, radec = True, ro = ro, vo = vo, solarmotion = solarmotion)
    orbit.integrate(ts, pot, method = 'leapfrog')
    
    # saving orbit data for plot
    x = orbit.x(ts) # kpc
    y = orbit.y(ts) # kpc
    z = orbit.z(ts) # kpc
    r = np.sqrt(x*x+y*y) # kpc
    
    np.save("./Orbits/Star_{}/x.npy".format(i), x)
    np.save("./Orbits/Star_{}/y.npy".format(i), y)
    np.save("./Orbits/Star_{}/z.npy".format(i), z)
    np.save("./Orbits/Star_{}/r.npy".format(i), r)
    
    e[i] = orbit.e()
    Rmax[i] = orbit.rap() # kpc
    Rmin[i] = orbit.rperi() # kpc
    zmax[i] = orbit.zmax() # kpc
    
# Save results
with open("GCPporbitas_orbparams.dat", "w") as f:
    # Save header
    f.write("#Estrela,ra,dec,dist,RV,pmra,pmdec,eGCS,e,zmaxGCS,zmax,RmaxGCS,Rmax,RminGCS,Rmin\n")
    
    # Iterate to save data for each star
    for i in range(Nstars):
        f.write("Star_{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(i,
                                                                        ra[i],
                                                                        dec[i],
                                                                        d[i],
                                                                        rv[i],
                                                                        pmra[i],
                                                                        pmde[i],
                                                                        eGCS[i],
                                                                        e[i],
                                                                        zmaxGCS[i],
                                                                        zmax[i],
                                                                        RmaxGCS[i],
                                                                        Rmax[i],
                                                                        RminGCS[i],
                                                                        Rmin[i]))

        
