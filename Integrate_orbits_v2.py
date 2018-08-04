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
star_id = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [0], dtype = "str")
ra = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [1]) # degrees
dec = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [2]) # degrees
d = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [3])/1000. # kpc
rv = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [4]) # km/s
pmra = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [5]) # mas/yr
pmde = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [6]) # mas/yr

Nstars = len(star_id)

# Prepare objects to save calculated data
e = np.empty(Nstars); e[:] = np.nan
Rmax = np.empty(Nstars); Rmax[:] = np.nan
Rmin = np.empty(Nstars); Rmin[:] = np.nan
zmax = np.empty(Nstars); zmax[:] = np.nan

# Integrate orbits for all stars
for i in range(Nstars):  
    print "Integrating orbit for {0}".format(star_id[i])

    # prepare data
    vxvv = [ra[i], dec[i], d[i], pmra[i], pmde[i], rv[i]]
    orbit = Orbit(vxvv=vxvv, radec = True, ro = ro, vo = vo, solarmotion = solarmotion)
    orbit.integrate(ts, pot, method = 'leapfrog')
    
    # saving orbit data for plot
    x = orbit.x(ts) # kpc
    y = orbit.y(ts) # kpc
    z = orbit.z(ts) # kpc
    r = np.sqrt(x*x+y*y) # kpc
    
    np.save("./Orbits/{}/x.npy".format(star_id[i]), x)
    np.save("./Orbits/{}/y.npy".format(star_id[i]), y)
    np.save("./Orbits/{}/z.npy".format(star_id[i]), z)
    np.save("./Orbits/{}/r.npy".format(star_id[i]), r)
    
    e[i] = orbit.e()
    Rmax[i] = orbit.rap() # kpc
    Rmin[i] = orbit.rperi() # kpc
    zmax[i] = orbit.zmax() # kpc
    
# Save results
with open("porbitas_orbparams.dat", "w") as f:
    # Save header
    f.write("#Estrela,ra,dec,dist,RV,pmra,pmdec,e,zmax,Rmax,Rmin\n")
    
    # Iterate to save data for each star
    for i in range(Nstars):
        f.write("{},{},{},{},{},{},{},{},{},{},{}\n".format(star_id[i],
                                                            ra[i],
                                                            dec[i],
                                                            d[i],
                                                            rv[i],
                                                            pmra[i],
                                                            pmde[i],
                                                            e[i],
                                                            zmax[i],
                                                            Rmax[i],
                                                            Rmin[i]))

        
