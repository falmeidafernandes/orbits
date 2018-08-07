# Calculate the velocities UVW and galactic coordinates XYZ for a stellar sample
# *U points towards Galactic Center
# X = d cos(b) cos(l)
# Y = d cos(b) sin(l)
# Z = d sin(b)

# Import libraries
import numpy as np
from druvw import uvw
from druvw import xyz

# Load data
star_id = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [0], dtype = "str")
ra = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [1]) # degrees
dec = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [2]) # degrees
d = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [3]) # pc
rv = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [4]) # km/s
pmra = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [5]) # mas/yr
pmde = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [6]) # mas/yr

# Calculate number of stars
Nstars = len(star_id)

# Define vector to hold UVWXYZ for each star'
U = np.empty(Nstars); U[:] = np.nan
V = np.empty(Nstars); V[:] = np.nan
W = np.empty(Nstars); W[:] = np.nan
X = np.empty(Nstars); X[:] = np.nan
Y = np.empty(Nstars); Y[:] = np.nan
Z = np.empty(Nstars); Z[:] = np.nan

# Iterate to calculate UVWXYZ for each star

for i in range(Nstars):
    U[i], V[i], W[i] = uvw(ra = ra[i],
                           dec = dec[i],
                           d = d[i],
                           pmra = pmra[i],
                           pmde = pmde[i],
                           rv = rv[i])
                           
    X[i], Y[i], Z[i] = xyz(ra = ra[i],
                           dec = dec[i],
                           d = d[i])

# Save data
with open("porbitas_UVWXYZ.dat", "w") as f:
    # Save header
    f.write("#Estrela,ra,dec,dist,RV,pmra,pmdec,U,V,W,X,Y,Z\n")
    
    # Iterate to save data for each star
    for i in range(Nstars):
        f.write("{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(star_id[i],
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
                                                                Z[i]))
