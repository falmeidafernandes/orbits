# Import libraries
import numpy as np
from matplotlib import pyplot as plt

# Load data
ra = np.loadtxt("GCSporbitas_orbparams.dat", delimiter = ',', usecols = [1]) # degrees
dec = np.loadtxt("GCSporbitas_orbparams.dat", delimiter = ',', usecols = [2]) # degrees
d = np.loadtxt("GCSporbitas_orbparams.dat", delimiter = ',', usecols = [3]) # kpc
pmra = np.loadtxt("GCSporbitas_orbparams.dat", delimiter = ',', usecols = [4]) # mas/yr
pmde = np.loadtxt("GCSporbitas_orbparams.dat", delimiter = ',', usecols = [5]) # mas/yr
rv = np.loadtxt("GCSporbitas_orbparams.dat", delimiter = ',', usecols = [6]) # km/s

e = np.loadtxt("GCSporbitas_orbparams.dat", delimiter = ',', usecols = [7])
eGCS = np.loadtxt("GCSporbitas_orbparams.dat", delimiter = ',', usecols = [8])
zmax = np.loadtxt("GCSporbitas_orbparams.dat", delimiter = ',', usecols = [9]) # kpc
zmaxGCS = np.loadtxt("GCSporbitas_orbparams.dat", delimiter = ',', usecols = [10]) # kpc
Rmax = np.loadtxt("GCSporbitas_orbparams.dat", delimiter = ',', usecols = [11]) # kpc
RmaxGCS = np.loadtxt("GCSporbitas_orbparams.dat", delimiter = ',', usecols = [12]) # kpc
Rmin = np.loadtxt("GCSporbitas_orbparams.dat", delimiter = ',', usecols = [13]) # kpc
RminGCS = np.loadtxt("GCSporbitas_orbparams.dat", delimiter = ',', usecols = [14]) # kpc


plt.plot(eGCS, e, 'ro')
plt.gca().set_xlabel('eGCS')
plt.gca().set_ylabel('e')
plt.savefig("comp_e.png")
plt.clf()
plt.close()

plt.plot(zmaxGCS, zmax, 'ro')
plt.gca().set_xlabel('zmaxGCS')
plt.gca().set_ylabel('zmax')
plt.savefig("comp_zmax.png")
plt.clf()
plt.close()

plt.plot(RmaxGCS, Rmax, 'ro')
plt.gca().set_xlabel('RmaxGCS')
plt.gca().set_ylabel('Rmax')
plt.savefig("comp_Rmax.png")
plt.clf()
plt.close()

plt.plot(RminGCS, Rmin, 'ro')
plt.gca().set_xlabel('RminGCS')
plt.gca().set_ylabel('Rmin')
plt.savefig("comp_Rmin.png")
plt.clf()
plt.close()
