    # Import libraries
import numpy as np
import galpy
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014 as MW14
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib import pyplot as plt

RUN = True # True, integrates orbit and saves results; False, loads results

# Galactic model parameters
ro = 8.0  # kpc
vo = 220.0  # km/s
ti = -0.10 # Gyr
tf = 0.10 # Gyr
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

colors = ["#000000", "#FF0000", "#FF0000", "#000000", "#FF0000", "#FF0000", "#000000", "#000000", "#FF0000", "#000000", "#FF0000", "#FF0000", "#000000"]
zorders = [2, 4, 4, 2, 4, 4, 2, 2, 4, 2, 4, 4, 2]
radconv = np.pi/180.


def plot_lb(i, ts, plot_extreme = False, save_name = 'fwd', zorder = 1, **kargs):
    # integrate orbit
    if RUN:
        vxvv = [ra[i], dec[i], d[i], pmra[i], pmde[i], rv[i]]
        o = Orbit(vxvv=vxvv, radec = True, ro = ro, vo = vo, solarmotion = solarmotion)
        o.integrate(ts, pot, method = 'leapfrog')
                    
        # get orbit l and b
        ll = o.ll(ts)
        bb = o.bb(ts)
        
        np.save("./Orbits/{}/l{}150Myr.npy".format(star_id[i], save_name), ll)
        np.save("./Orbits/{}/b{}150Myr.npy".format(star_id[i], save_name), bb)
    
    else:
        ll = np.load("./Orbits/{}/l{}150Myr.npy".format(star_id[i], save_name))
        bb = np.load("./Orbits/{}/b{}150Myr.npy".format(star_id[i], save_name))
    
    Npoints = len(ll)
    
    # Find boundary breaks at -180 and 180
    boundary_breaks = []
    
    # crossings for 180
    for j in range(1, Npoints):
        # crossings for 180+
        if (ll[j-1] <= 180) and (ll[j] > 180):
            boundary_breaks.append(j)
        # crossings for 180-
        elif (ll[j-1] >= 180) and (ll[j] < 180):
            boundary_breaks.append(j)
    
    boundary_breaks.append(Npoints)
    
    # fix ll > 180
    ll[ll > 180] = ll[ll > 180] - 360
    
    # plots
    j0 = 0
    for k in range(len(boundary_breaks)):
        jf = boundary_breaks[k]
        plt.plot(ll[j0:jf]*radconv, bb[j0:jf]*radconv, zorder = zorder, **kargs)
        j0 = jf   
    
    # plot extremes
    if plot_extreme == 'lower':
        lower = plt.plot(ll[-1]*radconv, bb[-1]*radconv, marker = 'x', markersize = 4, color = colors[i], alpha = 1, zorder = zorder + 1)
        return lower
        
    elif plot_extreme == 'upper':
        upper = plt.plot(ll[-1]*radconv, bb[-1]*radconv, marker = 'o', markersize = 4, markeredgecolor = colors[i], color = "#FFFFFF", alpha = 1, zorder = zorder + 1)
        return upper

def plot_grid(alpha1 = 0.2, alpha2 = 0.1, N = 1000, **kargs):
    lgrid = np.linspace(-180, 180, N)*radconv
    bgrid = np.linspace(-90, 90, N)*radconv
    
    for l, b in zip([-120, -60, 0, 60, 120], [-60, -30, 0, 30, 60]):
        plt.plot([l*radconv]*N, bgrid, alpha = alpha1, **kargs)
        plt.plot(lgrid, [b*radconv]*N, alpha = alpha1, **kargs)
    
    for l, b in zip([-150, -90, -30, 30,  90, 150], [-75, -45, -15, 15, 45, 75]):
        plt.plot([l*radconv]*N, bgrid, alpha = alpha2, **kargs)
        plt.plot(lgrid, [b*radconv]*N, alpha = alpha2, **kargs)
        

# Individual plots
if 0:
    for i in range(13):  
        print "Integrating and plotting orbit for {0}".format(star_id[i])
        plt.figure(figsize = (12,6.3))
        plt.subplot(111, projection="aitoff")
        plot_grid(alpha1 = 0.2, alpha2 = 0.2, color = '#666666', lw = 1, zorder = 0)
        current = plt.plot(l[i]*radconv, b[i]*radconv, marker = 'o', markeredgecolor = colors[i], color = colors[i], alpha = 1, zorder = zorders[i]+1)
        bwd = plot_lb(i, ts_bwd, plot_extreme = 'lower', save_name = 'bwd', color = colors[i], zorder = zorders[i], alpha = 0.5, ls = '--', lw = 2)
        fwd = plot_lb(i, ts_fwd, plot_extreme = 'upper', save_name = 'fwd', color = colors[i], zorder = zorders[i], alpha = 0.5, ls = '-', lw = 2)

        plt.tick_params(labelbottom=False, labelleft=False)
        plt.subplots_adjust(bottom = 0.01, top = 0.99, right = 0.99, left = 0.01)
        plt.savefig("./Orbits_plots/orbit_{0}_lb.png".format(star_id[i]))
        plt.savefig("./Orbits/{0}/orbit_{0}_lb.pdf".format(star_id[i]))
        plt.close()

# Plot all
if 1:
    plt.figure(figsize = (12,6.3))
    plt.subplot(111, projection="aitoff")
    plot_grid(alpha1 = 0.2, alpha2 = 0.2, color = '#666666', lw = 1, zorder = 0)
    
    for i in range(13):  
        print "Integrating and plotting orbit for {0}".format(star_id[i])
        current = plt.plot(l[i]*radconv, b[i]*radconv, marker = 'o', markeredgecolor = colors[i], color = colors[i], alpha = 1, zorder = zorders[i])
        bwd = plot_lb(i, ts_bwd, plot_extreme = 'lower', save_name = 'bwd', zorder = zorders[i], color = colors[i], alpha = 0.5, ls = '--', lw = 1)
        fwd = plot_lb(i, ts_fwd, plot_extreme = 'upper', save_name = 'fwd', zorder = zorders[i], color = colors[i], alpha = 0.5, ls = '-', lw = 1)

    plt.tick_params(labelbottom=False, labelleft=False)
    plt.subplots_adjust(bottom = 0.01, top = 0.99, right = 0.99, left = 0.01)
    plt.savefig("lb_plot_100Myr.png")
    #plt.savefig("lb_plot.pdf")
    plt.close()


