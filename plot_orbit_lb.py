import numpy as np
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord

star_ids = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [0], dtype = "str")
ra = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [1])
dec = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [2])
c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')

for star_id in star_ids:
    l = np.load("./Orbits_lb/{}/l.npy".format(star_id))
    b = np.load("./Orbits_lb/{}/b.npy".format(star_id))
    ts = np.load("./Orbits_lb/{}/ts.npy".format(star_id))

    plt.scatter(l, b)
    plt.gca().set_xlim((np.min((x.min(), y.min())), np.max((x.max(), y.max()))))
    plt.gca().set_ylim((np.min((x.min(), y.min())), np.max((x.max(), y.max()))))
    plt.gca().set_xlabel("x (kpc)")
    plt.gca().set_ylabel("y (kpc)")
    #plt.savefig("./Orbits/{0}/orbit_{0}_xy.png".format(star_id))
    plt.savefig("./Orbits_plots/orbit_{0}_xy.png".format(star_id))
    plt.close()

    plt.plot(x, z)
    plt.gca().set_xlim((np.min((x.min(), z.min())), np.max((x.max(), z.max()))))
    plt.gca().set_ylim((np.min((x.min(), z.min())), np.max((x.max(), z.max()))))
    plt.gca().set_xlabel("x (kpc)")
    plt.gca().set_ylabel("z (kpc)")
    #plt.savefig("./Orbits/{0}/orbit_{0}_xz.png".format(star_id))
    plt.savefig("./Orbits_plots/orbit_{0}_xz.png".format(star_id))
    plt.close()

radconv = np.pi/180.

star_id = star_ids[1]
l = np.load("./Orbits_lb/{}/l.npy".format(star_id))
# fix l > 180:
l[l > 180] = l[l > 180] - 360
b = np.load("./Orbits_lb/{}/b.npy".format(star_id))
ts = np.load("./Orbits_lb/{}/ts.npy".format(star_id))
plt.figure()
plt.subplot(111, projection="aitoff")
plt.plot(l*radconv, b*radconv)
plt.scatter(112.73112471*radconv, -28.01422132*radconv)
plt.show()
