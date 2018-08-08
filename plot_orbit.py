import numpy as np
from matplotlib import pyplot as plt

star_ids = np.loadtxt("porbitas.dat", delimiter = ',', usecols = [0], dtype = "str")
colors = ["#000000", "#FF0000", "#FF0000", "#000000", "#FF0000", "#FF0000", "#000000", "#000000", "#FF0000", "#000000", "#FF0000", "#FF0000", "#000000", "#0000AA", "#0000AA", "#0000AA", "#0000AA", "#0000AA"]

for i in range(len(star_ids)):
    star_id = star_ids[i]
    color = colors[i]
    
    x = np.load("./Orbits/{}/x.npy".format(star_id))
    y = np.load("./Orbits/{}/y.npy".format(star_id))
    z = np.load("./Orbits/{}/z.npy".format(star_id))
    r = np.load("./Orbits/{}/r.npy".format(star_id))

    plt.plot(x, y, color = color)
    plt.gca().set_xlim((np.min((x.min(), y.min())), np.max((x.max(), y.max()))))
    plt.gca().set_ylim((np.min((x.min(), y.min())), np.max((x.max(), y.max()))))
    plt.gca().set_xlabel("x (kpc)")
    plt.gca().set_ylabel("y (kpc)")
    #plt.savefig("./Orbits/{0}/orbit_{0}_xy.png".format(star_id))
    plt.savefig("./Orbits_plots/orbit_{0}_xy.png".format(star_id))
    plt.close()

    plt.plot(x, z, color = color)
    plt.gca().set_xlim((np.min((x.min(), z.min())), np.max((x.max(), z.max()))))
    plt.gca().set_ylim((np.min((x.min(), z.min())), np.max((x.max(), z.max()))))
    plt.gca().set_xlabel("x (kpc)")
    plt.gca().set_ylabel("z (kpc)")
    #plt.savefig("./Orbits/{0}/orbit_{0}_xz.png".format(star_id))
    plt.savefig("./Orbits_plots/orbit_{0}_xz.png".format(star_id))
    plt.close()

