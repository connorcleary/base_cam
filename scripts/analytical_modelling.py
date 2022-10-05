import numpy as np
import matplotlib.pyplot as plt

# parameters
k = 3 # hydraulic conductivity, m/d
zt = -60 # top of aquifer, m
zb = -76 # bottom of aquifer, m
alpha = 40 # alpha factor, -
clist = [250000] # three values of resistance of leaky layer, d
U = 0.003 # flow toward the coast, m^2/d
H = zt - zb # aquifer thickness, m

# solution
def interface(k, zt, zb, c, alpha, U):
    H = zt - zb
    # below sea
    hstar = -zt / alpha
    xtip = (18 * U * k * alpha * c**2) ** (1 / 3)
    xsea = np.linspace(0, xtip, 100)
    hsea = hstar + (xsea - xtip)**2 / (6 * k * alpha * c)
    h0 = hstar + xtip ** 2 / (6 * k * alpha * c)
    # below land
    phi0 = 0.5 * k * alpha * (h0 - hstar)**2
    phitoe = 0.5 * k * H ** 2 / alpha
    xtoe = -(phitoe - phi0) / U
    Cc = (0.5 * k * H ** 2 + k * H * zb) / alpha
    xland = np.linspace(xtoe, 0, 100)
    phi = -U * xland + phi0
    hland = np.sqrt(2 * phi/ (k * alpha)) + hstar
    # combine solution
    x = np.hstack((xland, xsea))
    zi = -alpha * np.hstack((hland, hsea))
    return x, zi

def main():
    # basic plot
    ax = plt.subplot(111, aspect=4, ylim=(-76, -60))
    for c in clist:
        x, zi = interface(k=k, zt=zt, zb=zb, c=c, alpha=alpha, U=U)
        plt.plot(x, zi, label=f'c={c:.0f} d')

    ax.set_aspect(100)
    plt.legend()
    plt.show()


if __name__=="__main__":
    main()