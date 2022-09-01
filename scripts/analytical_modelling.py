import numpy as np
import matplotlib.pyplot as plt

# parameters
k = 10 # hydraulic conductivity, m/d
zt = -10 # top of aquifer, m
zb = -30 # bottom of aquifer, m
alpha = 40 # alpha factor, -
clist = [1e-12, 5, 50] # three values of resistance of leaky layer, d
U = 0.4 # flow toward the coast, m^2/d
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
    plt.subplot(111, aspect=4, ylim=(-30, -10))
    for c in clist:
        x, zi = interface(k=k, zt=zt, zb=zb, c=c, alpha=alpha, U=U)
        plt.plot(x, zi, label=f'c={c:.0f} d')

    plt.legend()
    plt.show()

    c = 5
    x, zi = interface(k=k, zt=zt, zb=zb, c=5, alpha=alpha, U=U)
    x = np.hstack((-200, x))
    zi = np.hstack((zb, zi))
    Qx = U * np.ones_like(x)
    lab = np.sqrt(c * k * H)
    xtip = (18 * U * k * alpha * c**2) ** (1 / 3)
    Qx[x>0] = -(x[x>0] - xtip) ** 3 / (18 * k * alpha * c**2)
    xg = [x, x]
    zg = [zi, zt * np.ones_like(x)]
    psi = [np.zeros_like(x), Qx]

    plt.subplot(111, aspect=2)
    plt.contour(xg, zg, psi)

    plt.show()

if __name__=="__main__":
    main()