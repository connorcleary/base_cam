import numpy as np
import matplotlib.pyplot as plt

def tau_adv(b, K, n):
    return 4*b*n*1000/(25*K)

def tau_diff(b, D):
    return b**2/D

def v_adv(K, n, beta):
    rho_f = 1000
    rho_s = 1025
    return K*rho_s/rho_f/n*np.tan(beta)

def Gamma(K, n, beta, v_slr):
    rho_f = 1000
    rho_s = 1025
    return n*v_slr*rho_f/(K*rho_s*np.tan(beta)**2)


def main():
    bs = np.linspace(0, 80, 50)
    Ks = [a/(24*60*60)/100 for a in [10**-4, 10**-3, 10**-2, 10**-1, 10**0, 10**1]]
    D = 3e-11
    n = 0.2

    locs = ["Albany", "Aviero", "Bunbury", "Canterbury", "Carnavon", "Nantucket Is", "New Jersey", "Niger Delta", "Northern Florida", "Shanghai", "Suriname"]
    Ks = [5, 1, 20, 3, 11, 10, 15, 5, 42.5, 2.6, 25]
    Kls = [0.005, 0.0001, 1e-6, 0.0001, 0.0001,0.0095, 0.0003, 0.1, 0.0003, 0.004, 0.0001]
    Ds = [1, 80, 20, 25, 2, 34, 30, 15, 30, 80, 30]
    Ls = [1000*L for L in [42, 50, 125, 50, 65, 200, 140, 60, 110, 500, 130]]
    zs = [4, 130, 80, 35, 8, 36, 214, 110, 120, 80, 70] # dunno about shanghai

    labels = ['10^-6', '10^-5', '10^-4', '10^-3', '10^-2', '10^-1']
    diff = [1/(60*60*24*365*1000)*tau_diff(b, D) for b in bs]
    lengths=1000*np.linspace(20,200,100)

    f,ax = plt.subplots()
    # ax.plot(bs, diff, label="diff", color="r", ls=":")
    v_slr = 120/(20000*365)

    # for K, label in zip(Ks, labels):
    #     adv = [1/(60*60*24*365*1000)*tau_adv(b, K, n) for b in bs]
    #     ax.plot(bs, adv, label=label)

    for loc, Kl, D, K, L, z in zip(locs, Kls, Ds, Ks, Ls, zs):
        ax.vlines(x=Gamma(K, 0.35, np.arctan(z/L),v_slr), ymax=tau_adv(D, Kl, 0.8)/365/1000/10, ymin=tau_adv(D, Kl, 0.2)/365/1000/10, color="k")
        ax.hlines(y=tau_adv(D, Kl, 0.5)/365/1000/10, xmin=Gamma(K, 0.2, np.arctan(z/L),v_slr), xmax=Gamma(K, 0.5, np.arctan(z/L), v_slr), color="k")
        ax.text(y=tau_adv(D, Kl, 0.8)/365/1000/10, x=Gamma(K, 0.5, np.arctan(z/L),v_slr), s=loc, fontsize="x-small")

    ax.hlines(y=1, xmin=0, xmax=1000, color="r", ls=":")
    ax.vlines(x=1, ymin=0, ymax=1000, color="r", ls=":")
    # ax.legend(title="Kz [m/day]", loc="upper right")
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel("T_adv/10ka")
    ax.set_xlabel("v_slr/v_adv")
    
    plt.show()
        

if __name__=="__main__":
    main()