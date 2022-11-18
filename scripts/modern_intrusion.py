import numpy as np
import matplotlib.pyplot as plt

locs = ["Albany", "Aviero", "Bunbury", "Canterbury", "Carnavon", "Nantucket Is", "New Jersey", "Niger Delta", "Northern Florida", "Shanghai", "Suriname"]
Ks = [5, 1, 20, 3, 11, 10, 15, 5, 42.5, 2.6, 25]
Kls = [0.005, 0.0001, 1e-6, 0.0001, 0.0001,0.0095, 0.0003, 0.1, 0.0003, 0.004, 0.0001]
Ds = [1, 80, 20, 25, 2, 34, 30, 15, 30, 80, 30]
Ls = [1000*L for L in [42, 50, 125, 50, 65, 200, 140, 60, 110, 500, 130]]
zs = [4, 130, 80, 35, 8, 36, 214, 110, 120, 80, 70] # dunno about shanghai
x_bs = [1000*x_b for x_b in [8, 1, 4.4, 1.9, 6.1, 2, 3, 25, 5, 25, 5]]
x_tipss = [[frac*L for frac in np.linspace(0.05, 0.98, 40)] for L in Ls]
betas = [np.arctan(z/L) for z, L in zip(zs, Ls)]
Hs = [20, 200, 300, 16, 45, 220, 36, 200,  42.5, 2.6, 25]
z0s = [z+D+H for z, D, H in zip(zs, Ds, Hs)]
h_bs = [h for h in [0.6, 5, 1.6, -1.5, 1.1, 4, -24, 9, 0, -7, -10]]


def v_adv(h_b, h_f,  x_b, L, x_tip, K, n):
    return (h_b-h_f)/(x_b+L-x_tip)*K/n

def h_f(z0, x_tip, beta, H):
    return (1025-1000)/1000*(z0-H-x_tip*np.tan(beta))

h_fs = [h_f(z0, 0, beta, H) for z0, beta, H in  zip(z0s, betas, Hs)]
h_driving = [h_b - h_f for h_b, h_f in zip(h_bs, h_fs)]

def main():
    f, ax = plt.subplots()
    for loc, K, L, z0, x_b, x_tips, beta, H, h_b in zip(locs, Ks, Ls, z0s, x_bs, x_tipss, betas, Hs, h_bs):
        if loc in ["Aviero", "Bunbury", "Canterbury", "New Jersey", "Northern Florida", "Shanghai", "Suriname"]:
            ax.plot([x_tip/L for x_tip in x_tips],[-365*100*v_adv(h_b, h_f(z0, x_tip, beta, H), x_b, L, x_tip, K, 0.35) for x_tip in x_tips], label=loc)

    ax.set_xlabel("x_tip/L [-]")
    ax.set_ylabel("V_tip [log(m/100yr)]")
    ax.legend()
    ax.set_yscale("log")
    plt.show()

if __name__=="__main__":
    main()