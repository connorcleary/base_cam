import matplotlib.pyplot as plt
import numpy as np

def v_adv(K, n, beta):
    rho_f = 1000
    rho_s = 1025
    return K*rho_s/rho_f/n*np.tan(beta)

def Gamma(K, n, beta, v_slr):
    rho_f = 1000
    rho_s = 1025
    return n*v_slr*rho_f/(K*rho_s*np.tan(beta)**2)


def plot_advective_velocity():
    beta_shelf = np.tan(50/40000)
    beta_slope = None
    Kh_b = 0.1*365*1000
    n_b = 0.4
    beta = [np.arctan(z/40000) for z in [25, 50, 100]]
    v = v_adv(3, n_b, beta_shelf)
    f, ax = plt.subplots()
    for b in beta:
        v = [v_adv(K, n_b, b)/1000 for K in 365*1000*np.linspace(0.1, 100, 50)]
        ax.plot(np.linspace(0.1, 100, 50), v, label=f"{np.degrees(b):.2e}")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("velocity [km/ky]")
    ax.set_xlabel("K [m/day]")
    ax.set_title("Maximum advective velocity under SLR pressure perturbation")
    ax.legend(title="\u03B2 (\u00b0)")
    plt.show()

def main():
    locs = ["Albany", "Aviero", "Bunbury", "Canterbury", "Carnavon", "Nantucket Is", "New Jersey", "Niger Delta", "Northern Florida", "Shanghai", "Suriname"]
    Ks = [5, 1, 20, 3, 11, 10, 15, 5, 42.5, 2.6, 25]
    Ls = [1000*L for L in [42, 50, 125, 50, 65, 200, 140, 60, 110, 500, 130]]
    zs = [4, 130, 80, 35, 8, 36, 214, 110, 120, 80, 70] # dunno about shanghai
    lengths=1000*np.linspace(20,200,100) # m 
    Kh = [1, 10, 100] # m/d 
    betas = [np.arctan(120/length) for length in lengths]
    v_slr = 120/(20000*365) # m/day
    f, ax = plt.subplots()
    # for kh in Kh:
        # Gammas = Gamma(kh,  0.4, betas, v_slr)
        # ax.plot(lengths/1000, Gammas, label=f"{str(kh)}")
    for loc, K, L, z in zip(locs, Ks, Ls, zs):
        ax.vlines(x=L/1000, ymin=Gamma(K, 0.2, np.arctan(z/L),v_slr), ymax=Gamma(K, 0.5, np.arctan(z/L), v_slr), color="k")
        ax.text(s=loc, fontsize="x-small", x=L/1000, y=Gamma(0.8*K, 0.5, np.arctan(z/L), v_slr))
    ax.hlines(y=1, xmin=0, xmax=550, ls=":")
    ax.set_yscale("log")
    ax.set_ylabel("Non-dimensional transgression rate")
    ax.set_xlabel("Shelf length (km)")
    # ax.legend(title="K [m/day]")
    plt.show()
        

    
if __name__=="__main__":
    main()