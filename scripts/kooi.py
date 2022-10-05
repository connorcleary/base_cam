import matplotlib.pyplot as plt
import numpy as np

def v_adv(K, n, beta):
    rho_f = 1000
    rho_s = 1025
    return K*rho_s/rho_f/n*np.tan(beta)

def main():
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

if __name__=="__main__":
    main()