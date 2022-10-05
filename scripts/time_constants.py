import matplotlib.pyplot as plt
import numpy as np

def tau_c(Ss, L, Kh):
    return Ss*L**2/Kh

def t_ne(tau_c):
    return 3*4/np.pi**2*tau_c

def main():
    Sss = [1e-3, 1e-4, 1e-5]
    labels = ["1e-3", "1e-4", "1e-5"]
    Khs = np.linspace(0.1, 100, 50)
    L = 50000
    f, ax = plt.subplots()
    for Ss, label in zip(Sss, labels):
        t_nes = [t_ne(tau_c(Ss, L, Kh))/365/1000 for Kh in Khs]
        ax.plot(Khs, t_nes, label=label)

    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_ylabel("Time to near equillibrium [ky]")
    ax.set_xlabel("Hydraulic conductivity [m/day]")
    ax.set_title("Time to near equillibrium in 50km long gravel/sand aquifers")
    ax.legend(title="Ss")

    plt.show()
    

if __name__=="__main__":
    main()