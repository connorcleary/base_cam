import numpy as np
import matplotlib.pyplot as plt

def tau_adv(b, K, n):
    return 4*b*n*1000/(25*K)

def tau_diff(b, D):
    return b**2/D

def main():
    bs = np.linspace(0, 50, 50)
    Ks = [a/(24*60*60)/10 for a in [10**-4, 10**-3, 10**-2, 10**-1, 10**0, 10**1]]
    D = 3e-11
    n = 0.2

    labels = ['10^-4', '10^-3', '10^-2', '10^-1', '10^0', '10^1']
    diff = [1/(60*60*24*365*1000)*tau_diff(b, D) for b in bs]

    f,ax = plt.subplots()
    ax.plot(bs, diff, label="diff", color="r", ls=":")

    for K, label in zip(Ks, labels):
        adv = [1/(60*60*24*365*1000)*tau_adv(b, K, n) for b in bs]
        ax.plot(bs, adv, label=label)

    ax.hlines(y=20, xmin=0, xmax=50, color="k", ls=":")
    ax.legend(title="Kz [m/day]", loc="upper right")
    ax.set_yscale('log')
    ax.set_ylabel("Time [ky]")
    ax.set_xlabel("Thickness [m]")
    
    plt.show()
        

if __name__=="__main__":
    main()