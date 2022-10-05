# these are needed so that the notebook can also be run on Python 2
from __future__ import division, print_function
from pylab import *
from semi_interface import *

def main():
    # sc1 = SemiCoast(k=3, H=16, c=250000, grad=0.0003, 
    #             rhof=1000, rhos=1025, Ls=50000, 
    #             ztop=-130, sealevel=0)
    # print('toe of interface at:', sc1.toe())
    # print('tip of interface at:', sc1.tip())
    # sc1.plot(xmin=-10000, xmax=50000)

    sc1 = SemiCoast(k=3, H=16, c=250000, grad=0.003, 
                rhof=1000, rhos=1025, Ls=1000, 
                ztop=-130, sealevel=-130)
    print('toe of interface at:', sc1.toe())
    print('tip of interface at:', sc1.tip())
    sc1.plot(xmin=-60000, xmax=1000)


if __name__=="__main__":
    main()  