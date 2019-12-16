from numpy import *
import matplotlib.pyplot as plt

def main():
    ef = -2.2927
    enk = 13.605*loadtxt('./enk.dat') - ef
    nbnd, nk = shape(enk)
    wnk = loadtxt('wnk.dat')

    kpath = arange(nk)/nk
    ymin = -30
    ymax = 30
    plt.figure(figsize=(4.5,6))
    plt.subplots_adjust(left=.2)
    wmax = wnk.max()
    for ib, ek in enumerate(enk):
        c = zeros((nk,4))
        c[:,2] = 1
        c[:,3] = wnk[ib,:]/wmax
        plt.scatter(kpath, ek, c=c, edgecolor='none')
        plt.plot(kpath, ek, 'k-')
    plt.ylim(ymin,ymax)
    plt.ylabel('Energy (eV)')
    plt.savefig('unfold.png', dpi=300)
    plt.show()
main()
