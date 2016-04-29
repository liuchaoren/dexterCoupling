'''
plot the results
'''

from pylab import *
rcParams['axes.labelsize'] = 24
rcParams['axes.titlesize'] = 18
rcParams['legend.fontsize'] = 20
rcParams['lines.linewidth'] = 3
rcParams['xtick.labelsize'] = 24
rcParams['ytick.labelsize'] = 24
rcParams['legend.frameon'] = False

def plotBECT(filename):
    coups = genfromtxt(filename)
    bridgeNum = coups[:,0]

    allcoup = coups[:,1]
    BEcoup = coups[:,2]
    CTcoup = coups[:,3]
    fig = figure()
    plot(bridgeNum, allcoup, marker='o', label=r"$V_{tr}$", color='black')
    plot(bridgeNum, BEcoup, marker='o', label=r"$V_{tr}^{(be)}$", color='blue')
    plot(bridgeNum, CTcoup, marker='o', label=r"$V_{tr}^{(dae)}$", color='green')
    xlabel(r"number of bridge $\sigma$ bonds")
    ylabel("coupling/eV")
    legend()
    yscale('log')
    tight_layout()
    savefig("%s.png" % filename)
    close(fig)

def plotBECTContr(filename):
    fig = figure()
    coups = genfromtxt(filename)
    bridgeNum = coups[:,0]

    allcoup = coups[:,1]
    BEcoup = coups[:,2]
    CTcoup = coups[:,3]

    BEcontr = BEcoup/(BEcoup+CTcoup)
    CTcontr = CTcoup/(BEcoup+CTcoup)
    plot(bridgeNum, BEcontr, marker='o', label='BE contribution',lw=4)
    plot(bridgeNum, CTcontr, marker='o', label='DAE contribution',lw=4)
    xlabel(r"number of bridge $\sigma$ bonds",fontsize=34)
    ylabel(r"fraction of $V_{tr}$",fontsize=40)
    legend(loc=4,fontsize=30)
    xticks([4, 6, 8, 10, 12], fontsize=32)
    yticks([0.0, 0.5, 1.0], fontsize=32)
    tight_layout()
    savefig('%s_contr.png' % filename)
    close(fig)

def plotExchangeImportance(filename):
    fig=figure()
    coups = genfromtxt(filename)
    bridgeNum = coups[:,0]
    coups_allExchange= coups[:,1]
    coups_noExchange = coups[:, 2]
    error = abs(coups_allExchange - coups_noExchange) / coups_allExchange
    print error
    print bridgeNum
    plot(bridgeNum, error, marker='o')
    xlabel("number of bridge sites")
    ylabel("errors without exchanges")
    tight_layout()
    savefig("%s.png"% filename)

if __name__=='__main__':
    filename = "alkanes/coupling_allExchange"
    # plotBECT(filename)
    plotBECTContr(filename)
    # plotExchangeImportance(filename)
    filename = "alkanes/coupling_allExchange_lower0.19"
    # plotBECT(filename)
    plotBECTContr(filename)
