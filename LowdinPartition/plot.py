from pylab import *

rcParams['xtick.labelsize'] = 18
rcParams['ytick.labelsize'] = 18

basepath = "data/10ethylene/distance4.5/"
fullfile = basepath+"Veff_lenght_EDAgap0.1_Vx5"
NoExchangefile = basepath+"Veff_lenght_EDAgap0.1_Vx5_NoExchange"
No4CenterExchangefile = basepath+"Veff_lenght_EDAgap0.1_Vx5_No4CenterExchange"

fulldata = genfromtxt(fullfile)
NoExchangedata = genfromtxt(NoExchangefile)
No4CenterExchangedata = genfromtxt(No4CenterExchangefile)




plot(fulldata[:,0], fulldata[:,1], lw='2', label="full data", color='black')
plot(No4CenterExchangedata[:,0], No4CenterExchangedata[:,1], lw=2, label="4-center exchange=0", color='blue')
plot(NoExchangedata[:,0], NoExchangedata[:,1], lw='2', label="Exchange=0", color='red')

betafull = -(log(fulldata[-1, 1]) - log(fulldata[0, 1])) / (4.5 * 7)
betaNoExchange = -(log(NoExchangedata[-1, 1]) - log(NoExchangedata[0, 1])) / (4.5 * 7)
betaNo4CenterExchange = -(log(No4CenterExchangedata[-1, 1]) - log(No4CenterExchangedata[0, 1])) / (4.5 * 7)


yscale('log')
legend(loc=1)
xlabel("Number of Bridge Sites", fontsize=18)
ylabel("V (eV)", fontsize=18)

text(2.5, 1e-12, "%4.2f" % betafull + r"$\AA^{-1}$", color='black', fontsize=18)
text(2.5, 1e-11, "%4.2f" % betaNo4CenterExchange + r"$\AA^{-1}$", color='blue', fontsize=18)
text(2.5, 1e-10, "%4.2f" % betaNoExchange + r"$\AA^{-1}$", color='red', fontsize=18)

savefig(basepath+"Veff_length_EDAgap0.1_Vx5_all.png")
