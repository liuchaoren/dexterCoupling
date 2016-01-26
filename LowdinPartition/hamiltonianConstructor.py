import numpy as np
from pylab import *

rcParams['xtick.labelsize'] = 18
rcParams['ytick.labelsize'] = 18
# homoFock = basepath + 'homoFock.npy'
# lumoFock = basepath + 'lumoFock.npy'


def twoEmap(filename):
	twoEitems = np.genfromtxt(filename, dtype=[int, int, int, int, float])
	twoEmap = {}
	for eachitem in twoEitems:
		(i, j, x, y) = (eachitem[0], eachitem[1], eachitem[2], eachitem[3])
		value = eachitem[4]
		twoEmap[(i,j,x,y)] = value
	return twoEmap


def twoEhandler(i,j,x,y, twoElectronMap, signlist=((),())):
	hsign = signlist[0]
	lsign = signlist[1]
	counterchange = 0
	for eachhsign in hsign:
		for eachcal in (i,j):
			if eachhsign == eachcal:
				counterchange = counterchange + 1
	for eachlsign in lsign:
		for eachcal in (i,j):
			if eachlsign == eachcal:
				counterchange = counterchange + 1

	# return (ij|xy)
	newi = max(i,j)
	newj = min(i,j)
	newx = max(x,y)
	newy = min(x,y)
	return (-1)**counterchange * twoElectronMap[(newi, newj, newx, newy)]

def twoEhandlerNoExchange(i,j,x,y, twoElectronMap, signlist=((),())):
	hsign = signlist[0]
	lsign = signlist[1]
	counterchange = 0
	for eachhsign in hsign:
		for eachcal in (i,j):
			if eachhsign == eachcal:
				counterchange = counterchange + 1
	for eachlsign in lsign:
		for eachcal in (i,j):
			if eachlsign == eachcal:
				counterchange = counterchange + 1
	# return (ij|xy)
	if i==j and x==y:
		newi = max(i,j)
		newj = min(i,j)
		newx = max(x,y)
		newy = min(x,y)
		return (-1)**counterchange * twoElectronMap[(newi, newj, newx, newy)]
	else:
		return 0

def twoEhandlerNo4CenterExchange(i,j,x,y, twoElectronMap, signlist=((),())):
	hsign = signlist[0]
	lsign = signlist[1]
	counterchange = 0
	for eachhsign in hsign:
		for eachcal in (i,j):
			if eachhsign == eachcal:
				counterchange = counterchange + 1
	for eachlsign in lsign:
		for eachcal in (i,j):
			if eachlsign == eachcal:
				counterchange = counterchange + 1
	# return (ij|xy)
	if i==j or x==y:
		newi = max(i,j)
		newj = min(i,j)
		newx = max(x,y)
		newy = min(x,y)
		return (-1)**counterchange * twoElectronMap[(newi, newj, newx, newy)]
	else:
		return 0


def fockElement(Ehomo, Elumo, EDAhomo, EDAlumo, Vhomo, Vlumo, i, j, HLindicator, dim, signlist=()):
	if HLindicator == "homo":
		if i == j:
			if i == 0 or i == dim-1: 
				return EDAhomo
			else:
				return Ehomo
		elif abs(i-j) == 1:
			if ((i in signlist) and (not j in signlist)) or ((i not in signlist) and (j in signlist)):
				return -1*Vhomo
			else:
				return Vhomo
		else:
			return 0
	elif HLindicator == "lumo":
		if i == j:
			if i == 0 or i == dim-1: 
				return EDAlumo
			else:
				return Elumo
		elif abs(i-j) == 1:
			if ((i in signlist) and (not j in signlist)) or ((i not in signlist) and (j in signlist)):
				return -1*Vlumo
			else: 
				return Vlumo
		else:
			return 0

def lowdinPartition(H):
	dim = H.shape[0]
	Hpp = np.matrix([[H[0, 0], H[0, -1]], [H[-1, 0], H[-1,-1]]])
	# print Hpp

	Hpq = np.matrix(H[np.array([0, -1]), 1:-1])
	Hqp = np.matrix(H[1:-1, np.array([0, -1])])
	Hqq = np.matrix(H[1:-1, 1:-1])
	Etunnel = H[0,0]
	# print Hpp
	# print Hqq
	# print Hpq
	# print Hqp
	# print Etunnel * np.eye(dim-2)
	# print (Etunnel * np.eye(dim-2) - Hqq).I * (Etunnel * np.eye(dim-2) - Hqq)
	Heff = Hpp + Hpq * (Etunnel * np.eye(dim-2) - Hqq).I  * Hqp
	return Heff

def plotVeff(bridgeLen, Veffs, basepath, figurename):
	plot(bridgeLen, Veffs)
	yscale('log')
	xlabel("Number of Bridge Sites", fontsize=18)
	ylabel("V (eV)", fontsize=18)
	savefig(basepath+figurename)

def saveVeff(bridgeLen, Veffs, basepath, filename):
	fileobject = open(basepath+filename, 'w')
	for (bridgenum, Veff) in zip(bridgeLen, Veffs):
		fileobject.write("%d\t%e\n" % (bridgenum, Veff))



if __name__ == "__main__":
	basepath = 'data/10ethylene/distance4.5/'
	twoE = basepath + 'twoEtwoElectron'
	hartree = 27.21138602
	# homoFockM = np.load(homoFock)
	# lumoFockM = np.load(lumoFock)
	twoElectronMap = twoEmap(twoE)
	EDAgap = 0.1

	# Vhomo = -0.06199605
	# Vlumo = -0.02805976

	Vhomo = -0.00525939
	Vlumo = -0.00420136 

	Ehomo = -0.26449
	Elumo = 0.01314
	Egap = Elumo - Ehomo
	EDAhomo = Ehomo + Egap/2. - EDAgap/2
	EDAlumo = Ehomo + Egap/2. + EDAgap/2
	bridgeLen = []
	Veffs=[]

	for dim in range(3, 11):
		H = np.zeros([dim**2, dim**2])
		for i in range(dim):
			for x in range(dim):
				I = i * dim + x 
				for j in range(dim):
					for y in range(dim):
						J = j * dim + y
						Hixjy = 0.
						if i == j:
							Fxy = fockElement(Ehomo, Elumo, EDAhomo, EDAlumo, Vhomo, Vlumo, x, y, "lumo", dim, signlist=())
							Hixjy = Hixjy + Fxy
						if x == y: 
							Fij = fockElement(Ehomo, Elumo, EDAhomo, EDAlumo, Vhomo, Vlumo, i, j, "homo", dim, signlist=())
							Hixjy = Hixjy - Fij
						# twoEijxy = twoEhandlerNoExchange(i,j,x,y, twoElectronMap, signlist=((),()))
						twoEijxy = 0.0
						Hixjy = Hixjy - twoEijxy
						# H[I, J] = Hixjy
						if (i>0 and i<dim-1 and x>0 and x<dim-1) or (j>0 and j<dim-1 and y>0 and y<dim-1):
							H[I, J] = 0
						else:
							H[I, J] = Hixjy

		# np.save(basepath+"dexterH.npy", H)
		# np.savetxt("H3", H, fmt="%10.5e")
		Veffs.append(abs(lowdinPartition(H)[0,1]) * hartree)
		bridgeLen.append(dim-2)
		plotVeff(bridgeLen, Veffs, basepath, "Veff_length_EDAgap0.1_NoExchange_NoCoulomb.png")
		saveVeff(bridgeLen, Veffs, basepath, "Veff_lenght_EDAgap0.1_NoExchange_NoCoulomb")


	# for (bridgenum, Veff) in zip(bridgeLen, Veffs):






