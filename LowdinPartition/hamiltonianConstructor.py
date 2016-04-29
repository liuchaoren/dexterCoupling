import numpy as np
from pylab import *
from scanning import combineOpt
from numpy.linalg import norm

rcParams['xtick.labelsize'] = 18
rcParams['ytick.labelsize'] = 18
# homoFock = basepath + 'homoFock.npy'
# lumoFock = basepath + 'lumoFock.npy'
hartree = 27.21138602


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


def twoEhandlerNo4CenterExchangePlus(i,j,x,y, twoElectronMap, signlist=((),())):
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
	if i==j or x==y or (i==x and j==y) or (i==y and j==x):
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
	Etunnel = (H[0,0]+H[-1,-1])/2
	# print Hpp
	# print Hqq
	# print Hpq
	# print Hqp
	# print Etunnel * np.eye(dim-2)
	# print (Etunnel * np.eye(dim-2) - Hqq).I * (Etunnel * np.eye(dim-2) - Hqq)
	Heff = Hpp + Hpq * (Etunnel * np.eye(dim-2) - Hqq).I  * Hqp
	return Heff

def CISremoveBonds(H, bonds=[]):
	# bonds are zero-based
	dim = int(np.sqrt(H.shape[0]))
	keepindex= []
	for i in range(dim):
		for x in range(dim):
			if not (i in bonds or x in bonds):
				I = i * dim + x
				keepindex.append(I)

	return H[array(keepindex), :][:, array(keepindex)]

def removeBonds(H, bonds=[]):
	keepindex = []
	dim = H.shape[0]
	for i in range(dim):
		if not (i in bonds):
			keepindex.append(i)
	return H[array(keepindex), :][:, array(keepindex)]


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


def bridgeExcitonSize(H, d):
	indexkeep = []
	for i in range(0,dim):
		for x in range(0, dim):
			index = i * dim + x
			if i == 0 or i == dim - 1 or x == 0 or x == dim - 1:
				indexkeep.append(index)
			elif abs(i-x) <= d:
				indexkeep.append(index)

	return H[array(indexkeep), :][:,array(indexkeep)]

def CTExcitonHamiltonian(H, dim):
	indexkeep = []
	for i in range(dim):
		for x in range(dim):
			index =  i * dim + x
			if not (i in range(1,dim-1) and x in range(1,dim-1)):
				indexkeep.append(index)

	return H[array(indexkeep), :][:, array(indexkeep)]

def BridgeExcitonHamiltonian(H, dim):
	indexkeep = []
	for i in range(dim):
		for x in range(dim):
			index = i * dim + x
			if not ((i == 0 and x == dim -1) or (i== dim-1 and x==0)):
				indexkeep.append(index)

	return H[array(indexkeep), :][:, array(indexkeep)]
#
# def CISHamiltonian():
# 	H = np.zeros([dim**2,dim**2])
# 	for i in range(dim):
# 		for x in range(dim):
# 			I = i * dim + x
# 			for j in range(dim):
# 				for y in range(dim):
# 					J = j * dim + y
# 					Hixjy=0.0
# 					if i == j:
# 						Fxy = lumoFockM[x,y]
# 						Hixjy = Hixjy + Fxy
# 					if x == y:
# 						Fij = homoFockM[i,j]
# 						Hixjy = Hixjy - Fij
# 					# twoEijxy = twoEhandlerNoExchange(i, j, x, y, twoElectronMap)
# 					# twoEijxy = twoEhandlerNo4CenterExchange(i, j, x, y, twoElectronMap)
# 					# twoEijxy = twoEhandlerNo4CenterExchangePlus(i, j, x, y, twoElectronMap)
# 					twoEijxy = twoEhandler(i, j, x, y, twoElectronMap)
# 					Hixjy = Hixjy - twoEijxy
# 					H[I, J] = Hixjy
# 	return H

def EqualDA(H):
    ED = H[0,0]
    EA = H[-1,-1]
    aveDA = (ED+EA)/2
    H[0,0] = aveDA
    H[-1,-1] = aveDA
    return H

def lowerGap(H, dim, deltaE):
    # deltaE is the energy shift of bridge HOMO (up) and bridge LUMO (down)
    for i in range(dim):
        for x in range(dim):
            Iindex = i * dim + x
            for j in range(dim):
                for y in range(dim):
                    Jindex = j * dim + y
                    if (i==j and x==y):
                        if (not x in [0, dim-1]):
                            H[Iindex, Jindex] = H[Iindex, Jindex] - deltaE
                        if (not i in [0, dim-1]):
                            H[Iindex, Jindex] = H[Iindex, Jindex] - deltaE
    return H



#
# def HoleHamiltonian():
# 	return homoFockM
#
# def EHamiltonian():
# 	return lumoFockM


def coulomb(i, x, centers):
	r = norm(centers[i] - centers[x]) / 1e10
	e = 1.6e-19
	return 8.9875517873681764e9 * e**2 / r / 1.6e-19 / hartree

def CISThreeWayJ(Hcouplings, Lcouplings, centers):
	dim = 21
	H = np.zeros([dim**2,dim**2])
	order = [18] + range(18) + [20, 19]
	sequence = ['T', 'T', 'T', 'T', 'T', 'T'] + ['C', 'C', 'C', 'A', 'A', 'A'] + ['A', 'A', 'A', 'G', 'G', 'G'] + ['porphyrin', 'porphyrin', 'porphyrin']
	Ahomo = -7.8226 / hartree
	Alumo = 0.322 / hartree
	Chomo = -8.7128 / hartree
	Clumo = -0.3126 / hartree
	Ghomo = -7.5192 / hartree
	Glumo = 0.1484 / hartree
	Thomo = -9.082 / hartree
	Tlumo = -0.5052 /hartree
	Phomo = (-8.28415 - 3) / hartree
	Plumo = (-0.087 - 3) / hartree

	for t in range(21):
		if sequence[t] == 'T':
			Hcouplings[t,t] = Thomo
			Lcouplings[t,t] = Tlumo
		elif sequence[t] == 'A':
			Hcouplings[t,t] = Ahomo
			Lcouplings[t,t] = Alumo
		elif sequence[t] == 'G':
			Hcouplings[t,t] = Ghomo
			Lcouplings[t,t] = Glumo
		elif sequence[t] == 'C':
			Hcouplings[t,t] = Chomo
			Lcouplings[t,t] = Clumo
		elif sequence[t] == 'porphyrin':
			Hcouplings[t,t] = Phomo
			Lcouplings[t,t] = Plumo
		else:
			raise Exception

	for m in range(dim):
		i = order[m]
		for n in range(dim):
			x = order[n]
			I = m * dim + n
			for t in range(dim):
				j = order[t]
				for s in range(dim):
					y = order[s]
					J = t * dim + s
					Hixjy=0.0
					if i==j:
						Hixjy = Hixjy + Lcouplings[x,y]
					if x==y:
						Hixjy = Hixjy - Hcouplings[i,j]
					if i==j and x==y:
						Hixjy = Hixjy - coulomb(i,x, centers)
					H[I,J] = Hixjy
	return H




if __name__ == "__main__":
	lowerhomoslumos = 0.00 # move homos of bridge sites up by lowerhomoslumos and move lumos of bridge sites down by lowerhomoslumos
	centers = load("../g09Extract/data/3wj/zindo/centers.npy")
	homoFock = load("../g09Extract/data/3wj/zindo/homoFock.npy" )
	lumoFock = load("../g09Extract/data/3wj/zindo/lumoFock.npy" )

	CISH = CISThreeWayJ(homoFock, lumoFock, centers)
	save('CISH.npy', CISH)

    # output = open("results/alkanes/coupling_noExchange_lower%s" % lowerhomoslumos, "w")
    # output = open("results/polynorbornyl/parallel/coupling_allExchange_lower%s" % lowerhomoslumos, 'a')
    # output.write("# couplings of of polynorbornyl with 2 repeats and 3 repeats. Exchanges are all included\n\n# repeats num 	total     BE      CT\n")
    # iteratornum = 2000
    # scanarray1=arange(-1e-8, 1e-8, 1e-11)
    # scanarray2=arange(-1e-8, 1e-8, 1e-11)
    # scanarray3=arange(-1e-9, 1e-9, 1e-11)

    # for i in [4, 5, 6, 7, 8, 9, 10, 11, 12, 13]:
    #     if i in [7, 13]:
    #         twoE = "../G092E/data/alkyl/b%s/linear/symm/twoEtwoElectron" % i
    #         homoFockFile = "../g09Extract/data/alkyl/b%s/linear/symm/homoFock.npy" % i
    #         lumoFockFile = "../g09Extract/data/alkyl/b%s/linear/symm/lumoFock.npy" % i
    #     else:
    #         twoE = "../G092E/data/alkyl/b%s/symm/twoEtwoElectron" % i
    #         homoFockFile = "../g09Extract/data/alkyl/b%s/symm/homoFock.npy" % i
    #         lumoFockFile = "../g09Extract/data/alkyl/b%s/symm/lumoFock.npy" % i
    #
    # for i in [3]:
    #     twoE = "../G092E/data/alkyl/b13/folded/twoEtwoElectron%s" % i
    #     homoFockFile = "../g09Extract/data/alkyl/b13/folded/homoFock_%s.npy" % i
    #     lumoFockFile = "../g09Extract/data/alkyl/b13/folded/lumoFock_%s.npy" % i
    #
    # for i in [3]:
    #     twoE = "../G092E/data/polynorbornyl/parallel/repeats_%s/twoEtwoElectron" % i
    #     homoFockFile = "../g09Extract/data/polynorbornyl/parallel/repeats_%s/homoFock.npy" % i
    #     lumoFockFile = "../g09Extract/data/polynorbornyl/parallel/repeats_%s/lumoFock.npy" % i

        # homoFockM = np.load(homoFockFile)
        # lumoFockM = np.load(lumoFockFile)
        # twoElectronMap = twoEmap(twoE)
        # EDAgap = 0.1
        # dim =homoFockM.shape[0]

        # H = CISHamiltonian()
        # HoleH = HoleHamiltonian()
        # ElectronH = EHamiltonian()
        # H = lowerGap(H, dim, lowerhomoslumos)
        # print lowdinPartition(H)[0,1]
        # H = EqualDA(H)
        # HBE = BridgeExcitonHamiltonian(H, dim)
        # HCT = CTExcitonHamiltonian(H,dim)
        # save("H%s_0.0.npy" % i, H)
        # save("HBE%s_0.0.npy" % i, HBE)
        # save("HCT%s_0.0.npy" % i, HCT)
        # (evalH,evecH) = eig(H)
        # assert(abs(evecH[0,0]) > 0.6 and abs(evecH[-1,0]) > 0.6 and abs(evecH[0,1])> 0.6 and abs(evecH[-1,1]) > 0.6)
        # (evalHBE,evecHBE) = eig(HBE)
        # assert(abs(evecHBE[0,0]) > 0.6 and abs(evecHBE[-1,0]) > 0.6 and abs(evecHBE[0,1])> 0.6 and abs(evecHBE[-1,1]) > 0.6)
        # (evalHCT,evecHCT) = eig(HCT)
        # assert(abs(evecHCT[0,0]) > 0.6 and abs(evecHCT[-1,0]) > 0.6 and abs(evecHCT[0,1])> 0.6 and abs(evecHCT[-1,1]) > 0.6)

        # totalcoup = combineOpt(H, iteratornum, scanarray1) * hartree
        # print "total calculation for %s is done" % i
        # print "total CIS coupling is", totalcoup
        # BEcoup = combineOpt(HBE, iteratornum, scanarray2) * hartree
        # print "BE calculation for %s is done" % i
        # print "total BE CIS coupling is", BEcoup
        # CTcoup = combineOpt(HCT, iteratornum, scanarray3) * hartree
        # print "CT calculation for %s is done" % i
        # print "total CT CIS coupling is", CTcoup

        # totalcoup = abs((evalH[0]-evalH[1])/2 * hartree)
        # BEcoup = abs((evalHBE[0]-evalHBE[1])/2 * hartree)
        # CTcoup = abs((evalHCT[0]-evalHCT[1])/2 * hartree)
        # totalcoup = lowdinPartition(H)[0,1] * hartree
        # totalHoleCoup = lowdinPartition(HoleH)[0,1] * hartree
        # totalElectronCoup = lowdinPartition(ElectronH)[0,1] * hartree
        # BEcoup = lowdinPartition(HBE)[0,1] * hartree
        # CTcoup = lowdinPartition(HCT)[0,1] * hartree
        # H = EqualDA(H)
        # HBE = EqualDA(HBE)
        # HCT = EqualDA(HCT)
        # HoleH = EqualDA(HoleH)
        # ElectronH = EqualDA(ElectronH)
        # (eval, evec) = eig(H)
        # print evec[:,0]
        # print evec[:,1]
        # print "total CIS coupling is", totalcoup
        # print "total BE CIS coupling is", BEcoup
        # print "total CT CIS coupling is", CTcoup
        # print "total hole coupling is", totalHoleCoup
        # print "total electron coupling is", totalElectronCoup
        #
        # output.write("%d\t%e\t%e\t%e\n" % (i, totalcoup, BEcoup, CTcoup))
        # HCT = CTExcitonHamiltonian(H)
        # save("data/alkyl/b7/symm/H_allExchanges.npy", H)












