'''
The HOMO/LUMO orbitals I used in constructing the CIS excited states is non-orthogonal. Therefore, the formula for Hamiltonian elements
in Equatioin 3 in our PNAS paper does not hold. This script is trying to calculate the Hamiltonian elements from Slater determinant
which take the non-orthogonality into account.
'''

from numpy import *
from itertools import permutations

class orbital:
    def __init__(self, homo_lumo_i, orbital_i, spin):
        self.orbital_i = orbital_i # count from 0
        self.spin = spin # 0 - spin up, 1 - spin down
        self.homo_lumo_i = homo_lumo_i # 0 - homo, 1 - lumo

    def printOribital(self):
        print "I am a %s with index of %s and spin of %s" % ("HOMO" if self.homo_lumo_i==0 else "LUMO", self.orbital_i, "up" if self.spin == 0 else "down")


def arePermsEqualParity(perm0, perm1):
    """Check if 2 permutations are of equal parity.

    Assume that both permutation lists are of equal length
    and have the same elements. No need to check for these
    conditions.
    """
    perm1 = perm1[:] ## copy this list so we don't mutate the original

    transCount = 0
    for loc in range(len(perm0) - 1):                         # Do (len - 1) transpositions
        p0 = perm0[loc]
        p1 = perm1[loc]
        if p0 != p1:
            sloc = perm1[loc:].index(p0)+loc          # Find position in perm1
            perm1[loc], perm1[sloc] = p0, p1          # Swap in perm1
            transCount += 1

    # Even number of transpositions means equal parity
    if (transCount % 2) == 0:
        return True
    else:
        return False


def Slater(total_site_num, hole_i, hole_spin, electron_i, electron_spin, overlap):
    # print "I am generating one Slater"
    orbitalList = []
    for i in range(0,total_site_num):
        orbitalList.append(orbital(0, i, 0))
        orbitalList.append(orbital(0, i, 1))

    if (hole_spin==0):
        orbitalList[hole_i*2] = orbital(1, electron_i, electron_spin)
    else:
        orbitalList[hole_i*2+1] = orbital(1, electron_i, electron_spin)

    basePerm = range(0,total_site_num*2)
    slaterPerms = []
    slaterSigns = []

    for eachPerm in permutations(basePerm):
        eachFactor = []
        for eachO in eachPerm:
            eachFactor.append(orbitalList[eachO])

        slaterPerms.append(eachFactor)
        if arePermsEqualParity(eachPerm, basePerm):
            slaterSigns.append(1)
        else:
            slaterSigns.append(-1)


    sum = 0.
    for i in range(len(slaterSigns)):
        oneSign = slaterSigns[i]
        oneFactor = slaterPerms[i]
        for j in range(len(slaterSigns)):
            twoSign = slaterSigns[j]
            twoFactor = slaterPerms[j]
            product=1.
            for m in range(len(oneFactor)):
                if oneFactor[m].spin == twoFactor[m].spin:
                    product *= overlap[oneFactor[m].homo_lumo_i * total_site_num + oneFactor[m].orbital_i, twoFactor[m].homo_lumo_i * total_site_num + twoFactor[m].orbital_i]
                else:
                    product = 0.
                    break
            # print product

            sum += (product * oneSign * twoSign)

    return (slaterPerms, slaterSigns, 1/sqrt(sum))

class oneElectron:
    def __init__(self, electron_i):
        self.electron_i = electron_i

class twoElectron:
    def __init__(self, electron_one, electron_two):
        self.electron_one = electron_one
        self.electron_two = electron_two

class Hamiltonian:
    def __init__(self, electron_num):
        oneElectronPart = []
        twoElectronPart = []
        for i in range(electron_num):
            oneElectronPart.append(oneElectron(i))
            for j in range(electron_num):
                if i > j:
                    twoElectronPart.append(twoElectron(i,j))

        self.oneElectronPart = oneElectronPart
        self.twoElectronPart = twoElectronPart


def oneElectronElement(total_site_num, firstFactor, firstSign, secondFactor, secondSign, myOneElectron, overlap, orbitalOneElectron):
    # print "I am calculating one electron terms"
    product=1.
    electron_index = myOneElectron.electron_i
    for i in range(len(firstFactor)):
        if firstFactor[i].spin == secondFactor[i].spin:
            if i == electron_index:
                product *= orbitalOneElectron[firstFactor[i].homo_lumo_i * total_site_num + firstFactor[i].orbital_i, secondFactor[i].homo_lumo_i * total_site_num + secondFactor[i].orbital_i]
            else:
                product *= overlap[firstFactor[i].homo_lumo_i * total_site_num + firstFactor[i].orbital_i, secondFactor[i].homo_lumo_i * total_site_num + secondFactor[i].orbital_i]
        else:
            return 0.

    return product * firstSign * secondSign


def twoElectronElement(total_site_num, firstFactor, firstSign, secondFactor, secondSign, myTwoElectron, overlap, orbitalTwoElectron):
    # print "I am calculating two electron terms"
    product = 1.
    electron_one = myTwoElectron.electron_one
    electron_two = myTwoElectron.electron_two

    for i in range(len(firstFactor)):
        if firstFactor[i].spin != secondFactor[i].spin:
            return 0.
        else:
            if (i != electron_one) and (i != electron_two):
                product *= overlap[firstFactor[i].homo_lumo_i * total_site_num + firstFactor[i].orbital_i, secondFactor[i].homo_lumo_i * total_site_num + secondFactor[i].orbital_i]

    I = min(firstFactor[electron_one].homo_lumo_i * total_site_num + firstFactor[electron_one].orbital_i, secondFactor[electron_one].homo_lumo_i * total_site_num + secondFactor[electron_one].orbital_i)
    J = max(firstFactor[electron_one].homo_lumo_i * total_site_num + firstFactor[electron_one].orbital_i, secondFactor[electron_one].homo_lumo_i * total_site_num + secondFactor[electron_one].orbital_i)
    X = min(firstFactor[electron_two].homo_lumo_i * total_site_num + firstFactor[electron_two].orbital_i, secondFactor[electron_two].homo_lumo_i * total_site_num + secondFactor[electron_two].orbital_i)
    Y = max(firstFactor[electron_two].homo_lumo_i * total_site_num + firstFactor[electron_two].orbital_i, secondFactor[electron_two].homo_lumo_i * total_site_num + secondFactor[electron_two].orbital_i)
    if (I * total_site_num * 2 + J <= X * total_site_num * 2 + Y):
        product *= orbitalTwoElectron[(I, J, X, Y)]
    else:
        product *= orbitalTwoElectron[(X, Y, I, J)]
    return product * firstSign * secondSign


def HamiltonianElement(total_site_num, i, x, j, y, overlap, orbitalOneElectron, orbitalTwoElectron):
    (firstCIPerm, firstCISign, firstNorm) = Slater(total_site_num, i, 1, x, 0, overlap)
    (secondCIPerm, secondCISign, secondNorm) = Slater(total_site_num, j, 1, y, 0, overlap)
    myHamiltonian = Hamiltonian(total_site_num*2)
    sum=0.
    for i in range(len(firstCIPerm)):
        firstFactor = firstCIPerm[i]
        firstSign = firstCISign[i]
        for j in range(len(secondCIPerm)):
            secondFactor = secondCIPerm[j]
            secondSign = secondCISign[j]
            for eachOneElectron in myHamiltonian.oneElectronPart:
                sum += oneElectronElement(total_site_num, firstFactor, firstSign, secondFactor, secondSign, eachOneElectron, overlap, orbitalOneElectron)
            for eachTwoElectron in myHamiltonian.twoElectronPart:
                sum += twoElectronElement(total_site_num, firstFactor, firstSign, secondFactor, secondSign, eachTwoElectron, overlap, orbitalTwoElectron)

    return sum * firstNorm * secondNorm

def twoEmap(filename):
	twoEitems = genfromtxt(filename, dtype=[int, int, int, int, float])
	twoEmap = {}
	for eachitem in twoEitems:
		(i, j, x, y) = (eachitem[0], eachitem[1], eachitem[2], eachitem[3])
		value = eachitem[4]
		twoEmap[(i,j,x,y)] = value
	return twoEmap




if __name__ == "__main__":
    overlap = load("./AT/AT-pair/overlap_HOMOLUMO.npy")
    orbitalOneElectron = load("./AT/AT-pair/h_HOMOLUMO.npy")
    orbitalTwoElectron = twoEmap("./AT/AT-pair/twoElectron_HOMOLUMO")
    print HamiltonianElement(2, 0, 0, 0, 1, overlap, orbitalOneElectron, orbitalTwoElectron)




