'''
calculate two E integrals of four contracted gaussian functions
Author - Chaoren Liu
Date - Jan 25, 2016
'''
import numpy as np
# from scipy.sparse import csr_matrix

# from PyQuante.pyints import contr_coulomb as pycc
# from PyQuante.rys import contr_coulomb as ryscc
# from PyQuante.hgp import contr_coulomb as hgpcc
# from PyQuante.cints import contr_coulomb as ccc
from PyQuante.crys import contr_coulomb as cryscc
# from PyQuante.chgp import contr_coulomb as chgpcc
from PyQuante.Ints import getbasis
from PyQuante.Molecule import Molecule
from datetime import datetime

from PyQuante2E.converters import *

class sparseMatrix:
    def __init__(self, rows, cols, values, shape):
        self.rows = rows
        self.cols = cols
        self.values = values
        self.shape = shape


def coulomb(a, b, c, d, coul_func):
    # a, b, c, d are four basis in a basis set (for example 1s of a carbon (coordinate is known) in 6-31g)
    # calculate the two electron integral between 4 contracted gaussian functions
    Jij = coul_func(a.pexps,a.pcoefs,a.pnorms,a.origin,a.powers,
                    b.pexps,b.pcoefs,b.pnorms,b.origin,b.powers,
                    c.pexps,c.pcoefs,c.pnorms,c.origin,c.powers,
                    d.pexps,d.pcoefs,d.pnorms,d.origin,d.powers)
    return a.norm*b.norm*c.norm*d.norm*Jij


def crossProductSparse(a, b):
    # a, b are two arrays, return a sparse matrix m where m[i, j] = a[i] * b[j]
    basissize = len(a)
    tol = 1e-10
    rows = []
    cols = []
    values = []
    for i in range(basissize):
        for j in range(basissize):
            if abs(a[i]) > tol and abs(b[j]) > tol:
                rows.append(i)
                cols.append(j)
                values.append(a[i]*b[j])

    return sparseMatrix(rows, cols, values, (basissize, basissize))


def moCoulomb(I, J, X, Y, mol, basisset='6-31g', coul_func=cryscc):
    # I, J, X, Y are four MO's. I,J are HOMOs and X,Y are LUMOs
    # mol is the Molecule object, coul_func is the method choice
    bfs = getbasis(mol, basisset)
    IJmatrix = crossProductSparse(I, J)
    XYmatrix = crossProductSparse(X, Y)
    result = 0.0
    # print len(IJmatrix.cols)
    # print len(XYmatrix.cols)
    totalLoops = len(IJmatrix.cols) * len(XYmatrix.cols)

    c = 0
    for t in range(len(IJmatrix.cols)):
        for s in range(len(XYmatrix.cols)):
            i = IJmatrix.rows[t]
            j = IJmatrix.cols[t]
            k = XYmatrix.rows[s]
            l = XYmatrix.cols[s]
            twoE = coulomb(bfs[i], bfs[j], bfs[k], bfs[l], coul_func)
            result = result + twoE * IJmatrix.values[t] * XYmatrix.values[s]
            if c == 0:
                currentTime = datetime.now()
            elif c%100000 == 0:
                currentTimePre = currentTime
                currentTime = datetime.now()
                deltaTime = currentTime - currentTimePre
                print "%d loops in %d loops (%e) are done!" % (c, totalLoops, c*1.0/totalLoops)
                print "remaining calculation time is %s" % str(((totalLoops-c)/100000) * deltaTime)
            c += 1

    return result

def simplify(evecs,thr, overlap):
    returnEvecs = np.zeros(evecs.shape)
    syslen = evecs.shape[1]
    for i in range(syslen):
        oneEvec = evecs[:,range(i,i+1)]
        for j in range(oneEvec.shape[0]):
            if abs(oneEvec[j,0]) <  thr:
                oneEvec[j,0] = 0.
        renormal = sqrt((np.matrix(oneEvec.transpose()) * overlap * oneEvec)[0,0])
        oneEvec = oneEvec / renormal
        returnEvecs[:,range(i,i+1)] = oneEvec

    print returnEvecs.transpose() * matrix(overlap) * returnEvecs
    return returnEvecs


if __name__ == "__main__":
    homoevecs = np.load('data/AT/homos.npy')
    lumoevecs = np.load('data/AT/lumos.npy')
    overlap = np.load('../g09Extract/data/AT/overlap.npy')
    thr = 1e-2
    homoevecs = simplify(homoevecs, thr, overlap)
    lumoevecs = simplify(lumoevecs, thr, overlap)
    # b7
    # HOMOs = [38, 12, 18, 14, 16, 15, 17, 13, 37] # homos in order of D, B1, B2, ... A. zero-based
    # LUMOs = [40, 47, 41, 45, 43, 44, 42, 46, 39]

    # nlmos = np.load("../g09Extract/data/alkyl/nlmos.npy")
    gjffilename = 'data/AT/AT.com'
    outputname = 'data/AT/2EIntegrals'
    output = open(outputname, 'w')
    molecule = gjf2Molecule(gjffilename)
    print molecule
    # I = homoevecs[:, 1]
    # J = homoevecs[:, 1]
    # X = lumoevecs[:, 0]
    # Y = lumoevecs[:, 0]

    # print moCoulomb(I, J, X, Y, molecule)
    # syslen = len(HOMOs)
    syslen = homoevecs.shape[1]
    assert(syslen==2)

    for i in range(syslen):
        # I = nlmos[:, HOMOs[i]]
        I = homoevecs[:, i]
        for j in range(i, syslen):
            J = homoevecs[:, j]
            for k in range(syslen):
                # X = nlmos[:, LUMOs[k]]
                X = lumoevecs[:, k]
                for l in range(k, syslen):
                    Y = lumoevecs[:, l]
                    twoEvalue =  moCoulomb(I, J, X, Y, molecule)
                    # output.write("%d\t%d\t%d\t%d\t%e\n" % (min(HOMOs[i], HOMOs[j]), max(HOMOs[i], HOMOs[j]), min(LUMOs[k], LUMOs[l]), max(LUMOs[k], LUMOs[l]), twoEvalue))
                    output.write("%d\t%d\t%d\t%d\t%e\n" % (i, j, k, l, twoEvalue))
                    print "done with one 2e calculation!"




