"""
generate fockmatrix and 2E integrals between orthogonalized HOMO, LUMO states
"""
from numpy import *

def fockMatrix_orth(fock_AO, homos, lumos, orth_homolumo, syslen):
    homoslumos = hstack((homos, lumos))
    fock_local = matrix(homoslumos).transpose() * matrix(fock_AO) * matrix(homoslumos)
    fock_orth = matrix(orth_homolumo).transpose() * matrix(fock_local) * matrix(orth_homolumo)
    homoFock = fock_orth[:syslen, :syslen]
    lumoFock = fock_orth[syslen:, syslen:]
    return [homoFock, lumoFock]

def twoEmap(filename):
	twoEitems = genfromtxt(filename, dtype=[int, int, int, int, float])
	twoEmap = {}
	for eachitem in twoEitems:
		(i, j, x, y) = (eachitem[0], eachitem[1], eachitem[2], eachitem[3])
		value = eachitem[4]
		twoEmap[(i,j,x,y)] = value
	return twoEmap


def twoE_local_2d(twoE_local_filename, orth_homolumo, syslen):
    twoE_local = twoEmap(twoE_local_filename)

    twoE = zeros(((syslen * 2) ** 2, (syslen * 2) ** 2))
    for i in range(syslen*2):
        for j in range(syslen*2):
            i1 = min(i, j)
            j1 = max(i, j)
            I = i1 * syslen * 2 + j1
            for x in range(2*syslen):
                for y in range(2*syslen):
                    x1 = min(x, y)
                    y1 = max(x, y)
                    J = x1 * syslen * 2 + y1
                    if (I <= J):
                        twoE[I,J] = twoE_local[(i1, j1, x1, y1)]
                    else:
                        twoE[I,J] = twoE_local[(x1, y1, i1, j1)]
    return twoE


def crossProduct(a, b):
    # a, b are two arrays, return a sparse matrix m where m[i, j] = a[i] * b[j]
    # print len(a)
    basissize = len(a)
    result = zeros(basissize**2)
    for i in range(basissize):
        for j in range(basissize):
            result[i * basissize + j] = a[i] * b[j]
    return result


def twoE_orth(twoE_local_filename, orth_homolumo, syslen, output):
    twoE_local_2d_matrix = twoE_local_2d(twoE_local_filename, orth_homolumo, syslen)
    print twoE_local_2d_matrix.shape
    output = open(output, 'w')
    for i in range(syslen):
        I = orth_homolumo[:, i]
        for j in range(i, syslen):
            J = orth_homolumo[:, j]
            IJ = crossProduct(I, J)
            for x in range(syslen, 2*syslen):
                X = orth_homolumo[:, x]
                for y in range(x, 2*syslen):
                    Y = orth_homolumo[:, y]
                    XY = crossProduct(X, Y)
                    # print matrix(IJ).shape
                    oneTwoE = (matrix(IJ) * matrix(twoE_local_2d_matrix) * matrix(XY).transpose())[0,0]
                    output.write("%d\t%d\t%d\t%d\t\t%e\n" % (i, j, x-syslen, y-syslen, oneTwoE))
                    print "got one done"


if __name__ == "__main__":
    folder = "../g09Extract/data/porphyrin+DNA/porphyrin+AAATTT-A/"
    fockmatrix = load(folder + "fockmatrix.npy")
    homoevecs= load(folder + "homoevecs.npy")
    lumoevecs= load(folder + "lumoevecs.npy")
    orth_HOMOLUMO = load(folder + "orth_HOMOLUMO.npy")
    twoE_HOMOLUMOfilename = folder + "2EIntegral_all"
    syslen = 7

    [homoFock_orth, lumoFock_orth] = fockMatrix_orth(fockmatrix, homoevecs, lumoevecs, orth_HOMOLUMO, syslen)
    twoE_orth(twoE_HOMOLUMOfilename, orth_HOMOLUMO, syslen, folder + "twoE_orth")
    save(folder + "homoFock_orth.npy", homoFock_orth)
    save(folder + "lumoFock_orth.npy", lumoFock_orth)




