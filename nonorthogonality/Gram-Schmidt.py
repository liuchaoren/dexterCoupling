"""
using gram-schmidt orthogonality, make HOMOs and LUMOs (from block diagonalization method)
"""

from numpy import *


# the order is D, A, Bs, and HOMO first and then LUMO
def schmidt_orth(overlap, order, len):
    orthoBasis = zeros((len*2, len*2))
    orthoBasis[0,0] = 1
    for i in range(len*2)[1:]:
        currBasis = zeros(len*2)
        currBasis[order[i]] = 1
        nextOrth = currBasis.copy()
        for j in range(i):
            oneOrth = orthoBasis[:,j]
            nextOrth -= (((matrix(currBasis) * matrix(overlap) * matrix(oneOrth).transpose())[0,0]) / ((matrix(oneOrth) * overlap * matrix(oneOrth).transpose())[0,0]) * oneOrth)
        nextOrth = nextOrth / sqrt((matrix(nextOrth) * matrix(overlap) * matrix(nextOrth).transpose())[0,0])
        orthoBasis[:, i] = nextOrth

    idx = argsort(order)
    orthoBasis = orthoBasis[:,idx]
    return orthoBasis

def overlap_HOMOLUMO(homoevecs, lumoevecs, overlap):
    homospluslumos = hstack((homoevecs, lumoevecs))
    return matrix(homospluslumos.transpose()) * matrix(overlap) * matrix(homospluslumos)



if __name__ == "__main__":

    folder = "../g09Extract/data/porphyrin+DNA/porphyrin+AAATTT-T/"
    overlap = load(folder + "overlap.npy")
    homoeves = load(folder + "homoevecs.npy")
    lumoeves = load(folder + "lumoevecs.npy")
    overlap_HOMOLUMO = overlap_HOMOLUMO(homoeves, lumoeves, overlap)
    len = 7
    DAorder = [0, len, len-1, len*2-1]
    Bsorder = [[i, i+len] for i in range(1, len-1)]
    order = DAorder
    for each in Bsorder:
        order = order + each

    print order

    orth = schmidt_orth(overlap_HOMOLUMO, order, len)
    # orth is orthogonal homos and lumos in the original order, in the basis of the localized homos and lumos
    save(folder + "orth_HOMOLUMO.npy", orth)




