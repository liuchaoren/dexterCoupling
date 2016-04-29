'''
get the order of orbitals
'''

from numpy import *


from constants import b13_HOMOs1 as HOMOs
from constants import b13_LUMOs1 as LUMOs

nlmosref = load('data/alkyl/b13/folded/nlmos1.npy')
nlmos = load('data/alkyl/b13/folded/nlmos2.npy')
basisnum = nlmos.shape[0]

HOMOsnew = []
for each in HOMOs:
    reforbital = nlmosref[:, each]
    init=0
    mirrororbital = -1
    for i in range(basisnum):
        targetoribital = nlmos[:, i]
        product = abs(reforbital.dot(targetoribital))/linalg.norm(reforbital)/linalg.norm(targetoribital)
        if product > init:
            init = product
            mirrororbital = i

    HOMOsnew.append(mirrororbital)

print HOMOsnew


