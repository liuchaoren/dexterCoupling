'''
get the nlmos which are expanded in AOs, ordered from low energy to high energy
'''

import numpy as np
import re
from infoExtract import *
from constants import *


def aonlmos(logfile):
    sizeBasis = numOfBasisSet(logfile)
    # sizeEletron = numOfAlphaAndBetaElectrons(logfile)

    with  open(logfile) as log:
        for eachline in log:
            if re.search("NLMOs in the AO basis", eachline):
                break;

        nlmosNoFormat = []
        counter = 0
        numBasisTaken = 0
        numBasisLine = 0
        for eachline in log:
            ind = counter % (sizeBasis+3)
            counter += 1
            if ind == 0:
                numBasisTaken += numBasisLine
                if numBasisTaken == sizeBasis:
                    break
                continue
            elif ind == 2:
                continue
            elif ind == 1:
                numBasisLine = len(eachline.split()) - 1
            else:
                contents = eachline[16:]
                numitem = len(contents) / 8
                for i in range(numitem):
                    nlmosNoFormat.append(contents[i*8:(i+1)*8])

        reorder = []
        reorder += eachline[17:].split()
        for eachline in log:
            if re.search("Sorting of NBOs", eachline):
                reorder += eachline[17:].split()
            else:
                break

    assert len(nlmosNoFormat) == sizeBasis ** 2

    complete = (sizeBasis / 8) * 8
    residue = sizeBasis - complete

    nlmos=[]
    for i in range(sizeBasis):
        nlmos.append([])
    for i in range(sizeBasis**2):
        ind = i / (sizeBasis * 8)
        indexinBlock = i - ind * sizeBasis * 8
        if ind < complete / 8:
            belong = indexinBlock % 8 + ind * 8
        else:
            belong = indexinBlock % residue + ind * 8
        nlmos[belong].append(nlmosNoFormat[i])

    reorder = [int(i) - 1 for i in reorder]
    return np.array(nlmos, dtype='float').transpose()[:, reorder]


def aonlmosFchk(fchkfile):
    fchk = open(fchkfile)
    eachline = fchk.readline()
    while eachline[:21] != "Alpha MO coefficients":
        eachline = fchk.readline()

    N = int(eachline.split()[5])
    eachline = fchk.readline()
    count = 0
    allcoefs = []
    while count != N:
        linelen = len(eachline)
        itemnum = linelen / 16
        for i in range(itemnum):
            allcoefs.append(eachline[i*16:(i+1)*16])
            count+=1
        eachline = fchk.readline()
    numbasis = int(sqrt(N))
    return array(allcoefs,dtype='float').reshape(numbasis,numbasis).transpose()




def homolumos(nlmos, hindex):
    return nlmos[:, hindex]



if __name__ == "__main__":

    # for i in [4, 5, 6, 7, 8]:
    #     # filename = "data/alkyl/b12/symm/b12_opt.log"
    #     if i==13:
    #         fchkfile = "data/alkyl/b13/linear/symm/b13_opt.fchk"
    #         fockmatrix = np.load("data/alkyl/b13/linear/symm/fock.npy")
    #     else:
    #         fchkfile = "data/alkyl/b%s/symm/b%s_opt.fchk" % (i, i)
    #         fockmatrix = np.load("data/alkyl/b%s/symm/fock.npy" % i)
    for i in [1, 2, 3, 4, 5]:
        # filename = "data/alkyl/b12/symm/b12_opt.log"
        fchkfile = "data/alkyl/b7/folded/sample%s_opt.fchk" % i
        fockmatrix = np.load("data/alkyl/b7/folded/fock_%s.npy" % i)
        # filename = "data/alkyl/b12/symm/b12_opt.log"
        # fchkfile = "data/polynorbornyl/orthogonal/repeats_2/repeats_2_opt.fchk"
        # fockmatrix = np.load("data/polynorbornyl/orthogonal/repeats_2/fock.npy")

        nlmos = aonlmosFchk(fchkfile)

        # nlmos = aonlmos(filename)
        fockNLMOs=np.matrix(nlmos).transpose() * np.matrix(fockmatrix) * np.matrix(nlmos)
        # b12
        HOMOs = globals()["b7_HOMOs%s" % i]
        LUMOs = globals()["b7_LUMOs%s" % i]
        print fchkfile
        print HOMOs
        print LUMOs
        homos = homolumos(nlmos, HOMOs)
        lumos = homolumos(nlmos, LUMOs)
        homosfock = fockNLMOs[HOMOs][:, HOMOs]
        lumosfock = fockNLMOs[LUMOs][:, LUMOs]
        # np.save('data/alkyl/b12/symm/nlmos1.npy', nlmos)
        np.save('data/alkyl/b7/folded/homos_%s.npy' % i, homos)
        np.save('data/alkyl/b7/folded/lumos_%s.npy' % i, lumos)
        np.save('data/alkyl/b7/folded/homoFock_%s.npy' % i, homosfock)
        np.save('data/alkyl/b7/folded/lumoFock_%s.npy' % i, lumosfock)


