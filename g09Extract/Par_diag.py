#!/bin/py
# Partial diagonalize the Fock Matrix. Get the coupling as well as on-site energy
# conf1

import os
import sys
#sys.path.insert(1,'/opt/apps/sdg/lib/python2.7/site-packages')
import numpy as np
from numpy.linalg import inv
from numpy.linalg import norm
# from numpy import hstack
# import numpy.linalg as la
from scipy.linalg import block_diag
from scipy.linalg import eig


# hartree = 27.21138602
hartree = 1.
# transformation matrix, transform S to I. Beawere the transformation is not unitary transformation and resulting basis is not normalized. 
def A(S):
    val, vec = eig(S)
    return vec / np.sqrt(val)

def blockLoad(fullmatrix, sequence, basissize):
    basiscount=0
    segment = []
    unitnum = len(sequence)
    for i in range(unitnum):
        try:
            segment.append(fullmatrix[basiscount:basiscount+basissize[sequence[i]],basiscount:basiscount+basissize[sequence[i]]])
            basiscount=basiscount+basissize[sequence[i]]
        except KeyError as e:
            print "you sequence is wrong, or you need to update your dictionary"
            raise e

    return segment 

def homolumoIndex(sequence, basissize, electronsize):
    basiscount=0
    unitnum = len(sequence)
    homoindex = []
    lumoindex = []
    for i in range(unitnum):
        try:
            homoindex.append(basiscount+electronsize[sequence[i]]-1)
            lumoindex.append(basiscount+electronsize[sequence[i]])
            basiscount=basiscount+basissize[sequence[i]]
        except KeyError as e:
            print "your sequence is wrong, or you need to update your dictionary"

    return(homoindex, lumoindex)



def blockDiag(F, sequence, basissize):
    h_at_segment = blockLoad(F, sequence, basissize)
    evalall = np.array([])
    evecall = []
    for eachsegment in h_at_segment:
        eval_part, evec_part = eig(eachsegment)
        idx = eval_part.argsort()
        evalall = np.concatenate((evalall, eval_part[idx]))
        evecall.append(evec_part[:, idx])
    transformM = block_diag(*evecall)
    return (np.matrix(transformM).H * F * transformM, evalall, transformM)
    

if __name__ == "__main__":
    from constants import basissize_631g as basissize
    from constants import electronsize_631g as electronsize
    # input
    # sequence = ['T', 'T', 'T', 'T', 'T', 'T'] + ['C', 'C', 'C', 'A', 'A', 'A'] + ['A', 'A', 'A', 'G', 'G', 'G'] + ['porphyrin', 'porphyrin', 'porphyrin']
    sequence = ['A', 'T']
    # unitnum=len(sequence)
    totalbasis = sum([basissize[i] for i in sequence])
    path="data/AT/"


    # import fock matrix
    fockmatrix = np.load(path+"fock.npy")
    overlap = np.load(path+"overlap.npy")
    # overlapinv = inv(overlap)
    # h_at_segment = []
    # S_at_segment = []
    # homoindex=[]
    # lumoindex=[]
    S_at_segment = blockLoad(overlap, sequence, basissize)

    # get the diagonal unit matrix
    
    Strim = block_diag(*S_at_segment)
    A = A(Strim)
    # np.save("A", A)
    F_newBasis = np.matrix(A).H * fockmatrix * A 
    (blockDiaged, blockeval, blockevec) = blockDiag(F_newBasis, sequence, basissize)

    # np.save("blockeval", blockeval)
    blockevecTransformed = np.matrix(A)*blockevec 
    # blockevecTransformedNorm = np.array([norm(blockevecTransformed[:, i]) for i in range(totalbasis)])
    # blockevecTransformed = blockevecTransformed / blockevecTransformedNorm  # normalized
    (homoindex, lumoindex) = homolumoIndex(sequence, basissize, electronsize)
    homoFock = blockDiaged[np.array(homoindex)][:, np.array(homoindex)]
    lumoFock = blockDiaged[np.array(lumoindex)][:, np.array(lumoindex)]
    np.save(path+"homoFock.npy", homoFock.real)
    np.save(path+"lumoFock.npy", lumoFock.real)
    # np.save("blockevecTransformed", blockevecTransformed)
    # np.save("blockDiaged", blockDiaged)

    homoevecs = blockevecTransformed[:, np.array(homoindex)]
    lumoevecs = blockevecTransformed[:, np.array(lumoindex)]
    np.save(path+"homoevecs.npy", homoevecs.real)
    np.save(path+"lumoevecs.npy", lumoevecs.real)

