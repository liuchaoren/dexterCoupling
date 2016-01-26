__author__ = "Chaoren Liu"
'''
Given a two-electron integral file (output from g09) and the four MOs I, J, K, L (vector of AO), provide (IJ|KL). 
'''

from pyspark import SparkContext
#from functools import partial
#from collections import defaultdict
#import StringIO
#import csv
#import datetime
#import permute.interval as interval
#import pandas
import numpy as np
from scipy.sparse import csr_matrix
from parse import *
from numpy.random import rand
import itertools


twoElectron_file = 'twoElectronIntegral_part'
master = "local[4]"
numofPartitions = 100
# parser = "I={} J={} K={} L={} Int={}" # parse the twoElectron file
basisSetSize = 100

# I, J, K, L coefficient 
# coeff = np.genfromtxt("file")
# I = coeff[:, x]
# J = coeff[:, y]
# K = coeff[:, z]
# L = coeff[:. t]

I = (rand(basisSetSize) - 0.5)
I = I/np.linalg.norm(I)
J = (rand(basisSetSize) - 0.5)
J = J/np.linalg.norm(J)
K = (rand(basisSetSize) - 0.5)
K = K/np.linalg.norm(K)
L = (rand(basisSetSize) - 0.5)
L = L/np.linalg.norm(L)

leftVector = itertools.product(np.conjugate(I), J)
rightVector = itertools.product(np.conjugate(K), L)
t=[]
for i in leftVector:
    t.append(i[0]*i[1])
leftVector = np.array(t)
t = []
for i in rightVector:
    t.append(i[0]*i[1])
rightVector = np.array(t)


#master = "spark://ec2-54-160-210-33.compute-1.amazonaws.com:7077"

def load_twoElectron(context):
    twoElectronIntegralInAO = context.textFile(twoElectron_file, use_unicode=False).cache()
    return twoElectronIntegralInAO

def maptoSparseMatrix(line):
    i = int(line[3:6])
    j = int(line[9:12])
    k = int(line[15:18])
    l = int(line[21:24])
    value = float(line[29:].replace("D", "E"))
    # return parse(parser, line)
    # (i,j,k,l,value) = parse(parser, line)

    # i=int(i)
    # j=int(j)
    # k=int(k)
    # l=int(l)
    (m1, n1) = four2twoIndex(i,j,k,l)
    (m2, n2) = four2twoIndex(k,l,i,j)
    (m3, n3) = four2twoIndex(j,i,l,k)
    (m4, n4) = four2twoIndex(l,k,j,i)
    (m5, n5) = four2twoIndex(j,i,k,l)
    (m6, n6) = four2twoIndex(l,k,i,j)
    (m7, n7) = four2twoIndex(i,j,l,k)
    (m8, n8) = four2twoIndex(k,l,j,i)

    return csr_matrix((np.array([value]*8), (np.array([m1, m2, m3, m4, m5, m6, m7, m8]), np.array([n1, n2, n3, n4, n5, n6, n7, n8]))), shape=(basisSetSize**2, basisSetSize**2))

def maptoPartialResult(t):
    return np.dot((leftVector * t), rightVector.reshape(basisSetSize**2, 1))[0]
    # return np.dot(np.dot(leftVector, t), np.transpose(rightVector))

def reducetoResult(a, b):
    return a+b

def four2twoIndex(i,j,k,l):
    return ((i-1) * basisSetSize + (j-1), (k-1) * basisSetSize + (l-1))

def func(i):
    print i



def run(context):
    """ Data is in the following format: (bill_id, person_id, vote, type, chamber, year, date, session, status, extra).
        1. Key everything by bill_id, parse the date correctly, return (bill_id, (person_id, vote, date))
        2. Self join bills RDD to itself.
        3. Remove duplicate entries (via left.person_id < right.person_id) and comparisons with self.
        4. Map the joined data to a (person_id:person_id, (agree, disagree)) RDD.
        5. Iterate through the intervals, calculating (agree/(agree + disagree)) for each key in 4.
        6. Store a dictionary, for each key, store the lowest current percent and interval it was found in.
    """
    raw_twoElectron = load_twoElectron(context)
#    intervals = interval.interval_set('1/1/2012', '1/1/2014', freq='15D', max_delta=pandas.Timedelta(days=120))
    partialMatrix = raw_twoElectron.map(maptoSparseMatrix)
#    joined = bills.join(bills, 32)

    partialResults = partialMatrix.map(maptoPartialResult)
    result = partialResults.reduce(reducetoResult)
    return result
    # return partialMatrix
    # return partialResults

if __name__ == "__main__":
    # print leftVector
    # print rightVector
    context = SparkContext(master, "twoElectronCalculation")
    # run(context).foreach(func)
    result = run(context)
    print result
#    output = open('output.txt', 'w')
#    for result in results:
#        output.write("%s\t(%s\t\t%d\t\t%d)\n" % (result[0], result[1][0], result[1][1], result[1][2]))
#    output.close()
