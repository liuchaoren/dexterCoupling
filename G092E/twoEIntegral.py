__author__ = "Chaoren Liu"
'''
Given a two-electron integral file (output from g09), provide matrix (IJ|KL). 
'''

import numpy as np
import sys
import itertools
from scipy.sparse import csr_matrix


#twoElectron_file = 'part' + sys.argv[1]
twoElectron_file = "twoElectron"
basisSetSize = 126
trunksize = 100000
syslen=9

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
    indexUnique = [(m1, n1)]
    if not (m2, n2) in indexUnique:
        indexUnique.append((m2,n2)) 
    if not (m3, n3) in indexUnique:
        indexUnique.append((m3,n3)) 
    if not (m4, n4) in indexUnique:
        indexUnique.append((m4,n4)) 
    if not (m5, n5) in indexUnique:
        indexUnique.append((m5,n5)) 
    if not (m6, n6) in indexUnique:
        indexUnique.append((m6,n6)) 
    if not (m7, n7) in indexUnique:
        indexUnique.append((m7,n7)) 
    if not (m8, n8) in indexUnique:
        indexUnique.append((m8,n8)) 
    
    rows = [] 
    cols = []
    for (m, n) in indexUnique:
        rows.append(m)
        cols.append(n)
    numOfUnique = len(rows)

    return csr_matrix((np.array([value]*numOfUnique), (np.array(rows), np.array(cols))), shape=(basisSetSize**2, basisSetSize**2))

def maptoSparseMatrixAll(line):
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
    indexUnique = [(m1, n1)]
    if not (m2, n2) in indexUnique:
        indexUnique.append((m2,n2)) 
    if not (m3, n3) in indexUnique:
        indexUnique.append((m3,n3)) 
    if not (m4, n4) in indexUnique:
        indexUnique.append((m4,n4)) 
    if not (m5, n5) in indexUnique:
        indexUnique.append((m5,n5)) 
    if not (m6, n6) in indexUnique:
        indexUnique.append((m6,n6)) 
    if not (m7, n7) in indexUnique:
        indexUnique.append((m7,n7)) 
    if not (m8, n8) in indexUnique:
        indexUnique.append((m8,n8)) 
    
#    rows = [] 
#    cols = []
#    numOfUnique = len(indexUnique)
    for (m, n) in indexUnique:
        rows.append(m)
        cols.append(n)
        values.append(value)

#    return csr_matrix((np.array([value]*numOfUnique), (np.array(rows), np.array(cols))), shape=(basisSetSize**2, basisSetSize**2))

def maptoFullMatrix(line):
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
    
    indexUnique = [(m1, n1)]
    if not (m2, n2) in indexUnique:
        indexUnique.append((m2,n2)) 
    if not (m3, n3) in indexUnique:
        indexUnique.append((m3,n3)) 
    if not (m4, n4) in indexUnique:
        indexUnique.append((m4,n4)) 
    if not (m5, n5) in indexUnique:
        indexUnique.append((m5,n5)) 
    if not (m6, n6) in indexUnique:
        indexUnique.append((m6,n6)) 
    if not (m7, n7) in indexUnique:
        indexUnique.append((m7,n7)) 
    if not (m8, n8) in indexUnique:
        indexUnique.append((m8,n8)) 
    
    for (m, n) in indexUnique:
        fullMatrix[m,n] = value
#       fullMatrix[m,n] = fullMatrix[m, n] + 1
#    fullMatrix[m2,n2] = fullMatrix[m2, n2] + 1
#    fullMatrix[m3,n3] = fullMatrix[m3, n3] + 1
#    fullMatrix[m4,n4] = fullMatrix[m4, n4] + 1
#    fullMatrix[m5,n5] = fullMatrix[m5, n5] + 1
#    fullMatrix[m6,n6] = fullMatrix[m6, n6] + 1
#    fullMatrix[m7,n7] = fullMatrix[m7, n7] + 1
#    fullMatrix[m8,n8] = fullMatrix[m8, n8] + 1

def crossProduct(a, b):
    result = np.outer(a,b)
    return np.reshape(result, (1, len(a)*len(b)))


def maptoPartialResult(t):
    return np.dot((leftVector * t), rightVector.reshape(basisSetSize**2, 1))[0]
    # return np.dot(np.dot(leftVector, t), np.transpose(rightVector))

def reducetoResult(a, b):
    return a+b

def four2twoIndex(i,j,k,l):
    return ((i-1) * basisSetSize + (j-1), (k-1) * basisSetSize + (l-1))

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
    twoEfile = open(twoElectron_file, 'r')
    print "I am working on %s" % twoElectron_file
    homoevecs = np.load("homos.npy")
    lumoevecs = np.load("lumos.npy")
    ijkloutput = {}
    counter = trunksize
    while counter == trunksize:
        rows = [] 
        cols = []
        values = []
    #    fullMatrix = np.zeros([basisSetSize**2, basisSetSize**2], dtype=float)
        # duplicateMonitor = 0
    #    sumall = np.matrix([[0.]])
        counter = 0
        for eachline in twoEfile:
            maptoSparseMatrixAll(eachline)
            counter = counter + 1
            if counter == trunksize:
                break
        # print "rows, cols, values were constructed"
        matrixSparse = csr_matrix((np.array(values), (np.array(rows), np.array(cols))), shape=(basisSetSize**2, basisSetSize**2))
        for i in range(syslen):
            veci = homoevecs[:, i]
            for j in range(i+1):
                vecj = homoevecs[:, j]
                leftvec = crossProduct(veci, vecj)
                for k in range(syslen):
                    veck = lumoevecs[:, k]
                    for l in range(k+1):
                        vecl = lumoevecs[:, l]
                        rightvec = crossProduct(veck, vecl)
                        oneijkl = (np.matrix(leftvec) * matrixSparse  * np.matrix(rightvec).transpose())[0,0]
                        if ijkloutput.has_key((i,j,k,l)):
                            ijkloutput[(i,j,k,l)] = ijkloutput[(i,j,k,l)] + oneijkl
                        else:
                            ijkloutput[(i,j,k,l)] = oneijkl
        print "one trunk done"
    
# output 
    outputname = open("twoE"+twoElectron_file, "w")
    for (i,j,k,l), value in ijkloutput.iteritems():
        outputname.write("%8d\t%8d\t%8d\t%8d\t\t%10.6e\n" % (i,j,k,l,value))
                    

