'''
convert between formats, or interfaces of input and output
Author - Chaoren Liu
Date - Jan 25, 2016
'''


from PyQuante.Molecule import Molecule

from numpy import *
import Queue

def gjf2Molecule(filename):

    elementsdict = {}
    elements=open('../resource/elements').readlines()

    for eachele in elements:
        (elename, eleindex) = eachele.split()
        elementsdict[elename] = eleindex

    gjf = open(filename, 'r').readlines()
    moleculeQueue = []
    filelen=len(gjf)
    index = filelen
    indexline = gjf[index-1]


    while len(indexline.split()) == 0:
        index = index - 1
        indexline = gjf[index-1]


    while not indexline.split()[0].isdigit():
        try:
            (ele, x, y, z) = indexline.split()
        except IndexError as e:
            print ("gjf file format error!")
            raise e

        eleindex = elementsdict[ele.lower()]
        moleculeQueue.append((int(eleindex), (float(x),float(y),float(z))))

        index = index - 1
        indexline = gjf[index-1]

    (chg, multi) = indexline.split()


    moleculelist = []
    while len(moleculeQueue) != 0:
        moleculelist.append(moleculeQueue.pop())

    return Molecule(filename, moleculelist, units='Angstrom', charge=int(chg))




