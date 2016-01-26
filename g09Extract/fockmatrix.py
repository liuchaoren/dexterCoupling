__author__="Chaoren Liu"
'''
extract fock matrix and overlap from g09 dumped file
'''

from pylab import *
import re
from overlap import StartLine, matrixFromLinearFormat
from infoExtract import numOfBasisSet

fockfile = "data/2ethylene/distance4.5/fock"
# outputfile = "data/pPt+AT/pPt+AT.log"
fockmatrixOutput = "data/2ethylene/distance4.5/fock.npy"


if __name__ == "__main__":
	fockFileObject = open(fockfile, "r").readlines()
	# outputFileObject = open(outputfile, 'r').readlines()
	forckPatern = "Dump"
	startLineNum = StartLine(forckPatern, fockFileObject)
	linearFock = []
	for eachline in fockFileObject[startLineNum+1:]:
		linearFock = linearFock + eachline.replace("D", "E").split()
	# numOfBasisSet = numOfBasisSet(outputFileObject)
	# print numOfBasisSet
	numOfBasisSet = 26 * 2
	print numOfBasisSet

	matrix = matrixFromLinearFormat(numOfBasisSet, linearFock)
	save(fockmatrixOutput, matrix)

