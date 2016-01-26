__author__="Chaoren Liu"
'''
extract overlap matrix from g09 output file
'''

from pylab import *
import re
from infoExtract import numOfBasisSet, numOfAlphaAndBetaElectrons

outputfile="pPt+AATT.log"

def StartLine(patern, fileObject):
# find the line number where overlap matrix starts, count from 0
	counter = 0
	for eachline in fileObject:
		if re.search(patern, eachline):
			return counter
		counter = counter + 1

def OverlapLineNum(basisNum):
# sum is the number of lines which contains overlap matrix, rippedList is the line index which should be removec
	sum = 0
	rippedList = [sum]
	i = basisNum
	while i != 0:
		sum = sum + (i+1)
		rippedList.append(sum)
		i = i - 5
	return (sum, rippedList)


# outputtmp = open("overlap", 'w')
# for eachline in OverlapMatrixLinearFormat:
# 	outputtmp.write("%f\n" % eachline)

def matrixFromLinearFormat(basisNum, MatrixLinearFormat):
	# convert a linear format symmetric matrix (only the lower part is stored) into 2d matrix
	matrix = zeros([basisNum, basisNum])
	counter = 0
	for i in range(basisNum):
		for j in range(i+1):
			matrix[i, j] = MatrixLinearFormat[counter]
			counter = counter + 1
	for j in range(basisNum):
		for i in range(j+1):
			matrix[i, j] = matrix[j, i]
	return matrix

if __name__ == "__main__":
	fileObject = open(outputfile,'r').readlines()
	basisNum = numOfBasisSet(fileObject)

	numofElectrons = sum(numOfAlphaAndBetaElectrons(fileObject))
	overlapPatern = "\*\*\* Overlap \*\*\*"
	startLineNum = StartLine(overlapPatern, fileObject) + 1
	(totalNum, rippedList) = OverlapLineNum(basisNum)


	OverlapMatrixLines = array(fileObject[startLineNum:startLineNum+totalNum])
	filter = array([True] * totalNum)
	for eachindex in rippedList[:-1]:
		filter[eachindex] = False

	OverlapMatrixLines = OverlapMatrixLines[filter]

	# outputtmp = open("overlap", 'w')
	# for eachline in OverlapMatrixLines:
	# 	outputtmp.write(eachline)

	OverlapMatrixLinearFormat = []
	for eachline in OverlapMatrixLines:
		OverlapMatrixLinearFormat = OverlapMatrixLinearFormat + eachline.replace("D", "E").split()[1:]

	OverlapMatrixLinearFormat = array(OverlapMatrixLinearFormat, dtype='float')

	overlapMatrix = matrixFromLinearFormat(basisNum, OverlapMatrixLinearFormat)
	save("overlap.npy", overlapMatrix)

