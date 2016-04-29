__author__="Chaoren Liu"
'''
extract fock matrix and overlap from g09 dumped file
'''

from numpy import *
import re
from overlap import StartLine, matrixFromLinearFormat
from infoExtract import numOfBasisSet
from constants import basissize_sto3g as basissize



if __name__ == "__main__":
	# for i in [4, 5, 6, 7, 8]:
	# 	if i==13:
	# 		fockfile = "data/alkyl/b%s/linear/symm/fockmatrix" % i
	# 		# outputfile = "data/alkyl/b%s/linear/symm/b%s_opt.log" % (i, i)
	# 		fockmatrixoutput = "data/alkyl/b%s/linear/symm/fock.npy" % i
	# 	else:
	# 		fockfile = "data/alkyl/b%s/symm/fockmatrix" % i
	# 		# outputfile = "data/alkyl/b%s/symm/b%s_opt.log" % (i, i)
	# 		fockmatrixoutput = "data/alkyl/b%s/symm/fock.npy" % i
	# for i in [1, 2, 3, 4, 5]:
		fockfile = "data/AT/overlap"
		# outputfile = "data/alkyl/b%s/linear/symm/b%s_opt.log" % (i, i)
		fockmatrixoutput = "data/AT/overlap.npy"

		# fockfile = "data/polynorbornyl/orthogonal/repeats_2/fockmatrix"
		# # outputfile = "data/alkyl/b%s/linear/symm/b%s_opt.log" % (i, i)
		# fockmatrixoutput = "data/polynorbornyl/orthogonal/repeats_2/fock.npy"
		# print fockfile
		# print outputfile
		# print fockmatrixoutput

		fockFileObject = open(fockfile, "r").readlines()
		# outputFileObject = open(outputfile, 'r').readlines()
		fockPatern = "Dump"
		startLineNum = StartLine(fockPatern, fockFileObject)
		linearFock = []
		for eachline in fockFileObject[startLineNum+1:]:
			linearFock = linearFock + eachline.replace("D", "E").split()
			# linearFock = linearFock + eachline.split()
			# print "one"
		# nb = numOfBasisSet(outputfile)
		# nb = basissize["alkane%s" % i]
		nb = 193
		# numOfBasisSet = 204
		# print nb
		# numOfBasisSet = 26 * 2
		# print numOfBasisSet

		# linearFock = load("data/3wj/zindo/fock1.npy")
		# print linearFock
		matrix = matrixFromLinearFormat(nb, linearFock)
		save(fockmatrixoutput, matrix)

