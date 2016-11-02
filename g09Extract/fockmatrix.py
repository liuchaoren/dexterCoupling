__author__="Chaoren Liu"
'''
extract fock matrix and overlap from g09 dumped file
'''

from numpy import *
import re
from os import system
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
		fockfile = "data/porphyrin+DNA/porphyrin+AAATTT-T/fockmatrix"
		# nb = 874
		fockmatrixoutput = fockfile + ".npy"
		system(r'sed "s/D/E/g" ' + fockfile + r' > ' + fockfile + r'_tmp')
		fockfile_tmp = fockfile + '_tmp'
		# outputfile = "data/alkyl/b%s/linear/symm/b%s_opt.log" % (i, i)

		# fockfile = "data/polynorbornyl/orthogonal/repeats_2/fockmatrix"
		# # outputfile = "data/alkyl/b%s/linear/symm/b%s_opt.log" % (i, i)
		# fockmatrixoutput = "data/polynorbornyl/orthogonal/repeats_2/fock.npy"
		# print fockfile
		# print outputfile
		# print fockmatrixoutput

		fockFileObject = open(fockfile, "r").readlines()
		# outputFileObject = open(outputfile, 'r').readlines()
	 	(startLineNum, elementNum) = StartLine(fockFileObject)
		nb = int((-1 + sqrt(1+8*elementNum))/2)
		linearFock = ravel(genfromtxt(fockfile_tmp, dtype='float', skip_header=startLineNum, skip_footer=1)).tolist()
		linearFock = linearFock + array(fockFileObject[-1].replace("D", "E").split(), dtype='float').tolist()
		# print linearFock[-10:]

		# linearFock = []
		# for eachline in fockFileObject[startLineNum+1:]:
		# 	linearFock = linearFock + eachline.replace("D", "E").split()
			# linearFock = linearFock + eachline.split()
			# print "one"
		# nb = numOfBasisSet(outputfile)
		# nb = basissize["alkane%s" % i]
		assert(elementNum == len(linearFock))
		# numOfBasisSet = 204
		# print nb
		# numOfBasisSet = 26 * 2
		# print numOfBasisSet

		# linearFock = load("data/3wj/zindo/fock1.npy")
		# print linearFock
		matrix = matrixFromLinearFormat(nb, linearFock)
		save(fockmatrixoutput, matrix)
		system("rm " + fockfile_tmp)

