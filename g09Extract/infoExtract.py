__author__="Chaoren Liu"

from pylab import *
import re

def numOfBasisSet(fileObject):
	patern = "basis functions.*primitive gaussians"
	for eachline in fileObject:
		if re.search(patern, eachline):
			return int(eachline.split()[0])


def numOfAlphaAndBetaElectrons(fileObject):
	# return (a,b), a is the number of alpha electrons, b is the number of beta electrons
	patern = "alpha electrons"
	for eachline in fileObject:
		if re.search(patern, eachline):
			# return (int(eachline.split()[0]), int(eachline.split()[3]))
			return (int(eachline.split()[0]), int(eachline.split()[3]))



