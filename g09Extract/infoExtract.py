__author__="Chaoren Liu"

from pylab import *
import re

def numOfBasisSet(filename):
	patern = "basis functions.*primitive gaussians"
	fileObject = open(filename)
	for eachline in fileObject:
		if re.search(patern, eachline):
			return int(eachline.split()[0])


def numOfAlphaAndBetaElectrons(filename):
	# return (a,b), a is the number of alpha electrons, b is the number of beta electrons
	patern = "alpha electrons"
	fileObject = open(filename)
	for eachline in fileObject:
		if re.search(patern, eachline):
			# return (int(eachline.split()[0]), int(eachline.split()[3]))
			return (int(eachline.split()[0]), int(eachline.split()[3]))



