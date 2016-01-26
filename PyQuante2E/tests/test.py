'''
test the 2e integrals of benzene

'''



# from PyQuante.pyints import contr_coulomb as pycc
# from PyQuante.rys import contr_coulomb as ryscc
# from PyQuante.hgp import contr_coulomb as hgpcc
# from PyQuante.cints import contr_coulomb as ccc
from PyQuante.crys import contr_coulomb as cryscc
# from PyQuante.chgp import contr_coulomb as chgpcc
# from PyQuante.CGBF import CGBF
# from PyQuante.cints import ijkl2intindex
from PyQuante.Ints import getbasis
from PyQuante.Molecule import Molecule

from src.converters import *
from src.twoEInts import *


gjffilename='../tests/benzene.gjf'

tmpname = '../tests/tmp'

if __name__ == "__main__":
    myMol = gjf2Molecule(gjffilename)
    # print myMol.atoms
    bfs = getbasis(myMol, '6-31g')

    tmp = open(tmpname).readlines()
    newcal = open('../tests/tmp_new', 'w')

    for eachline in tmp[:]:
        (i, j, k, l) = (int(eachline[3:6])-1, int(eachline[9:12])-1, int(eachline[15:18])-1, int(eachline[21:24])-1)

        # newcal.write('%d\t%d\t%d\t%d\t%20.12e\n' % (i, j, k, l, coulomb(bfs[i], bfs[j], bfs[k], bfs[l], cryscc)))
        coulomb(bfs[i], bfs[j], bfs[k], bfs[l], cryscc)



    # print '%20.12e' % coulomb(bfs[12], bfs[1], bfs[1], bfs[1], cryscc)



