'''
test the 2e integrals of benzene

'''

# from PyQuante.pyints import contr_coulomb as pycc
# from PyQuante.rys import contr_coulomb as ryscc
# from PyQuante.hgp import contr_coulomb as hgpcc
# from PyQuante.cints import contr_coulomb as ccc
from PyQuante.crys import contr_coulomb as cryscc
# from PyQuante.chgp import contr_coulomb as chgpcc
from PyQuante.Ints import getbasis
from PyQuante.Molecule import Molecule

from PyQuante2E.converters import *
from PyQuante2E.twoEInts import *


gjffilename='../data/2benzene/2benzene_symm.gjf'

if __name__ == "__main__":
    myMol = gjf2Molecule(gjffilename)
    # print myMol.atoms
    bfs = getbasis(myMol, '6-31g')

    print coulomb(bfs[73], bfs[39], bfs[38], bfs[27], cryscc)






