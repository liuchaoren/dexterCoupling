'''
in two electron integral calculations, we calculated (i,j|t,s) where i, j, t, s can be any of HOMOs and LUMOs
'''

from LowdinPartition.hamiltonianConstructor import twoEmap



if __name__ == '__main__':
    dim = 7
    path = "data/porphyrin+DNA/porphyrin+AAATTT-T/"
    twoEfile = path + "2EIntegral_all"
    twoEall = twoEmap(twoEfile)
    partTwoE = {}
    for onekey in twoEall.keys():
        if (onekey[0] < dim and onekey[1] < dim and onekey[2] >= dim and onekey[3] >= dim):
            partTwoE[(onekey[0], onekey[1], onekey[2]-dim, onekey[3]-dim)] = twoEall[onekey];
        # twoEall.pop(onekey)


    partTwoEfile = open(path + "twoE_nonorth", 'w')

    for i in range(dim):
        for j in range(i, dim):
            for k in range(dim):
                for l in range(k, dim):
                    partTwoEfile.write("%d\t\t%d\t\t%d\t\t%d\t\t%e\n" % (i, j, k, l, partTwoE[(i, j, k, l)]))

