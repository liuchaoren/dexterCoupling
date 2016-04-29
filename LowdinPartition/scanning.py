'''
scan the energy to get a delocalized D/A state
'''
from numpy import *
from numpy.linalg import eig
hartree = 27.21138602

def initialOpt(H, iteratornum):
    EDAstart = (H[0,0]+H[-1,-1])/2
    EDnow = EDAstart
    EAnow = EDAstart
    for i in range(iteratornum):
        H[0,0] = EDnow
        H[-1,-1] = EAnow
        (eval, evec) = eig(H)
        (minindex, secondminindex) = eval.argsort()[:2]
        if (abs(evec[0,minindex]) > abs(evec[-1,minindex])):
            EDsplit = eval[minindex]
            EAsplit = eval[secondminindex]
        else:
            EDsplit = eval[secondminindex]
            EAsplit = eval[minindex]
        error = (EDsplit - EAsplit) / 2
        EDnow = EDnow - error
        EAnow = EAnow + error
        print error
        if abs(error) < 2e-9:
            return (EDnow,EAnow)
    return (EDnow, EAnow)

def fineTunning(H, EDnow, EAnow, energyArray):
    estimator = 10
    for eachE in energyArray:
        H[0,0] = EDnow + eachE
        H[-1,-1] = EAnow - eachE
        (eval, evec) = eig(H)
        (minindex, secondminindex) = eval.argsort()[:2]
        trialEstimator = abs(abs(evec[0, minindex]) -  abs(evec[-1, minindex]))
        if trialEstimator < estimator:
            estimator = trialEstimator
            finalED = H[0,0]
            finalEA = H[-1,-1]
            # evec1 = evec[:,0]
            # evec2 = evec[:,1]
            # saveEval = eval
            print estimator
    # print evec1, evec2
    # print "coupling is", (saveEval[0]-saveEval[1])/2
    return (finalED, finalEA)

def combineOpt(H, iteratornum, energyArray):
    (EDnow, EAnow) = initialOpt(H, iteratornum)
    # print "I am working on fine scanning"
    (finalED, finalEA) = fineTunning(H, EDnow, EAnow, energyArray)
    H[0,0] = finalED
    H[-1,-1] = finalEA
    savetxt("HBE_b7_lower0.19.npy", H, fmt='%24.16e')
    (eval,evec) = eig(H)
    (minindex, secondminindex) = eval.argsort()[:2]
    assert(abs(evec[0,minindex]) > 0.6 and abs(evec[-1,minindex])>0.6 and abs(evec[0, secondminindex])>0.6 and abs(evec[-1,secondminindex])>0.6)
    return abs((eval[minindex] - eval[secondminindex])/2)

if __name__ == "__main__":

    for i in range(2,3):
        H = load('H%s_0.0.npy' % i)
        HBE = load('HBE%s_0.0.npy' % i)
        HCT = load('HCT%s_0.0.npy' % i)
        dim = int(sqrt(H.shape[0]))
    # H=CTExcitonHamiltonian(H, dim)

        iteratornum = 2000
        print "start to calculate H"
        (EDnow, EAnow) = initialOpt(H, iteratornum)
        scanarray=arange(-1e-9, 1e-9, 1e-12)
        (finalED, finalEA) = fineTunning(H, EDnow, EAnow, scanarray)
        H[0,0] = finalED
        H[-1,-1] = finalEA
        # save("data/alkyl/b13/folded/sample2/H_allExchanges_opt.npy", H)
        (eval,evec) = eig(H)
        print 'total coupling is', (eval[1] - eval[0])/2 * hartree
        print finalED, finalEA
        print evec[0,0], evec[-1,0], evec[0,1], evec[-1,1]

        print "start to calculate HBE"
        (EDnow, EAnow) = initialOpt(HBE)
        scanarray=arange(-1e-9, 1e-9, 1e-12)
        (finalED, finalEA) = fineTunning(HBE, EDnow, EAnow, scanarray)
        HBE[0,0] = finalED
        HBE[-1,-1] = finalEA
        # save("data/alkyl/b13/folded/sample2/H_allExchanges_opt.npy", H)
        (eval,evec) = eig(HBE)
        print 'BE coupling is', (eval[1] - eval[0])/2 * hartree
        print finalED, finalEA
        print evec[0,0], evec[-1,0], evec[0,1], evec[-1,1]

        print "start to calculate HCT"
        (EDnow, EAnow) = initialOpt(HCT)
        scanarray=arange(-1e-9, 1e-9, 1e-12)
        (finalED, finalEA) = fineTunning(HCT, EDnow, EAnow, scanarray)
        HCT[0,0] = finalED
        HCT[-1,-1] = finalEA
        # save("data/alkyl/b13/folded/sample2/H_allExchanges_opt.npy", H)
        (eval,evec) = eig(HCT)
        print 'CT coupling is', (eval[1] - eval[0])/2 * hartree
        print finalED, finalEA
        print evec[0,0], evec[-1,0], evec[0,1], evec[-1,1]


