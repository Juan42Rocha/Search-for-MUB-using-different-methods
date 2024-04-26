import numpy as np
import cmath



def ChechSol(fname,d,n,k):
    #"CsDHoHu2n1k4itr1state.txt"
    st = np.genfromtxt(fname)
    qbv = k**d
    sts = np.reshape(st,(n*d,qbv))


    for icont in range(n*d):
        aux = np.sum(sts[icont])
        if (aux !=  qbv-2):
            print("No valid sol\n")

    vector = []

            