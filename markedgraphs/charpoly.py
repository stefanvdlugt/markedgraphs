import numpy as np

def charpoly(A):
    """
    Takes a matrix A with integer coefficients and outputs the coefficients of its characteristic polynomial.
    Uses the Faddeev--LeVerrier algorithm.
    """
    n=len(A)
    M=np.zeros((n,n),dtype=int)
    AM=np.dot(A,M)
    I=np.identity(n,dtype=int)
    c=[1]
    for k in range(1,n+1):
        M=AM+c[-1]*I
        AM=np.dot(A,M)
        c.append(-np.trace(AM)//k)
    return tuple(c)
