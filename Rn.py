import numpy as np
def Rn(n,u,N):
    rn1 = np.zeros(N+1)
    rn2 = np.zeros(N+1)
    rn3 = np.zeros(N+1)
    for m in np.arange(n+1):
        rn1[m] = 0.5*u[n-m]
    for m in np.arange(n+1,N+1):
        rn2[m] = 0.5*u[m-n]
    for m in np.arange(1,N-n+1):
        rn3[m] = 0.5*u[n+m]
    rn = rn1 + rn2 +rn3

    return rn