import numpy as np
def Dn(n,N):
    D = np.zeros(N+1)
    for p in np.arange(n + 1, N + 1):
        if (p + n) % 2 == 1:
            D[p] = 2 * p
    return D