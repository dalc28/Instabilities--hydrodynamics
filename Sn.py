import numpy as np
def Sn(n,u2,N):
    s1 = np.zeros(N+1)
    s2 = np.zeros(N+1)
    s3 = np.zeros(N+1)
    for m in np.arange(n + 1):
        s1[m] = 0.5 * u2[n - m]
    for m in np.arange(n + 1, N + 1):
        s2[m] = 0.5 * u2[m - n]
    for m in np.arange(1, N - n + 1):
        s3[m] = 0.5 * u2[n + m]
    sn = s1 + s2 + s3
    return sn