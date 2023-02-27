#Análisis espacial de inestabilidades hidrodinámicas
##Software diseñado por Diego Armando Landinez Capacho
##Universidad de Valle, colombia,2023.
## Paquetes a usar
import numpy as np
import scipy.linalg as la
import numpy.polynomial.chebyshev as npcheby
import pandas as pd
import matplotlib.pyplot as plt
from Dn import Dn
from Qn import Qn
from Rn import Rn
from Sn import Sn
from U import U
N = 100
#alphar = 0.925761
bethar = 0.2
#yinf = (2*np.pi)/alphar
yinf = 6
ymax = yinf
ymin = -ymax
#Aplicamos una transformación para mapear los valores de y a los de z que corresponden al dominio de los polinomios de
#Chebyshev
# r es el factor de escala de la transformación recomendado como 2.0 para el método de colocación.
r = 2.0
y = np.linspace(ymin,ymax,N+1)
z = y/np.sqrt((r**2)+(y**2))
zmax = max(z)
zmin = -zmax
m = ((1-(z**2))**(3/2))/r
# %% Interpolación de Chebyshev del perfil de velocidades
vel_profile_cheby = npcheby.Chebyshev.interpolate(U, N, [zmin,zmax])
u = vel_profile_cheby.coef

# %% Interpolación de Chebychev para la segunda derivada del perfil de velocidades

vel_coef_d2 = npcheby.chebder(vel_profile_cheby.coef,2)
d2_list = vel_coef_d2.tolist(); d2_tuple = tuple(d2_list)
vel_profile_d2_cheby = npcheby.Chebyshev(d2_list, [ymin,ymax])
u2 = np.append(vel_coef_d2, np.zeros(N+1 - len(vel_coef_d2)))

## Determinación de las matrices

qn = np.zeros((N+1,N+1))
rn = np.zeros((N+1,N+1))
sn = np.zeros((N+1,N+1))
ln = np.zeros((N+1,N+1))
D = np.zeros((N+1,N+1))

for n in np.arange(N+1):
    ##Matriz L_n ecuación 42.1
    D[n] = Dn(n,N,m)
    ##Matriz Rn ecuación 42.3
    rn[n] = Rn(n,u,N)
    ##Matriz Sn ecuación 42.4
    sn[n] = Sn(n,u2,N)
ln = D @ D


for n in np.arange(N+1):
    qn[n] = Qn(ln,n,u,N)


#Condiciones de frontera

# J_1: ϕ(y_inf) = 0
J_1 = np.ones(N+1)

# J_2: ϕ(-y_inf) = 0
J_2 = np.ones(N+1)
J_2[1::2] = -1
'''
# J_3: Slope of ϕ(y) at y = y_inf
J_3 = np.sum(D, axis=0)

# J_4: Slope of ϕ(y) at y = -y_inf
J_4 = np.zeros(N+1)

for j in range(0,N+1):
    for i in range(0,N+1):
       J_4[j] += D[i,j] * (-1)**i
'''
#J_5: ϕ'(y=y_inf) = -αϕ
#J_5 = J_3 + alphar
# J_6: ϕ'(y=-y_inf) = αϕ
#J_6 = J_4 - alphar

# Construcciones de las matrices C
C0 = np.zeros([N+1,N+1])
C0 = (-1)*rn
C0i = la.inv((-1)*rn)
C1 = (-bethar)*(np.identity(N+1))
C2 = (sn-qn)
C3 = (bethar)*ln
# incorporar las condiciones de frontera a las matrices
C3[N-1] = J_1
C3[N] = J_2
C0[N-1] = np.zeros(N+1)
C0[N] = np.zeros(N+1)
C1[N-1] = np.zeros(N+1)
C1[N] = np.zeros(N+1)
C2[N-1] = np.zeros(N+1)
C2[N] = np.zeros(N+1)
#contrucción de la matriz Lambda
I = np.identity(N+1)
ceros = np.zeros([N+1,N+1])
A = np.block([[-C0i@C1,-C0i@C2,-C0i@C3],[I,ceros,ceros],[ceros,I,ceros]])
B = np.block([[C0,ceros,ceros],[ceros,I,ceros],[ceros,ceros,I]])

## Resolviendo la ecuación matricial
# eigenvalues = spla.eigvals(C, B*1j, check_finite=Fals
eigenvalues, eigenvectors = la.eig(A,B,check_finite=False)
#alpha_i_complex = np.amax(np.extract(eigenvalues != np.inf, eigenvalues))
#alpha_i = np.real(alpha_i_complex)

eigenv = 1/eigenvalues
eigenv1 = np.zeros(3*(N+1),dtype='complex_')
for i in np.arange(3*(N+1)):
    if np.real(eigenvalues[i])>0:
        if np.real(eigenvalues[i])>1:
            if np.imag(eigenvalues[i])>0:
                eigenv1[i] = eigenv[i]
Alpha = np.extract(np.imag(eigenv1) != 0 , eigenv1)
alphamax = np.max(Alpha)
print(alphamax)
#print(np.shape(m))
#print(np.linalg.det(rn))
#detC0 = np.linalg.det(rn)
#invR = np.linalg.inv(C0)
#print(eigenvalues)


'''
phi_r_list = np.real(eigenvectors[:,1]).tolist()
phi_r_tuple = tuple(phi_r_list)
phi_r_cheby = npcheby.Chebyshev(phi_r_tuple, [ymin,ymax])

#phi_i_list = np.imag(eigenvectors[:,1]).tolist()
#phi_i_tuple = tuple(phi_i_list)
#phi_i_cheby = npcheby.Chebyshev(phi_i_tuple, [ymin,ymax])# - np.imag(c_i_complex)
#x = np.linspace(-1,1,101)
#print(c_i)
print(np.real(eigenvalues))

print("")
print("N:", N)
print("alpha:", alpha)
print("c_i:",c_i.round(4))
print("c_i complex:",c_i_complex.round(4))
print("ω_i:",(alpha*c_i).round(5))
print("ω_i complex:",(alpha*c_i_complex).round(5))

plt.figure(1, dpi=200)
plt.plot(y,phi_r_cheby(y))
plt.title("$\phi_r(y)$") # plt.title(str(N))
plt.xlabel("$y$"); plt.ylabel("$\phi_r(y)$")
plt.grid(linestyle='--')
plt.show()
'''