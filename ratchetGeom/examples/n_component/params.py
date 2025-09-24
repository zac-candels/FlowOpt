import numpy as np
from sympy import *

# This script solves the alpha, Lambda and Theta matrices needed in the N component free energy in the following paper
#  Boyer, Franck & Minjeaud, Sebastian. (2014). Hierarchy of consistent n-component Cahn–Hilliard systems. Mathematical
#  Models and Methods in Applied Sciences. 24. 2885-2928.

# Number of immiscible fluid components.
N=6
assert N>1,"Number of components must be greater than 1"

# identity matrix
id = eye(N)

# Matrix of surface tensions - MUST SPECIFY THIS, symmetric, diagonal 0.
#s = MatrixSymbol('s', N,N).as_explicit()
s=0.005*ones(N,N)
#s = Matrix([[0,0.007,0.005],[0.007,0,0.005],[0.005,0.005,0]])
#s = Matrix([[0,0.005,0.004,0.006,0.006],[0,0,0.0055,0.0045,0.0045],[0,0,0,0.005,0.005],[0,0,0,0,0.005],[0,0,0,0,0]])
ssym = Matrix(N, N, lambda i, j: s[min(i,j),max(i,j)] if (i!=j) else 0)

# Vector of gamma parameter from [1]
gamma = MatrixSymbol('gamma', N,1).as_explicit()

# Kronecker product of gamma vector with vector of ones in equation 2.7.
gammaxone = tensorproduct(gamma, ones(N,1))
# Need to convert to N x N matrix.
gammaxone = Matrix(N, N, lambda i, j: gammaxone[i,0,j,0])

# Solution for alpha matrix (gamma vector still unknown) from equation 2.7.
alpha = (id+gammaxone)*ssym.inv()

# Solution to gamma vector from equation 2.7.
gammasol = solve(alpha*ones(N,1),gamma)

# Solution for alpha only in terms of surface tensions.
alphasol = Matrix(N, N, lambda i, j: alpha[i,j].subs(gammasol))

print("Alpha (equation 2.7):")
f = open("params.txt", "w")
for i in range(N):
    print([alphasol[i,j] for j in range(N)])
    f.write(str([alphasol[i,j] for j in range(N)]).strip('[]').replace(',','')+"\n")
print("")

# Capital gamma from equation 3.16.
def calcGamma(i,j,k,l):
    return round(-6*(alphasol[i,j]*(ssym[j,k]+ssym[j,l])
                     +alphasol[i,k]*(ssym[k,j]+ssym[k,l])
                     +alphasol[i,l]*(ssym[l,j]+ssym[l,k])-gammasol[gamma[i,0]]),15)*3

# This returns a list of Lambdas (equation 3.18) for a given triplet of components. These Lambdas are sorted from
# lowest to highest index.
def calcLambda(j,k,l):
    idcs = np.array([],dtype=int)
    for m in range(N):
        if (m!=j and m!=k and m!=l):
            idcs=np.append(idcs,m)
    asubmat = Matrix(N-3, N-3, lambda m, n: alphasol[idcs[m],idcs[n]])
    GammaMat = Matrix(N-3, 1, lambda m, n: calcGamma(idcs[m],j,k,l))
    
    return linsolve((asubmat,GammaMat))

# Fills a dictionary of Lambdas with the triplet indices 'jkl' as the key and the Lambdas of each index (sorted in 
# increasing order) not equal to i, j or k as the value e.g for N=5 and the triplet of components 1,2,3, the key is
# '123' and the value is a list of the Lambdas for components 4 and 5 in that order.
Lambdas={}
for j in range(N):
    for k in range(j+1,N):
        for l in range(k+1,N):
            idcs = symbols(str(j)+str(k)+str(l))
            Lambdas[idcs]=calcLambda(j,k,l)

print("Lambda (equation 3.18):")
for key in Lambdas:
    j=0
    for i in range(N):
        if (i!=int(str(key)[0]) and i!=int(str(key)[1]) and i!=int(str(key)[2])):
            print("Lambda_"+str(i)+"_"+str(key),Lambdas[key].args[0][j])
            f.write(str(i)+"_"+str(key)+' '+str(Lambdas[key].args[0][j])+'\n')
            j+=1
print("")

# This returns a list of Thetas (equation 3.19) for a given triplet of components in a particular order (only the last
# index l0 given matters for the order, see section 3.2.3). These Thetas are sorted from lowest to highest index.
def calcTheta(j,k,l):
    idcs = np.array([],dtype=int)
    for m in range(N):
        if (m!=j and m!=k and m!=l):
            idcs=np.append(idcs,m)
    asubmat_lhs = Matrix(N-3, N-3, lambda m, n: alphasol[idcs[m],idcs[n]])
    asubmat_rhs = Matrix(N-3, 1, lambda m, n: 2*alphasol[idcs[m],l])
    
    return linsolve((asubmat_lhs,asubmat_rhs))

# Fills a dictionary of Thetas with the triplet indices 'jk_l' as the key and the Thetas of each index (sorted in
# increasing order) not equal to i, j or k as the value e.g for N=5 and the triplet of components 1,2,3, with l=3, the
# key is '12_3' and the value is a list of the Lambdas for components 4 and 5 in that order
Thetas={}
for j in range(N):
    for k in range(j+1,N):
        for l in range(k+1,N):
            if (j!=k and j!=l and k!=l):
                idcs = symbols(str(j)+str(k)+'_'+str(l))
                Thetas[idcs]=calcTheta(j,k,l)
                idcs = symbols(str(j)+str(l)+'_'+str(k))
                Thetas[idcs]=calcTheta(j,l,k)
                idcs = symbols(str(k)+str(l)+'_'+str(j))
                Thetas[idcs]=calcTheta(k,l,j)

print("Theta (equation 3.19):")
for key in Thetas:
    j=0
    for i in range(N):
        if (i!=int(str(key)[0]) and i!=int(str(key)[1]) and i!=int(str(key)[3])):
            print("Theta_"+str(key)[0]+str(key)[1]+"_"+str(i)+"_"+str(key)[3],Thetas[key].args[0][j])
            f.write(str(key)[0]+str(key)[1]+"_"+str(i)+"_"+str(key)[3]+' '+str(Thetas[key].args[0][j])+'\n')
            j+=1
f.close()