# -*- coding: utf-8 -*-
"""

"""
#%%
import numpy as np
from scipy import interpolate

#%%
S0=100
sig=0.2
q1=0.0177
q2=0.0077
r=0.015 - 0.0177
rho=0.5
K=100
T=1

#%% Numerical PDE parameter
smax=200
M=100
hs=smax/M
N=500
ht=T/N

#%%
si=np.array([hs*i for i in range(1,M+1)])
tj=np.array([ht*j for j in range(1,N+1)])
a=np.array([1-sig**2*i**2*ht/hs**2-r*ht for i in si])
l=np.array([sig**2*i**2*ht/(2*hs**2)-r*i*ht/(2*hs) for i in si])
u=np.array([sig**2*i**2*ht/(2*hs**2)+r*i*ht/(2*hs) for i in si])
A=np.diag(a[:(M-1)])+np.diag(l[1:(M-1)],-1)+np.diag(u[:(M-2)],1)

#%%
def find_price(S,C0,si):
    f=interpolate.interp1d(si,C0)
    p=f(S)
    return p

#%%
C=np.zeros((M-1,N))
C[:,N-1]=np.array([max(i-K,0)for i in si])[:(M-1)]
for j in range(N-2,-1,-1):
    Cj=np.dot(A,C[:,j+1])+np.array(np.append(np.repeat(0,M-2),u[M-2]*(smax-K)))
    C[:,j]=np.maximum(Cj,[max(i-K,0) for i in si[:(M-1)]])    
print(find_price(S0,C[:,0],si[:-1]))
