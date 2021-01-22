import numpy as np
import pandas as pd
from scipy.linalg import toeplitz
from scipy.stats import norm



class NumPDE:

    def __init__(self,s0,Ns,Nt,smin,smax,sig,r,tt,k1,k2,contigent,option_type='European'):
        self.s0 = s0
        self.smin = smin
        self.smax = smax
        self.Ns = Ns
        self.Nt = Nt
        self.sig = sig
        self.r = r
        self.tt = tt
        self.type = option_type
        self.k1 = k1
        self.k2 = k2
        self.contigent = contigent


        self.ht = tt / Nt
        self.hs = (smax - smin) / Ns

        self.s = np.linspace(self.smin,self.smax,self.Ns+1)

    def matrix_a(self):

        # explicit euler discretizing of matrix A

        self.a = 1 - (self.sig ** 2) * (self.s ** 2) * self.ht / (self.hs ** 2) - self.r * self.ht
        self.l = 0.5 * (self.sig ** 2) * (self.s ** 2) * self.ht / (self.hs ** 2) - self.r * self.s * self.ht / (2 * self.hs)
        self.u = 0.5 * (self.sig ** 2) * (self.s ** 2) * self.ht / (self.hs ** 2) + self.r * self.s * self.ht / (2 * self.hs)

        A = np.diag(self.a[1:self.Ns])  # main diagonal

        # replace with u
        upper_diag = self.u[1:self.Ns-1]
        for i in range(len(upper_diag)):
            A[i+1][i] = upper_diag[i]

        # replace with l
        lower_diag = self.l[2:self.Ns]
        for i in range(len(lower_diag)):
            A[i][i+1] = lower_diag[i]

        return A



    def eigen(self):
        # compute eigen value of A
        eigen = np.linalg.eigvals(self.matrix_a())
        return abs(eigen)


    def price_call_spread(self):
        A = self.matrix_a()

        long = self.s - self.k1
        long[long < 0] = 0
        short = self.s - self.k2
        short[short < 0] = 0
        payoff = long - short
        price = payoff[1:Ns]


        for i in range(Nt):
            price = A.dot(price)
            price[0] = price[0] + self.l[1] * (self.k2 - self.k1) * np.exp(-self.r * i * self.ht)

            if self.type == 'American':
                for i in range(len(price)):
                    if price[i] < payoff[i+1]:      # only use payoff from 1:(Ns-1)
                        price[i] = payoff[i+1]

        return np.interp(s0, self.s[1:Ns], price)


    def price_put_spread(self):
        A = self.matrix_a()

        long =  self.k1 - self.s
        long[long < 0] = 0
        short =  self.k2 - self.s
        short[short < 0] = 0
        payoff = long - short
        price = payoff[1:Ns]


        for i in range(Nt):
            price = A.dot(price)
            price[0] = price[0] + self.l[1] * (self.k2 - self.k1) * np.exp(-self.r * i * self.ht)

            if self.type == 'American':
                for i in range(len(price)):
                    if price[i] < payoff[i+1]:
                        price[i] = payoff[i+1]

        return np.interp(s0, self.s[1:Ns], price)


    def price_call(self):

        A = self.matrix_a()

        k = self.k1
        payoff = self.s - self.k
        payoff[payoff < 0] = 0
        price = payoff[1:Ns]

        # solving bs pde using explicit finite difference method
        for i in range(Nt):
            price = A.dot(price)
            price[0] = price[0] + self.l[1] * (self.k2 - self.k1) * np.exp(-self.r * i * self.ht)

            if self.type == 'American':
                for i in range(len(price)):
                    if price[i] < payoff[i+1]:
                        price[i] = payoff[i+1]

        return np.interp(s0, self.s[1:Ns], price)



def BS_formula(s0,k,r,sig,tt,contigent):
    d1 = (np.log(s0 / k) + (r + 0.5 * sig ** 2) * tt) / (sig * np.sqrt(tt))
    d2 = d1 - sig * np.sqrt(tt)

    if contigent == "call":
        return s0 * norm.cdf(d1) - k * np.exp(-r * tt) * norm.cdf(d2)

    elif contigent == "put":
        return k * np.exp(-r * tt) * norm.cdf(-d2) - s0 * norm.cdf(-d1)



if __name__ == "__main__":

    s0 = 277.33

    r = 0.0247

    smin = 0
    smax = 500
    Ns = 250
    tt = 7/12
    Nt = 1000
    sig = 0.14    # VIX data
    k1 = 285
    k2 = 290


    contigent = 'call_spread'
    numpde_call_spread_euro = NumPDE(s0,Ns,Nt,smin,smax,sig,r,tt,k1,k2,contigent,option_type='European')
    a = numpde_call_spread_euro.matrix_a()
    print(max(np.linalg.eigvals(a)))

    print("The price of call spread without early exercise is  " ,(numpde_call_spread_euro.price_call_spread()))

    contigent = 'call_spread'
    numpde_call_spread_amer = NumPDE(s0, Ns, Nt, smin, smax, sig, r, tt, k1, k2, contigent, option_type='American')
    print("The price of call spread with early exercise is: " ,(numpde_call_spread_amer.price_call_spread()))











