from math import exp,log,sqrt,pi
from scipy.stats import norm,multivariate_normal
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



def vinilla(s0,k,r,sig1,t1,LB,UB,nodes):

    x,w = np.polynomial.legendre.leggauss(nodes)

    intervals = 0.5*(UB-LB)*x + (UB + LB)/2
    payoff = intervals - k
    integrand = payoff * w * norm.pdf(intervals,loc = s0,scale = sig1)

    E_integral = sum(integrand) * 0.5 * (UB - LB) / (1+r)
    return E_integral



def contigent_payoff(s2,s1,k1,k2):
    if s2 < k2:
        return max(s1 - k1,0)
    else:
        return 0


def density(x,y,sig1,sig2,mean1,mean2,pho):
    a1 = 1/(2 * pi * sig1 * sig2 * sqrt(1-pho**2))
    a2 = (-1/(2 * (1- pho ** 2))) * ( (x-mean1)**2 /sig1**2  +  (y-mean2)**2/sig2**2  -  2*pho*(x-mean1)*(y-mean2)/(sig1*sig2))
    return a1 * exp(a2)

def bi_normal_pdf(S1, S2,sig1,sig2,S0,rho):
    coef = 1 / (2 * pi * sig2 * sig1 * (1 - rho ** 2) ** 0.5)
    part1 = (S1 - S0) ** 2 / sig1 ** 2
    part2 = (S2 - S0) ** 2 / sig2 ** 2
    part3 = 2 * rho * (S1 - S0) * (S2 - S0) / sig2 / sig1
    core = np.exp(-0.5 / (1 - rho ** 2) * (part1 + part2 - part3))

    return coef * core



def contigent(s0,r,sig1,sig2,pho,k1,k2):

    LB1 = s0 + 5 * sig1
    UB1 = s0 - 5 * sig1
    LB2 = s0 + 5 * sig2
    UB2 = s0 - 5 * sig2

    var = multivariate_normal(mean=[s0, s0], cov=[[sig1 ** 2, pho * sig1 * sig2], [pho * sig1 * sig2, sig2 ** 2]])

    interval1 = np.linspace(LB1, UB1, 1000)      # outer
    interval2 = np.linspace(LB2, UB2, 1000)      # inner

    delta1 = interval1[1] - interval1[0]
    delta2 = interval2[1] - interval2[0]

    # outer_integrand = []
    integral = 0
    for i in range(len(interval1)-1):

        for j in range(len(interval2)-1):
            # p_x = var.pdf([x1,x2])
            payoff = contigent_payoff(interval2[j], interval1[i], k1, k2)

            integral +=  payoff * bi_normal_pdf(interval1[i], interval2[j],sig1,sig2,s0,pho) * (interval1[i+1] - interval1[i]) * (interval2[j+1] - interval2[j])
            # mid_price.append(p_x * payoff * delta1 * delta2)
            # p_x * delta1 * (x1 - k) * delta2
        #outer_integrand.append(sum(mid_price * delta2))

    # E_integral = sum(mid_price)
    return  integral





if __name__ == "__main__":

    #### parameters
    s0 = 271
    k1 = 260
    k2 = 250

    sig1 = 20
    sig2 = 15
    r = 0
    t1 = 12/12
    t2 = 6/12
    pho = 0.95

    option_type = "call"
    LB1 = k1
    UB1 = s0 + 4 * sig1
    nodes = 5

    print()
    print("1. The price of the vinilla option is: %.5f " %vinilla(s0, k1, r, sig1, t1, LB1, UB1, nodes))
    print()

    # 2
    print("2. The price of contingent option is:  %.8f" %contigent(s0, r, sig1, sig2, pho, k1, k2))
    print()

    # 3
    phos = [0.8,0.5,0.2]
    for pho in phos:
        print("3. The price of contingent option with rho = "+ str(pho) +" is: %.8f " % contigent(s0, r, sig1, sig2, pho, k1, k2))
    print()


    # 5
    prices = [240,230,220]
    for k2 in prices:
        print("5. The price of contingent option with 6 months SPY = "+ str(k2) +" is: %.8f " % contigent(s0, r, sig1, sig2, pho, k1, k2))
    print()