from math import exp,log,sqrt
from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



def BS(s0,k,t,r,sig,option_type):

    d1 = (log(s0 / k) + (r + 0.5 * sig ** 2) * t) / (sig * sqrt(t))
    d2 = d1 - sig * sqrt(t)

    if option_type.lower() == 'call':
        ec = s0 * norm.cdf(d1) - k * exp(-r * t) * norm.cdf(d2)
        return ec
    if option_type.lower() == 'put':
        ep = k * exp(-r * t) * norm.cdf(-d2) - s0 * norm.cdf(d1)
        return ep



def left_Riemann(s0,k,t,r,sig,nodes,LB,UB):


    #log_mean =  (r - 0.5 * (sig **2)) * t
    # log_sd = sig * (sqrt(t))
    e = exp(1)
    delta = (UB - LB)/nodes

    integral = 0
    intervals = np.linspace(LB,UB-delta,nodes)
    payoff = s0 * e ** ((r - 0.5 * sig ** 2) * t + (sig * t ** 0.5 * intervals)) - k


    '''
    for i in range(len(intervals)):
        e = exp(1)
        payoff = s0 * e ** ((r - 0.5 * sig ** 2) * t + (sig * t ** 0.5 * intervals[i])) - k

        integrand = payoff * norm.pdf(intervals[i]) * delta   # payoff * pdf(x) * dx
        # integrand = payoff * (norm.cdf(intervals[i+1]) - norm.cdf(intervals[i]))

        integral += integrand
    '''

    integrand = payoff * norm.pdf(intervals) * delta
    E_integral = sum(integrand) * exp(-r * t)

    return E_integral


def midpoint(s0,k,t,r,sig,nodes,LB,UB):
    d1 = (log(s0 / k) + (r + 0.5 * sig ** 2) * t) / (sig * sqrt(t))
    d2 = d1 - sig * sqrt(t)

    e = exp(1)

    # log_mean =  (r - 0.5 * (sig **2)) * t
    # log_sd = sig * (sqrt(t))

    delta = (UB - LB) / nodes

    integral = 0
    intervals = np.linspace(LB + 0.5 * delta, UB - 0.5 * delta, nodes)

    payoff = s0 * e ** ((r - 0.5 * sig ** 2) * t + (sig * t ** 0.5 * intervals)) - k
    integrand = payoff * norm.pdf(intervals) * delta
    E_integral = sum(integrand) * exp(-r * t)

    return E_integral



def guass_nodes(s0,k,t,r,sig,nodes,LB,UB):
    e = exp(1)
    x, w = np.polynomial.legendre.leggauss(nodes)

    intervals = 0.5 * (x + 1) * (UB - LB) + LB
    payoff = s0 * e ** ((r - 0.5 * sig ** 2) * t + (sig * t ** 0.5 * intervals)) - k
    integrand = payoff * w * norm.pdf(intervals)
    E_integral = sum(integrand) * 0.5 * (UB - LB) * exp(-r * t)


    return E_integral



if __name__ == "__main__":
    s0 = 10
    k = 12
    sig = 0.2
    r = 0.04
    t = 3 / 12
    option_type = "call"


    # 1
    d1 = (log(s0 / k) + (r + 0.5 * sig ** 2) * t) / (sig * sqrt(t))
    d2 = d1 - sig * sqrt(t)



    LB = -d2
    UB = 4

    c1 = BS(s0, k, t, r, sig, option_type)
    print("The price of the European " + option_type + " option from BS formula is: %.5f" %c1)
    print()




    # 2
    N = [5, 10, 50, 100]
    tbl = pd.DataFrame(0.0000000, index=["N=" + str(x) for x in N], columns=["Left", "Mid", "Gauss"])




    for nodes in N:
        c_lf_riemann = left_Riemann(s0,k,t,r,sig,nodes,LB,UB)
        print("When node is "+str(nodes))
        print("The price of the European " + option_type + " option using Left Riemann Sum is: %.5f" %c_lf_riemann)
        print("error is: %.5f" %(c1 - c_lf_riemann))
        tbl["Left"][str(nodes)] = c_lf_riemann

        c_midpoint = midpoint(s0,k,t,r,sig,nodes,LB,UB)
        print("The price of the European " + option_type + " option using Midpoint is: %.5f" %c_midpoint)
        print("error is: %.5f" %(c1 - c_midpoint))
        tbl["Mid"][str(nodes)] = c_midpoint


        c_gauss = guass_nodes(s0,k,t,r,sig,nodes,LB,UB)
        print("The price of the European " + option_type + " option using Gauss Nodes is: %.5f" %c_gauss)
        print("error is: %.5f" %(c1 -  c_gauss))
        print()
        tbl["Gauss"][str(nodes)] = c_gauss



    # 3
    err1 = []
    err2 = []
    err3 = []
    for i in np.linspace(2,20,10):
        err1.append(abs(c1 - left_Riemann(s0,k,t,r,sig,int(i),LB,UB)))
        err2.append(abs(c1 - midpoint(s0,k,t,r,sig,int(i),LB,UB)))
        err3.append(abs(c1 - guass_nodes(s0,k,t,r,sig,int(i),LB,UB)))
    plt.plot(err1)
    plt.plot(err2)
    plt.plot(err3)

    plt.ylabel("Error")
    plt.show()

