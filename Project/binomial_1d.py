import matplotlib.pyplot as plt
import numpy as np


def Binomial(n, S, K, r,div, v, t, PutCall):
    At = t / n
    u = np.exp((r-div-0.5*v**2) * At +v * np.sqrt(At))
    d = np.exp((r-div-0.5*v**2) * At -v * np.sqrt(At))
    p = 0.5

    # Binomial price tree
    stockvalue = np.zeros((n + 1, n + 1))
    stockvalue[0, 0] = S
    for i in range(1, n + 1):
        stockvalue[i, 0] = stockvalue[i - 1, 0] * u
        for j in range(1, i + 1):
            stockvalue[i, j] = stockvalue[i - 1, j - 1] * d

    # option value at final node
    optionvalue = np.zeros((n + 1, n + 1))
    for j in range(n + 1):
        if PutCall == "C":  # Call
            optionvalue[n, j] = max(0, stockvalue[n, j] - K)
        elif PutCall == "P":  # Put
            optionvalue[n, j] = max(0, K - stockvalue[n, j])

    # backward calculation for option price
    for i in range(n - 1, -1, -1):
        for j in range(i + 1):
            if PutCall == "P":
                optionvalue[i, j] = max(0, K - stockvalue[i, j], np.exp(-r * At) * (
                            p * optionvalue[i + 1, j] + (1 - p) * optionvalue[i + 1, j + 1]))
            elif PutCall == "C":
                optionvalue[i, j] = max(0, stockvalue[i, j] - K, np.exp(-r * At) * (
                            p * optionvalue[i + 1, j] + (1 - p) * optionvalue[i + 1, j + 1]))
    return optionvalue[0, 0]

    # Inputs


n = 100  # input("Enter number of binomial steps: ")           #number of steps
S = 100  # input("Enter the initial underlying asset price: ") #initial underlying asset price
r = 0.015  # input("Enter the risk-free interest rate: ")        #risk-free interest rate
div = 0.0177
K = 100  # input("Enter the option strike price: ")            #strike price
v = 0.35 # input("Enter the volatility factor: ")              #volatility
t = 1.

# Graphs and results for the Option prices

y = [-Binomial(n, S, K, r, div,v, t, "C")] * (K)
y += [x - Binomial(n, S, K, r,div, v, t, "C") for x in range(K)]

print("American Call Price: %s" % (Binomial(n, S, K, r,div, v, t, PutCall="C")))
