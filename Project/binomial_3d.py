import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, multivariate_normal


def call_put_flag(call_put='call'):
    """ This method is used as a call or put flag.
        For call option, return 1.
        For put option, return -1.
    """

    if call_put == "call":
        return 1  # call option
    else:
        return -1  # put option


def payoff(s1, s2, num1, num2, K, type=1, call_put='call'):
    """ This function will return Payoffs of rainbow options. """

    if call_put == "call":  # z is the call put flag
        z = 1
    else:
        z = -1

    if type == 1:  # Option on the maximum of two assets
        return max(0, z * max(num1 * s1, num2 * s2) - z * K)
    elif type == 2:  # Option on the minimum of two assets
        return max(0, z * min(num1 * s1, num2 * s2) - z * K)
    elif type == 3:  # Spread Option
        return max(0, z * (num1 * s1 - num2 * s2) - z * K)


def binomial(s1, s2, num1, num2, b1, b2, q1, q2, sig1, sig2, rho, expT, n, r, K, type=1, call_put='call',
             euro_amer="e"):
    """ This method will return option prices."""

    dt = expT / n
    mu1 = b1 - (sig1 ** 2) / 2
    mu2 = b2 - (sig2 ** 2) / 2
    u = np.exp(mu1 * dt + sig1 * np.sqrt(dt))
    d = np.exp(mu1 * dt - sig1 * np.sqrt(dt))

    OptionValue = np.zeros(n ** 2).reshape(n, n)

    for j in range(0, n):
        Y1 = (2 * j - n) * np.sqrt(dt)
        NodesValueS1 = s1 * u ** j * d ** (n - j)

        for i in range(0, n):
            NodesValueS2 = s2 * np.exp(mu2 * n * dt) * np.exp(
                sig2 * (rho * Y1 + np.sqrt(1 - rho ** 2) * (2 * i - n) * np.sqrt(dt)))

            OptionValue[j][i] = payoff(NodesValueS1, NodesValueS2, num1, num2, K, type, call_put)

    for m in range(n - 1, 0, -1):
        for j in range(0, m):
            Y1 = (2 * j - m) * np.sqrt(dt)
            NodesValueS1 = s1 * u ** j * d ** (m - j)

            for i in range(0, m):
                Y2 = rho * Y1 + np.sqrt(1 - rho ** 2) * (2 * i - m) * np.sqrt(dt)
                NodesValueS2 = s2 * np.exp(mu2 * m * dt) * np.exp(sig2 * Y2)
                OptionValue[j][i] = 0.25 * (
                        OptionValue[j][i] + OptionValue[j][i + 1] + OptionValue[j + 1][i] + OptionValue[j + 1][
                    i + 1]) * np.exp(-r * dt)

                if euro_amer == "a":
                    OptionValue[j][i] = max(OptionValue[j][i],
                                            payoff(NodesValueS1, NodesValueS2, num1, num2, K, type=1, call_put="call"))
    return OptionValue[0][0]


if __name__ == '__main__':
    s1 = 100
    s2 = 100

    num1 = 1
    num2 = 1

    q1 = 0.22
    q2 = 0.2
    sig1 = 0.2657
    sig2 = 0.2578

    rho = 0.895998
    expT = 1
    n = 100
    r = 0.0237
    K = 156.3
    type = 1
    call_put = 'call'

    b1 = r - q1
    b2 = r - q2

    print(binomial(s1, s2, num1, num2, b1, b2, q1, q2, sig1, sig2, rho, expT, n, r, K, type=1, call_put='call',
                   euro_amer="e"))
    print(binomial(s1, s2, num1, num2, b1, b2, q1, q2, sig1, sig2, rho, expT, n, r, K, type=1, call_put='call',
                   euro_amer="a"))

    euro = []
    amer = []
    rhos = np.linspace(-1, 1, 50)
    for rho in rhos:
        e = binomial(s1, s2, num1, num2, b1, b2, q1, q2, sig1, sig2, rho, expT, n, r, K, type=1, call_put='call',
                     euro_amer="e")
        a = binomial(s1, s2, num1, num2, b1, b2, q1, q2, sig1, sig2, rho, expT, n, r, K, type=1, call_put='call',
                     euro_amer="a")
        euro.append(e)
        amer.append(a)

    plt.plot(rhos, euro, label="European")
    plt.plot(rhos, amer, label="American")
    plt.title("Option Prices vs Correlation")
    plt.xlabel("Rho")
    plt.ylabel("Option Prices")
    plt.legend()

    plt.show()
