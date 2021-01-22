import pandas as pd
import numpy as np
from sklearn.decomposition import PCA


# Problem 1

tickers = ['ATVI', 'ADBE', 'AMD', 'ALGN', 'ALXN', 'AMZN', 'AMGN', 'AAL', 'AAPL', 'AMAT', 'ASML', 'ADSK', 'ADP', 'AVGO', 'BIDU', 'BIIB', 'BMRN', 'CDNS', 'CELG', 'CERN', 'CHKP', 'CHTR', 'CTRP', 'CTAS', 'CSCO', 'CTXS', 'CMCSA', 'COST', 'CSX', 'CTSH', 'DLTR', 'EA', 'EBAY', 'EXPE', 'FAST', 'FB', 'FISV', 'FOX', 'FOXA', 'GILD', 'GOOG', 'GOOGL', 'HAS', 'HSIC', 'ILMN', 'INCY', 'INTC', 'INTU', 'ISRG', 'IDXX', 'JBHT', 'KLAC', 'LRCX', 'LBTYA', 'LBTYK', 'LULU', 'MELI', 'MAR', 'MCHP', 'MDLZ', 'MNST', 'MSFT', 'MU', 'MXIM', 'MYL', 'NTAP', 'NFLX', 'NTES', 'NVDA', 'NXPI', 'ORLY', 'PAYX', 'PCAR', 'BKNG', 'PEP', 'QCOM', 'REGN', 'ROST', 'SIRI', 'SWKS', 'SBUX', 'SYMC', 'SNPS', 'SPY', 'TTWO', 'TSLA', 'TXN', 'TMUS', 'ULTA', 'UAL', 'VRSN', 'VRSK', 'VRTX', 'WBA', 'WDC', 'WDAY', 'WLTW', 'WYNN', 'XEL', 'XLNX']

tf = pd.read_csv("./stock_data/" + tickers[0] + ".csv")
tf.index = tf['Date']
tf.fillna(method='ffill')
df = pd.DataFrame(index=tf.index)
df[tickers[0]] = tf['Adj Close']



for i in range(1,len(tickers)):
    tf = pd.read_csv("./stock_data/" + tickers[i] + ".csv")
    tf.index = tf['Date']
    tf.fillna(method='ffill')
    df[tickers[i]] = tf['Adj Close']


ret = pd.DataFrame(index=df.index)
ret.to_csv("./daily_return.csv")


for ticker in tickers:
    ret[ticker] = df[ticker] / df[ticker].shift(1) -1


# print(ret)
# compute the covaraince matrix
cov_mat = ret.cov()
cov_mat.to_csv("./covariance_matrix.csv")

w,v = np.linalg.eig(cov_mat)  # w: vector of eigenvalues, v:eigenvector




idx = w.argsort()[::-1]
w = w[idx]
v = v[:,idx]

# count how many positive eigenvalues in the matrix decomposition
num = [1 for i in range(len(w)) if w[i]>0]

explianed_var = [np.abs(i)/np.sum(w) for i in w]
sum_var = np.cumsum(explianed_var)
def find(sort_var,l):  #a is required variance
    i=0
    while sort_var[i]<l:
        i+=1
    return i

n1=find(sum_var,0.5)+1     # 8 of eigenvalues explianed 50% of the variance
n2=find(sum_var,0.9)+1     # 57 of eigenvalues explianed 90% of the variance


"""
pca = PCA(n_components=22)
pca.fit(cov_mat)
print(pca.explained_variance_ratio_)
print(sum(pca.explained_variance_ratio_))
"""

########### compute the vectors from svd   ##########################
svd1,svd_diag, svd2 = np.linalg.svd(cov_mat)
#print(svd1)
svd_diag_stable = list(svd_diag)[:n2]
svd_diag_recirocal = [1/i for i in svd_diag_stable]
# print(len(svd_diag_recirocal))
rest = np.zeros(len(tickers) - len(svd_diag_recirocal))
svd_all = np.append(svd_diag_recirocal,rest)
# print(svd_all)
diag_stable_cov = np.diag(svd_all)
# print(diag_stable_cov)



# Problem 2
##################################################################

# Parameters:
# Matrix:
# G: constriant matrix, size: N*2
# C: covariance matrix of returns, size: N*N
# R: expected return matrix of N companies
# c: the RHS of Constriants, size: 2*1
# a: parameter in a<w,Cw>, scalar

##################################################################

N = len(tickers)  # number of stocks

#  Form G matrix
g1 = np.ones(N)        # g1: budget constriant
g2_1 = np.ones(17)     # g2: constriant on the first 17 parameters
g2_2 = np.zeros(N-17)
G = np.matrix([g1,np.append(g2_1,g2_2)])
# print(G)


# compute the inverse of G(C.I)G.T
C_inverse = svd1 * diag_stable_cov * svd2
inverse = np.linalg.inv(G * C_inverse * G.T)
# print()
# print(inverse)




a = 1 # the parameter in a<w,Cw>

R = (np.matrix(np.mean(ret))).T   # R is the expected return matrix of the N companies, R size: N*1
c = np.matrix([[1],[0.1]])    # RHS of constriants

# compute the value of lambda
part1 = G * C_inverse * R - 2 * a * c
# print(part1)
# print()
lam = inverse * part1
# print(lam)


# compute weights of portfolio
weights = (1/(2*a)) * C_inverse * (R - G.T * lam)
# print(weights)

up = max(weights)
low = min(weights)
print(up)   # maximum weight
print(low)  # minimum weight
