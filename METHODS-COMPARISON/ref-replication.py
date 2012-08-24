from numpy import exp as nexp
from numpy import array
from math import log, exp, sqrt, pi
from scipy.optimize import fsolve
from operator import mul

#
# WU
#

def get_npsi(s):
    xs = [0.24534,
          0.73747,
          1.23407,
          1.73853,
          2.25497,
          2.78880,
          3.34785,
          3.94476,
          4.60369,
          5.38748]
    ws = [4.62243 * 10**(-1),
          2.86675 * 10**(-1),
          1.09017 * 10**(-1),
          2.48105 * 10**(-2),
          3.24377 * 10**(-3),
          2.28338 * 10**(-4),
          7.80255 * 10**(-6),
          1.08606 * 10**(-7),
          4.39934 * 10**(-10),
          2.22939 * 10**(-13)]
    xs = array(xs + map(lambda x: -x, xs))
    ws = array(2 * ws)

    def _f(mu, sigsq):
        return sum(list( ws / sqrt(pi) * nexp(-s * nexp(sqrt(2*sigsq) * xs + mu))))

    return _f


def get_nfit(mus, sigsqs, s1, s2):
    mgfi1 = get_npsi(s1)
    mgfi2 = get_npsi(s2)

    mgfi1s = reduce(mul, [mgfi1(mu1, sigsq1) for mu1, sigsq1 in zip(mus, sigsqs)])
    mgfi2s = reduce(mul, [mgfi2(mu2, sigsq2) for mu2, sigsq2 in zip(mus, sigsqs)])

    def _f(x):
        mu, sigsq = x
        return [mgfi1(mu, sigsq) - mgfi1s, mgfi2(mu, sigsq) - mgfi2s]
    return _f

def get_OR_WU(mus,sigsqs, s1, s2, initv=[5,1], maxit=1000):
    f = get_nfit(mus, sigsqs, s1, s2)
    return fsolve(f, initv, maxfev=maxit)


print get_OR_WU((8, 4), (2,1), 0.01, 0.05)
# res = [ 7.64510201  1.13001788]

print get_OR_WU( (2,2), (1,1), 0.01, 0.05 )
# res = [ 2.8993171   0.57862103]

#
# SY
#

AG1 = [-0.1153,
       -0.5667,
       0.1151,
       -0.71624,
       0.13129,
       -0.16113,
       0.20842,
       -0.44997,
       0.27562,
       -0.51097,
       0.13451,
       -0.26701,
       0.59191,
       -0.36953,
       0.69127,
       0.80570,
       0.14297,
       -0.32259,
       0.20556,
       -0.38886,
       -0.31453,
       -0.27300,
       0.62442,
       -0.40482,
       0.77467]

AG1P = [0,1,2,1,1,0,2,2,2,1,0,2,2,2,1,-1,2,2,2,1,-1,1,1,1,0]

AG2 = [0.40128,
       -0.44832,
       0.73917,
       -0.37721,
       0.52622,
       -0.15791,
       0.16372,
       -0.29949,
       0.13378,
       -0.17474,
       0.21745,
       -0.21141,
       0.39767,
       -0.17768,
       0.21953,
       -0.75253,
       0.11471,
       -0.21895,
       0.10067,
       -0.12249,
       0.84479,
       -0.22228,
       0.42895,
       -0.20357,
       0.24995]

AG2P = [-1,1,1,1,0,1,2,2,2,1,1,2,2,2,1,0,2,2,2,1,-1,1,1,1,0]

AG3 = [-0.39586,
       -0.54549,
       0.11392,
       -0.71163,
       0.13152,
       0.68399,
       0.19625,
       -0.43140,
       0.26477,
       -0.49405,
       -0.47172,
       -0.24868,
       0.56185,
       -0.35018,
       0.65609,
       0.18193,
       0.13235,
       -0.30500,
       0.19377,
       -0.36583,
       -0.28174,
       -0.25182,
       0.58933,
       -0.38094,
       0.72601]

AG3P = [1,1,2,1,1,1,2,2,2,1,1,2,2,2,1,1,2,2,2,1,0,1,1,1,0]

def buildCoeffDict():
    res = {'G1':[], 'G2':[], 'G3':[]}
    for i in range(0,5):
        tmpg1, tmpg2, tmpg3 = [], [], []
        for j in range(0,5):
            tmpg1.append(AG1[5*i+j] * 10**(AG1P[5*i + j]))
            tmpg2.append(AG2[5*i+j] * 10**(AG2P[5*i + j]))
            tmpg3.append(AG3[5*i+j] * 10**(AG3P[5*i + j]))
        res['G1'].append(tmpg1)
        res['G2'].append(tmpg2)
        res['G3'].append(tmpg3)
    return res

Aij = buildCoeffDict()

def G1(m, ssq):
    return 10 ** sum( [ Aij['G1'][j][k] * sqrt(ssq)**(j/float(2)) * abs(m) ** (k/float(2)) for j in range(0,5) for k in range(0,5)] )

def G2(m, ssq):
    return 10 ** sum( [ Aij['G2'][j][k] * sqrt(ssq)**(j/float(2)) * abs(m) ** (k/float(2)) for j in range(0,5) for k in range(0,5)] )

def G3(m, ssq):
    return 10 ** sum( [ Aij['G3'][j][k] * sqrt(ssq)**(j/float(2)) * abs(m) ** (k/float(2)) for j in range(0,5) for k in range(0,5)] )

def _get_OR_SY(mus, sigsqs, nterm=40):
    mu1, mu2 = map(float, mus)
    ssq1, ssq2 = map(float, sigsqs)
    if mu1 < mu2:
        mu2, mu1 = mu1, mu2
        ssq1, ssq2 = ssq2, ssq1
    mw = mu2 - mu1
    ssqw = ssq1 + ssq2
    #print G1(mw, ssqw)
    #print G3(mw,ssqw)
    #print G2(mw,ssqw)
    VZ = ssq1 - G1(mw, ssqw) ** 2 - 2 * ssq1 / ssqw * G3(mw, ssqw) + G2(mw, ssqw)
    EZ = mu1 + G1(mw, ssqw)
    return (EZ, VZ)


print(_get_OR_SY([2,2], [4,4]))

#print(sum(AG1))
#print(sum(AG1P))
#print(sum(AG2))
#print(sum(AG2P))
#print(sum(AG3))
#print(sum(AG3P))
#print C(5)
#print b(5)
#print b(3)

#
# FW
#

def get_u1(mus, sigmasqs):
    return sum( list( nexp( array(mus) + array(sigmasqs) / 2 ) ) )

def get_u2(mus, sigmasqs):
    mus = mus
    sigmasqs = sigmasqs
    t1 = sum( list( nexp(2* array(mus) + 2* array(sigmasqs) ) ) )
    t2 = 0
    for mu1, sigsq1 in zip(mus[:-1], sigmasqs[:-1]):
        for mu2, sigsq2 in zip(mus[1:], sigmasqs[1:]):
            t2 += 2 * exp(mu1 + mu2) * exp( 0.5 * (sigsq1 + sigsq2) )
    return(t1 + t2)

def get_FW(mus, sigmasqs):
    u1 = get_u1(mus, sigmasqs)
    print u1
    u2 = get_u2(mus, sigmasqs)
    print u2
    return ( 2 * log(u1) - 0.5 * log(u2), log(u2) - 2 * log(u1))

#print(get_FW([2,2], [4,4]))
