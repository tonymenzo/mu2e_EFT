#
# Translation from Mathematica code in ../mathematica/mu2e_eft_6_2_2022.nb
# Ken McElvain   16June2022
#
# Update for lower components 28Nov2022
# Relativistic version
#
# Need to write Mathematica code to generate python test code.
#

# See https://martinheinz.dev/blog/50 for discussion of function overloading
from multipledispatch import dispatch
import os
import copy
import argparse
import pathlib
import pprint
import matplotlib.pyplot as plt
from collections import namedtuple # sort of a struct
import numpy as np
import math
import cmath
import py3nj  # angular momentum brackets
import scipy.special # hypergeometric and bessel functions
import yaml


# Set up arguments
parser = argparse.ArgumentParser(prog='mu2e', description='Analyse Mu->e transitions')
parser.add_argument("-test", action="store_true", help="Run internal tests")
parser.add_argument("-v", action="store_true", help="Enable more logging")
parser.add_argument("-pelements", action="store_true", help="Report element/isotope data")
elasticenvname = "MU2E_ELASTIC"
if elasticenvname in os.environ:
    elasticdirdflt = os.environ[elasticenvname]
else:
    elasticdirdflt = "./Elastic/"
parser.add_argument("-edir", action="store", type=str, default=elasticdirdflt,
    help="Directory with elastic scattering density matrices")
parser.add_argument('path', nargs='*', help='Paths to analysis input files')
args = parser.parse_args()

# We will use floating point numbers for half integer
# values.  The Mathematica equivalent uses rationals.
# This code relies on the exact representation of integer/2.0
# in floating point.

# some constants
log2 = math.log(2.0)
logsqrtpi = 0.5 * math.log(math.pi)

auMeV = 931.494028      #  (12C mass / 12) in MeV
mumass = 105.6583755    #  current value
alpha=1.0/137.035999082
emass = 0.511           # 0.5109989 MeV
hbar = 6.5821e-22       # 6.582119569e-22 MeV sec
hc = 197.3269718        # hbar * c (Mev fm)
Nmass=939.0             # Neutron Mass: 939.56542052 MeV
                        # Proton Mass:  938.272 088 16
mN=939.0                # Also used for nucleon mass
mv = 246200             # Weak scale in MeV, sqrt(1 / (sqrt(2) GF))
                        # Used to normalize EFT's LECs

# Supports testing code
def check(c, msg):
    global errcnt
    if not c:
        errcnt = errcnt + 1
    print(msg, " ", "success" if c else "fail")

# v is a value in 0, 0.5, 1, 1.5, 2, ...
# convert integers to strings directly
# half integers go to 2*v/2
def halfstr(v):
    if IntegerQ(v):
        return str(v)
    x = round(v * 2.0)
    return f"{x}/2"
#
# We pass half integer values around as floats
# This is because of half integer angular momentum quantum numbers
# Double and they become "integers"
#
def IntegerQ(v):
    x = abs(v - round(v))
    return x == 0.0

# Double factorial j!!
# I am restricting the domain to non-negative integers
def DFact(j):
    if j < 0 or not IntegerQ(j):
        raise ValueError(f"DFact: domain error, got {j}")
    j = int(j)
    if j & 1:  # odd
        # note:  gamma(n + 1) = n!
        k = (j + 1) >> 1 # j = 2 k - 1
        t = k * log2 + scipy.special.loggamma(k + 0.5) - logsqrtpi
        return math.exp(t)
    else:
        k = j >> 1
        return 2**k * math.factorial(k)

def testDFact():
    z = round(DFact(11))
    zc = 10395
    check(z == zc, f"DFact(11) = {z}, expecting {zc}")
    z = round(DFact(21))
    zc = 13749310575
    check(z == zc, f"DFact(21) = {z}, expecting {zc}")

    z = round(DFact(20))
    zc = 3715891200
    check(z == zc, f"DFact(20) = {z}, expecting {zc}")

    z = round(DFact(0))
    zc = 1
    check(z == zc, f"DFact(0) = {z}, expecting {zc}")

#
# Get number of z-projection states for angular momentum j
#
def JNorm(j):
    return 2 * j + 1

def QNorm(j):
    return math.sqrt(2 * j + 1)


#
# see if it is possible to couple angular momentum x, and y to form z
# The order of x, y, z doesn't matter
#
def Triangular(x,y,z):
    if z > x + y:
        return False
    if z < abs(x - y):
        return False
    return True

#
# check if z projection (m) is in allowed range for
# total angular momentum j
#
def RangeM(j, m):
    return abs(m) <= j

#
# Test if the three angular momentum states |j,m> can generate
# a non-zero result from 3J symbol/bracket
#
def ThreeJCondition(a, b, c):
    j1,m1 = a
    j2,m2 = b
    j3,m3 = c
    return (Triangular(j1, j2, j3) and
           RangeM(j1, m1) and
           RangeM(j2, m2) and
           RangeM(j3, m3) and
           (m1 + m2 + m3 == 0))

#
# Implement 3J symbol
#   / j1  j2  j3 \      (-1)^(j1-j2-m3)
#  |              | =   --------------  < j1,m1,  j2,m2 | j3,-m3>
#   \ m1  m2  m3 /      sqrt(2 j3 + 1)
#
def ThreeJ(a, b, c):
    if not ThreeJCondition(a, b, c):
        return 0.0
    j1,m1 = a
    j2,m2 = b
    j3,m3 = c
    j1 = round(j1 * 2.0)
    j2 = round(j2 * 2.0)
    j3 = round(j3 * 2.0)
    m1 = round(m1 * 2.0)
    m2 = round(m2 * 2.0)
    m3 = round(m3 * 2.0)
    return py3nj.wigner3j(j1, j2, j3, m1, m2, m3)

#
# Used to check triples (internal angular momentum sums) of angular momentum values from 6J
# for consistancy.   Basically the Triangle rule and sum to integer.
#
def SixJConditionTriad(j1, j2, j3):
    jsum = j1 + j2 + j3
    if not IntegerQ(jsum):
        return False
    return Triangular(j1, j2, j3)

#
# Check all triples
#
def SixJCondition(j1, j2, j3, J1, J2, J3):
    return (SixJConditionTriad(j1, j2, j3) and
           SixJConditionTriad(j1, J2, J3) and
           SixJConditionTriad(J1, j2, J3) and
           SixJConditionTriad(J1, J2, j3))

#
# 6J symbol
# Connected to Racah W-coefficients which express recoupling of 3 angular momenta
#
def SixJ(ju, jl):
    j1,j2,j3 = ju
    J1,J2,J3 = jl
    if not SixJCondition(j1, j2, j3, J1, J2, J3):
        return 0.0
    j1 = round(j1 * 2.0)
    j2 = round(j2 * 2.0)
    j3 = round(j3 * 2.0)
    J1 = round(J1 * 2.0)
    J2 = round(J2 * 2.0)
    J3 = round(J3 * 2.0)
    rslt = py3nj.wigner6j(j1, j2, j3, J1, J2, J3)
    return rslt

#
# 9J symbol
# Connected to Racha W-coefficients which express recoupling of 4 angular momenta
#
def NineJSymbol(a, b, c):
    j1,j2,j12=a
    j3,j4,j34=b
    j13,j24,j=c
    # convert to doubled integer values
    j1 = round(j1 * 2.0)
    j2 = round(j2 * 2.0)
    j12 = round(j12 * 2.0)
    j3 = round(j3 * 2.0)
    j4 = round(j4 * 2.0)
    j34 = round(j34 * 2.0)
    j13 = round(j13 * 2.0)
    j24 = round(j24 * 2.0)
    j = round(j * 2.0)
    return py3nj.wigner9j(j1,j2,j12, j3,j4,j34, j13,j24,j)

#
# Add a consistant name for 9J with the name for 6J and 3J
# that cleanly returns 0 for out of bounds values
#
def NineJ(a, b, c):
    return NineJSymbol(a, b, c)


########################################
# Ancillary Functions:   Bessel Function
########################################
def BesselFactor1(y, oscp, osc, lcap):
    np,lp=oscp
    n,l=osc
    r = 2 ** lcap
    r /= DFact(JNorm(lcap))
    r *=  y**(lcap/2.0)
    r *= math.exp(- y)
    r *= math.sqrt( math.factorial(np - 1) * math.factorial(n - 1) )
    return r

def testBesselFactor1():
    z = BesselFactor1(0.625, (2, 1), (1, 2), 2)
    zc = 0.08921023808649838
    check(abs(z-zc) < 1e-8, f"BesselFactor1(0.625, [2,1], [1,2], 2) = {z}, expecting {zc}")

def BesselFactor2(y, oscp, osc, lcap):
    np,lp=oscp
    n,l=osc
    r = math.sqrt(scipy.special.gamma(np + lp + 0.5) * scipy.special.gamma(n + l + 0.5))
    return r

def testBesselFactor2():
    z = BesselFactor2(1.625, (2, 1), (1, 2), 2)
    zc = 15.0 * math.sqrt(math.pi) / 8.0
    check(abs(z-zc) < 1e-8, f"BesselFactor2(1.625, [2,1], [1,2], 2) = {z}, expecting {zc}")

def Summand1(y, oscp, osc, lcap):
    np, lp, mp = oscp
    n, l, m = osc
    if m < 0:
        raise ValueError(f"Summand1: m = {m} < 0")
    if mp < 0:
        raise ValueError(f"Summand1: mp = {mp} < 0")
    d = math.factorial(m) * math.factorial(mp)
    x = n - 1 - m;
    if x < 0:
        raise ValueError(f"Summand1: n - 1 - m = {x} < 0")
    d *= math.factorial(x)
    x = np - 1 - mp
    if x < 0:
        raise ValueError(f"Summand1: np - 1 - mp = {x} < 0")
    d *= math.factorial(x)
    return (-1)**(m + mp) / d

def testSummand1():
    y = 0.375
    oscp = (4,2,3)
    osc = (4,2,1)
    lcap = 4
    z = Summand1(y, oscp, osc, lcap)
    zc = 1.0/12.0
    check(abs(z-zc) < 1e-8, f"Summand1({y}, {oscp}, {osc}, {lcap}) = {z}, expecting {zc}")

def Summand2(y, oscp, osc, lcap):
    np, lp, mp = oscp
    n, l, m = osc
    r = scipy.special.loggamma( (l + lp + lcap + 2*m + 2*mp + 3.0)/2.0);
    r -= scipy.special.loggamma(l + m + 1.5);
    r -= scipy.special.loggamma(lp + mp + 1.5);
    return math.exp(r);

def testSummand2():
    z = Summand2(0.375, (4, 2, 3), (4, 2, 1), 4)
    zc = 442 / (7.0 * math.sqrt(math.pi))
    check(abs(z-zc) < 1e-8, f"Summand2(0.375, (4,2,3), (4,2,1), 4) = {z}, expecting {zc}")

def Summand3(y, oscp, osc, lcap):
    np, lp, mp = oscp
    n, l, m = osc
    a = (lcap - lp - l - 2*mp - 2*m) * 0.5
    b = lcap + 1.5
    c = y
    z = scipy.special.hyp1f1(a, b, c)
    return z

def testSummand3():
    z = Summand3(0.375, (4, 2, 3), (4, 2, 1), 4)
    zc = 0.7500960895721925
    check(abs(z-zc) < 1e-8, f"Summand3(0.375, (4,2,3), (4,2,1), 4) = {z}, expecting {zc}")

def BesselFactor3(y, Oscp, Osc, lcap):
    np,lp=Oscp
    n,l=Osc
    tot = 0.0
    for m in range(n):
        for mp in range(np):
            oscp = [np, lp, mp]
            osc = [n, l, m]
            tot += (Summand1(y, oscp, osc, lcap) *
                    Summand2(y, oscp, osc, lcap) *
                    Summand3(y, oscp, osc, lcap))
    return tot

def testBesselFactor3():
    z = BesselFactor3(1.625, (2, 1), (1, 2), 2)
    zc = 0.2056882827554204
    check(abs(z-zc) < 1e-8, f"BesselFactor3(1.625, (2,1), (1,2), 2) = {z}, expecting {zc}")

def BesselElement(y, oscp, osc, lcap):
    np,lp=oscp
    n,l=osc
    return (BesselFactor1(y, oscp, osc, lcap) *
            BesselFactor2(y, oscp, osc, lcap) *
            BesselFactor3(y, oscp, osc, lcap))

def testBesselElement():
    z = BesselElement(1.625, (2, 1), (1, 2), 2)
    zc = 0.05832830085048889
    check(abs(z-zc) < 1e-8, f"BesselElement(1.625, (2,1), (1,2), 2) = {z}, expecting {zc}")
    z = BesselElement(0.125, (3, 0), (2, 2), 3)
    zc = -0.03924307088428147
    check(abs(z-zc) < 1e-8, f"BesselElement(0.125, (3,0), (2,2), 3) = {z}, expecting {zc}")

def BesselFactor1A(y, oscp, osc, lcap):
    np,lp=oscp
    n,l=osc
    r = 2.0**(lcap - 1)
    r /= DFact(JNorm(lcap))
    r *= y**((lcap - 1)/2.0)
    r *= math.exp(-y)
    r *= math.sqrt( math.factorial(np-1) * math.factorial(n-1))
    return r

def testBesselFactor1A():
    z = BesselFactor1A(0.125, (3, 0), (2, 2), 3)
    zc = 0.005943043278035158
    check(abs(z-zc) < 1e-8, f"BesselFactor1A(0.125, (3,0), (2,2), 3) = {z}, expecting {zc}")

def Summand2A(y, oscp, osc, lcap):
    np, lp, mp = oscp
    n, l, m = osc
    r = scipy.special.loggamma((l + lp + lcap + 2*m + 2*mp + 2)*0.5)
    r -= scipy.special.loggamma(l + m + 1.5);
    r -= scipy.special.loggamma(lp + mp + 1.5);
    return math.exp(r)

def testSummand2A():
    q = 0.375
    oscp = (4, 2, 3)
    osc = (4, 2, 1)
    jcap = 4
    z = Summand2A(q, oscp, osc, jcap)
    zc = 131072 / (3465 * math.pi)
    check(abs(z-zc) < 1e-8, f"Summand2A({q}, {oscp}, {osc}, {jcap}) = {z}, expecting {zc}")
    
def Summand3A(y, oscp, osc, lcap):
    np, lp, mp = oscp
    n, l, m = osc
    r = -(l + lp + lcap + 2*m + 2*mp + 2) * 0.5
    a = (lcap - lp - l - 2*mp - 2*m - 1) * 0.5
    b = lcap + 1.5
    c = y
    r *= scipy.special.hyp1f1(a, b, c)
    r += 2*m*scipy.special.hyp1f1(a+1.0,b,c)
    return r

def testSummand3A():
    y = 0.375
    oscp = (6, 2, 3)
    osc = (3, 2, 1)
    lcap = 4
    z = Summand3A(y, oscp, osc, lcap)
    zc = -4.949712716520789
    check(abs(z-zc) < 1e-8, f"Summand3A({y}, {oscp}, {osc}, {lcap}) = {z}, expecting {zc}")

def BesselFactor3A(y, oscp, osc, lcap):
    np,lp=oscp
    n,l=osc
    tot = 0.0
    for m in range(n):
        for mp in range(np):
            oscp = [np, lp, mp]
            osc = [n, l, m]
            tot += (Summand1(y, oscp, osc, lcap) *
                    Summand2A(y, oscp, osc, lcap) *
                    Summand3A(y, oscp, osc, lcap))
    return tot

def testBesselFactor3A():
    y = 0.125
    Oscp = (3, 1)
    Osc = (2, 2)
    lcap = 7
    z = BesselFactor3A(y, Oscp, Osc, lcap)
    zc = 218.41272758911202
    check(abs(z-zc) < 1e-8, f"BesselFactor3A(y, Oscp, Osc, {lcap}) = {z}, expecting {zc}")

def BesselElementMinus(y, oscp, osc, lcap):
    np,lp=oscp
    n,l=osc
    return (BesselFactor1A(y, oscp, osc, lcap) *
            BesselFactor2(y, oscp, osc, lcap) *
            BesselFactor3A(y, oscp, osc, lcap))

def testBesselElementMinus():
    z = BesselElementMinus(1.5, [3, 1], [2, 2], 2)
    zc = 0.011962252530986942
    check(abs(z-zc) < 1e-8, f"BesselElementMinus(1.5, [3,1], [2,2], 2) = {z}, expecting {zc}")
    z = BesselElementMinus(1.75, [4, 1], [3, 0], 1)
    zc = -0.09011903502102633
    check(abs(z-zc) < 1e-8, f"BesselElementMinus(1.75, [4,1], [3,0], 1) = {z}, expecting {zc}")

def Summand4A(y, oscp, osc, lcap):
    np, lp, mp = oscp
    n, l, m = osc
    r = - (l + lp + lcap + 2*m + 2*mp + 2) * 0.5
    a = (lcap - lp - l - 2*mp - 2*m - 1) * 0.5
    b = lcap + 1.5
    c = y
    r *= scipy.special.hyp1f1(a, b, c)
    r += (2*l + 2*m + 1)*scipy.special.hyp1f1(a+1.0,b,c)
    return r

def testSummand4A():
    z = Summand4A(0.375, [6, 2, 3], [3, 2, 1], 4)
    zc = -1.0589870277956415
    check(abs(z-zc) < 1e-8, f"Summand4A(0.375, [6,2,3], [3,2,1], 4) = {z}, expecting {zc}")

def BesselFactor4A(y, oscp, osc, lcap):
    np,lp=oscp
    n,l=osc
    tot = 0.0
    for m in range(n):
        for mp in range(np):
            oscp = [np, lp, mp]
            osc = [n, l, m]
            tot += (Summand1(y, oscp, osc, lcap) *
                    Summand2A(y, oscp, osc, lcap) *
                    Summand4A(y, oscp, osc, lcap))
    return tot

def testBesselFactor4A():
    z = BesselFactor4A(0.975, [3, 3], [3, 2], 7)
    zc = -3.5116176611351797
    check(abs(z-zc) < 1e-8, f"BesselFactor4A(0.975, [3,3], [3,2], 7) = {z}, expecting {zc}")

def BesselElementPlus(y, oscp, osc, lcap):
    np,lp=oscp
    n,l=osc
    return (BesselFactor1A(y, oscp, osc, lcap) *
            BesselFactor2(y, oscp, osc, lcap) *
            BesselFactor4A(y, oscp, osc, lcap))

# REL
# define lower-component Bessel function matrix elements (See Eq. D.16, D26, and D.27 in Evan's thesis)
# y = (q b / 2)^2,  which is positive unless q can be imaginary
def IMat(L, m, y):
    a = math.sqrt(math.pi) / 4.0
    if y < 0:
        print(f"L={L}, y={y}")
        raise ValueError("y is negative")
    a *= math.pow(y, L/2.0) * math.exp(- y)
    b = scipy.special.gamma((L + m + 1.0)/2.0) / scipy.special.gamma(L + 3.0/2.0)
    c = scipy.special.hyp1f1(1 + (L - m)/2.0, L+3.0/2.0, y)
    return a * b * c

# REL
def BesselElementLower1(y, oscp, osc, lcap):
    np,lp=oscp
    n,l=osc
    a = (1.0 / (2.0 * math.sqrt(y)))
    a *= math.sqrt(2.0 * scipy.special.gamma(np) * 2 * scipy.special.gamma(n) / (scipy.special.gamma(np + lp + 0.5) * scipy.special.gamma(n + l + 0.5)))
    tot = 0.0
    for i in range(np):
        for j in range(n):
            tx = scipy.special.binom(np + lp - 0.5, np - i - 1.0) * scipy.special.binom(n + l - 0.5, n - j - 1.0)
            ty = math.pow(-1, i + j) / (math.factorial(i) * math.factorial(j))
            tz = IMat(lcap, 1 + 2*i + 2*j + lp + l, y)
            tot += tx * ty * tz
    return a * tot

# REL
def BesselElementLower2(y, oscp, osc, lcap):
    np,lp=oscp
    n,l=osc
    a = (1.0 / (2.0 * math.sqrt(y)))
    a *= math.sqrt(2.0 * scipy.special.gamma(np) * 2 * scipy.special.gamma(n) / (scipy.special.gamma(np + lp + 0.5) * scipy.special.gamma(n + l + 0.5)))
    tot = 0.0
    for i in range(np):
        for j in range(n):
            tx = scipy.special.binom(np + lp - 0.5, np - i - 1.0) * scipy.special.binom(n + l - 0.5, n - j - 1.0)
            ty = math.pow(-1, i + j) / (math.factorial(i) * math.factorial(j) * (2 * lcap + 1))
            tz = lcap * IMat(lcap-1.0, 2 + 2*i + 2*j + lp + l, y) - (lcap + 1) * IMat(lcap + 1, 2 + 2*i+2*j+lp+l, y)
            tot += tx * ty * tz
    return a * tot

def testBesselElementPlus():
    z = BesselElementPlus(1.5, [3, 1], [2, 2], 2)
    zc = 0.09937871333435762
    check(abs(z-zc) < 1e-8, f"BesselElementPlus(1.5, [3,1], [2,2], 2) = {z}, expecting {zc}")
    z = BesselElementPlus(1.75, [4, 1], [3, 0], 1)
    zc = -0.05094060299913639
    check(abs(z-zc) < 1e-8, f"BesselElementPlus(1.75, [4,1], [3,0], 1) = {z}, expecting {zc}")

def testBessel():
    testBesselFactor1()
    testBesselFactor2()
    testSummand1()
    testSummand2()
    testSummand3()
    testBesselFactor3()
    testBesselElement()
    testBesselFactor1A()
    testSummand2A()
    testSummand3A()
    testBesselFactor3A()
    testBesselElementMinus()
    testSummand4A()
    testBesselFactor4A()
    testBesselElementPlus()

##################################################
# Operator support
##################################################

#
# This is for a single particle state.
# Purpose: convert from NPrincipal (quanta) and j back to L
# j = L \pm 1/2
# but that makes L even or odd based on the choice of \pm
# We know that   NPrincipal - L = 2 * (nodal - 1) must be even
# and that fixes the choice for \pm .
#
def Lnumber(NPrincipal, j):
    a = NPrincipal - (j + 0.5)
    if not IntegerQ(a):
        raise ValueError("Lnumber j+0.5 is not an integer")
    a = round(a)
    return round(j + 0.5 - (a & 1))

# convert from principle (quanta) back to nodal given j
#  NPrincipal = 2 * (nodal - 1) + L
# if we find L, then we can get nodal
#
def Nodal(NPrincipal, j):
    n = (NPrincipal - Lnumber(NPrincipal, j))*0.5 + 1
    return round(n)

#
# Convert from principal quantum, j  oscillator spec to [n, L, j] spec
#
def Osc2n(Osc):
    N,j=Osc
    return (Nodal(N,j), Lnumber(N,j), j)

def ParityState(NPrincipal, j):
    return (-1)**Lnumber(NPrincipal, j)

def ParityNormal(jcap):
    return (-1)**jcap

def ParityConsNormal(NPrincipalp, jp, NPrincipal, j, jcap):
    return ParityState(NPrincipal, j) * ParityState(NPrincipalp, jp) * ParityNormal(jcap)

def testParityConsNormal():
    z = ParityConsNormal(2, 0.5, 2, 1.5, 2)
    zc = 1
    check(z == zc, f"ParityConsNormal(2,1/2,2,3/2,2)={z}, expecting {zc}")

@dispatch(int, float, int, float, int)
def PhysicalConditionsNormal(NPrincipalp, jp, NPrincipal, j, jcap):
    return (ParityConsNormal(NPrincipalp, jp, NPrincipal, j, jcap)==1 and
            (Nodal(NPrincipal, j)>0) and (Nodal(NPrincipalp, jp)>0) and
            (abs(j-jp) <= jcap) and (jcap <= j + jp) )

@dispatch(list, list, int)
def PhysicalConditionsNormal(Oscp, Osc, jcap):
    Np,jp = Oscp
    N,j = Osc
    return PhysicalConditionsNormal(Np, jp, N, j, jcap)

@dispatch(tuple, tuple, int)
def PhysicalConditionsNormal(Oscp, Osc, jcap):
    Np,jp = Oscp
    N,j = Osc
    return PhysicalConditionsNormal(Np, jp, N, j, jcap)

def ParityAbnormal(jcap):
    return (-1)**(jcap + 1)

# look for abnormal parity
def ParityConsAbnormal(NPrincipalp, jp, NPrincipal, j, jcap):
    return ParityState(NPrincipal, j) * ParityState(NPrincipalp, jp) * ParityAbnormal(jcap)

@dispatch(int, float, int, float, int)
def PhysicalConditionsAbnormal(NPrincipalp, jp, NPrincipal, j, jcap):
    return (ParityConsAbnormal(NPrincipalp, jp, NPrincipal, j, jcap)==1 and
            (Nodal(NPrincipal, j)>0) and (Nodal(NPrincipalp, jp)>0) and
            (abs(j-jp) <= jcap) and (jcap <= j + jp) )

@dispatch(list, list, int)
def PhysicalConditionsAbnormal(Oscp, Osc, jcap):
    Np,jp = Oscp
    N,j = Osc
    return PhysicalConditionsAbnormal(Np, jp, N, j, jcap)

@dispatch(tuple, tuple, int)
def PhysicalConditionsAbnormal(Oscp, Osc, jcap):
    Np,jp = Oscp
    N,j = Osc
    return PhysicalConditionsAbnormal(Np, jp, N, j, jcap)

# Operators: Normal parity 
# MJ:
def mjelement(y, oscp, osc, jcap):
    np,lp,jp=oscp
    n,l,j=osc
    r = (-1)**(0.5+j+jcap)
    r *= math.sqrt(JNorm(j)*JNorm(jp)*JNorm(l)*JNorm(lp)*JNorm(jcap)/(4 * math.pi))
    r *= ThreeJ([lp,0], [jcap,0], [l,0])
    r *= SixJ( (lp,jp,0.5), (j,l,jcap))
    r *= BesselElement(y, [np,lp], [n,l], jcap)
    return r

def testmjelement():
    y = 0.625
    oscp = (2, 0, 0.5)
    osc = (1, 2, 1.5)
    jcap = 2
    z = mjelement(y, oscp, osc, jcap)
    zc = -0.1094239693097734
    check(abs(z - zc) < 1e-8, f"mjelement({y}, {oscp}, {osc},{jcap}) = {z}, expecting {zc}")

# rewrite in terms of N,N',j,j'
def MJE(y, Oscp, Osc, jcap):
    return mjelement(y, Osc2n(Oscp), Osc2n(Osc), jcap)

# Explicity force a 0 if physical conditions are not satisfied.
def MJ(y, Oscp, Osc, jcap):
    if PhysicalConditionsNormal(Oscp, Osc, jcap):
        return MJE(y, Oscp, Osc, jcap)
    return 0.0

def testMJ():
    z = MJ(0.625, [2,0.5], [2, 1.5], 2)
    zc = -0.1094239693097734
    check(abs(z - zc) < 1e-8, f"MJ(0.625,[2,0.5], [2,1.5],2) = {z}, expecting {zc}")

def MJLSigma(y, oscp, osc, jcap, lcap):
    np,lp,jp = oscp
    n,l,j = osc
    r = (-1)**lp
    r *= math.sqrt(JNorm(j) * JNorm(jp) * JNorm(l) * JNorm(lp) * JNorm(jcap) * JNorm(lcap) / (4 * math.pi))
    r *= math.sqrt(6.0)
    r *= ThreeJ([lp,0], [lcap, 0], [l,0])
    r *= NineJ([lp,l,lcap], [0.5,0.5,1], [jp, j, jcap])
    r *= BesselElement(y, [np,lp], [n,l], lcap)
    return r

def testMJLSigma():
    z = MJLSigma(2.625, [3,1,1.5], [4, 1, 1.5], 3, 2)
    zc = -0.033672957119242494
    check(abs(z - zc) < 1e-8, f"MJSigma(0.625,[3,1,1.5], [4,1,1.5], 3, 2) = {z}, expecting {zc}")

# Now implement MJLSigma in terms of principal quantum number
def SigmaJE(y, Oscp, Osc, jcap):
    return MJLSigma(y, Osc2n(Oscp), Osc2n(Osc), jcap, jcap)

# Explicitly require that the function is 0 if physical conditions are not met
def SigmaJ(y, Oscp, Osc, jcap):
    if PhysicalConditionsNormal(Oscp, Osc, jcap):
        return SigmaJE(y, Oscp, Osc, jcap)
    return 0.0

def testSigmaJ():
    y = 0.825
    Oscp = (3, 0.5)
    Osc = (5, 1.5)
    jcap = 3
    z = SigmaJ(y, Oscp, Osc, jcap)
    zc = 0
    check(z == zc, f"SigmaJ({y}, {Oscp}, {Osc}, {jcap}) = {z}, expecting {zc}")

    y = 2.625
    Oscp = (3, 1.5)
    Osc = (4, 1.5)
    jcap = 3
    z = SigmaJ(y, Oscp, Osc, jcap)
    zc = 0.010195423929554344
    check(abs(z-zc)<1e-8, f"SigmaJ({y}, {Oscp}, {Osc}, {jcap}) = {z}, expecting {zc}")

# DeltaPJ
def MJLDivQoverall(y, oscp, osc, jcap, lcap):
    np,lp,jp = oscp
    n,l,j=osc
    r = (-1)**(lcap + j + 0.5)
    r *= QNorm(lp) * QNorm(jp) * QNorm(j) * QNorm(jcap) * QNorm(lcap)
    r *= SixJ( (lp,jp,0.5), (j, l, jcap) ) / math.sqrt(4 * math.pi)
    return r

def MJLDivQsummand1(y, oscp, osc, jcap, lcap):
    np,lp,jp = oscp
    n,l,j=osc
    r = -math.sqrt(l + 1) * QNorm(l + 1)
    r *= SixJ( (lcap, 1, jcap), (l, lp, l+1) )
    r *= ThreeJ( [lp, 0], [lcap, 0], [l+1, 0] )
    r *= BesselElementMinus(y, [np, lp], [n, l], lcap)
    return r

def MJLDivQsummand2(y, oscp, osc, jcap, lcap):
    np,lp,jp = oscp
    n,l,j=osc
    if l == 0:
        return 0.0
    r = math.sqrt(l) * QNorm(l - 1)
    r *= SixJ((lcap, 1, jcap), (l, lp, l-1))
    r *= ThreeJ( [lp, 0], [lcap, 0], [l-1, 0])
    r *= BesselElementPlus(y, [np, lp], [n, l], lcap)
    return r

def MJLDivQ(y, oscp, osc, jcap, lcap):
    return ( MJLDivQoverall(y, oscp, osc, jcap, lcap) *
             (MJLDivQsummand1(y, oscp, osc, jcap, lcap) +
             MJLDivQsummand2(y, oscp, osc, jcap, lcap) ))

def DeltaPrime(y, oscp, osc, jcap):
    qc = QNorm(jcap)
    r = - math.sqrt(jcap) / qc
    r *= MJLDivQ(y, oscp, osc, jcap, jcap+1)
    s = math.sqrt(jcap + 1) / qc
    s *= MJLDivQ(y, oscp, osc, jcap, jcap-1)
    return r + s

def testDeltaPrime():
    z = DeltaPrime(0.75, [1,1,0.5], [1,3,2.5], 2)
    zc = 0.10135675506386227
    check(abs(z - zc) < 1e-8, f"DeltaPrime(0.75,[1,1,0.5], [1,3,2.5], 2) = {z}, expecting {zc}")

#
# Same in terms of Principal quantum number
#
def DeltaPrimeJE(y, Oscp, Osc, jcap):
    return DeltaPrime(y, Osc2n(Oscp), Osc2n(Osc), jcap)

#
# Now require that physcial conditions are met to return a non-zero value
def DeltaPJ(y, Oscp, Osc, jcap):
    if PhysicalConditionsNormal(Oscp, Osc, jcap):
        return DeltaPrimeJE(y, Oscp, Osc, jcap)
    else:
        return 0.0

def testDeltaPJ():
    z = DeltaPJ(0.125, [3, 2.5], [7, 4.5], 4)
    zc = 0.01977120073009146
    check(abs(z - zc) < 1e-8, f"DeltaPJ(0.125,[1,3,2.5], [2,5,4.5], 4) = {z}, expecting {zc}")

def DeltaPrimePrime(y, Oscp, Osc, jcap):
    qc = QNorm(jcap)
    r = math.sqrt(jcap + 1) / qc
    r *= MJLDivQ(y, Oscp, Osc, jcap, jcap + 1)
    s = math.sqrt(jcap) / qc
    s *= MJLDivQ(y, Oscp, Osc, jcap, jcap - 1)
    return r + s

def DeltaPrimePrimeJE(y, Oscp, Osc, jcap):
    return DeltaPrimePrime(y, Osc2n(Oscp), Osc2n(Osc),  jcap)

def DeltaPPJ(y, Oscp, Osc, jcap):
    if PhysicalConditionsNormal(Oscp, Osc, jcap):
        return DeltaPrimePrimeJE(y, Oscp, Osc, jcap)
    else:
        return 0

def DeltaTPPJ(y, Oscp, Osc, jcap):
    return DeltaPPJ(y, Oscp, Osc, jcap) - MJ(y, Oscp, Osc, jcap)*0.5

def testDeltaTPPJ():
    y = 0.25
    Oscp = ( 3, 5./2 )
    Osc = ( 7, 9./2 )
    jcap = 4
    z = DeltaTPPJ(y, Oscp, Osc, jcap)
    zc = 0.02148595350676322
    check(abs(z - zc) < 1e-8, f"DeltaTPPJ({y},{Oscp}, {Osc}, {jcap}) = {z}, expecting {zc}")

#
# PhiPPJ
#
def PhiPPoverall(y, oscp, osc, jcap):
    np,lp,jp=oscp
    n,l,j=osc
    return (-1)**(lp + 1) * 6 * QNorm(lp) * QNorm(jp) * QNorm(j) / math.sqrt(4.0 * math.pi)

def PhiPPsummand1(y, oscp, osc, jcap):
    np,lp,jp=oscp
    n,l,j=osc
    r = QNorm(l+1) * math.sqrt(l+1)
    r *= ThreeJ( (lp,0), (jcap+1,0), (l+1,0) )
    r *= BesselElementMinus(y, (np,lp), (n,l), jcap+1)
    if r == 0.0:
        return 0.0
    tot = 0.0
    for lcap in range(jcap, jcap+2): # jcap and jcap+1
        t = (-1)**(jcap - lcap + 1) * JNorm(lcap)
        t *= SixJ( (jcap + 1, 1, lcap), (1, jcap, 1) )
        t *= SixJ( (jcap + 1, 1, lcap), (l, lp, l+1) )
        t *= NineJ( (lp, l, lcap), (0.5, 0.5, 1), (jp, j, jcap) )
        tot += t
    return r * tot

def PhiPPsummand2(y, oscp, osc, jcap):
    np,lp,jp=oscp
    n,l,j=osc
    if l == 0:
        return 0.0
    r = QNorm(l - 1) * math.sqrt(l)
    r *= ThreeJ( (lp, 0), (jcap+1, 0), (l-1, 0) )
    r *= BesselElementPlus(y, (np, lp), (n, l), jcap + 1)
    if r == 0.0:
        return 0.0
    tot = 0.0
    for lcap in range(jcap, jcap + 2):  # jcap and jcap+1
        t = (-1)**(jcap-lcap) * JNorm(lcap)
        t *= SixJ( (jcap + 1, 1, lcap), (1, jcap, 1) )
        t *= SixJ( (jcap + 1, 1, lcap), (l, lp, l-1) )
        t *= NineJ( (lp, l, lcap), (0.5, 0.5, 1), (jp, j, jcap) )
        tot += t
    return r * tot

def PhiPPsummand3(y, oscp, osc, jcap):
    np,lp,jp=oscp
    n,l,j=osc
    if jcap == 0:
        return 0.0
    r = QNorm(l + 1) * math.sqrt(l + 1.0)
    r *= ThreeJ((lp,0),(jcap-1,0),(l+1,0))
    r *= BesselElementMinus(y,(np,lp),(n,l),jcap-1)
    if r == 0.0:
        return 0.0
    tot = 0.0
    for lcap in range(jcap-1, jcap+1): # jcap-1 and jcap
        t = (-1)**(jcap-lcap+1)*JNorm(lcap)
        t *= SixJ( (jcap-1, 1, lcap), (1, jcap, 1) )
        t *= SixJ( (jcap-1, 1, lcap), (l, lp, l+1) )
        t *= NineJ( (lp, l, lcap), (0.5, 0.5, 1), (jp, j, jcap) )
        tot += t
    return r * tot

def PhiPPsummand4(y, oscp, osc, jcap):
    np,lp,jp=oscp
    n,l,j=osc
    if jcap == 0:
        return 0.0
    if l == 0:
        return 0.0
    r = QNorm(l-1) * math.sqrt(l)
    r *= ThreeJ( (lp, 0), (jcap-1,0), (l-1,0) )
    r *= BesselElementPlus(y, (np,lp), (n,l), jcap-1)
    if r == 0.0:
        return 0.0
    tot = 0.0
    for lcap in range(jcap-1, jcap+1):  # jcap-1 and jcap
        t = (-1)**(jcap-lcap) * JNorm(lcap)
        t *= SixJ((jcap-1, 1, lcap), (1, jcap, 1) )
        t *= SixJ((jcap-1, 1, lcap), (l, lp, l-1) )
        t *= NineJ( (lp, l, lcap), (0.5, 0.5, 1), (jp, j, jcap) )
        tot += t
    return r * tot

def PhiPPJX(y, oscp, osc, jcap):
    a = QNorm(jcap+1) * math.sqrt(jcap+1)
    if a != 0.0:
        a *= PhiPPsummand1(y, oscp, osc, jcap) + PhiPPsummand2(y, oscp, osc, jcap)
    b = 0.0
    if jcap != 0:
        b = QNorm(jcap-1) * math.sqrt(jcap)
        if b != 0.0:
            b *= (PhiPPsummand3(y, oscp, osc, jcap) + PhiPPsummand4(y, oscp, osc, jcap))
    return PhiPPoverall(y, oscp, osc, jcap) * (a + b)

def PhiPPJE(y, Oscp, Osc, jcap):
    return PhiPPJX(y, Osc2n(Oscp), Osc2n(Osc), jcap)

def PhiPPJ(y, Oscp, Osc, jcap):
    if PhysicalConditionsNormal(Oscp, Osc, jcap):
        return PhiPPJE(y, Oscp, Osc, jcap)
    else:
        return 0.0

def testPhiPPJ():
    y = 0.25
    Oscp = ( 3, 5./2 )
    Osc = ( 7, 9./2 )
    jcap = 4
    z = PhiPPJ(y, Oscp, Osc, jcap)
    zc = -0.00811043512782808
    check(abs(z - zc) < 1e-8, f"PhiPPJ({y},{Oscp}, {Osc}, {jcap}) = {z}, expecting {zc}")
    if False:
        oscp = Osc2n(Oscp)
        osc = Osc2n(Osc)
        z = PhiPPsummand1(y, oscp, osc, jcap)
        print(f"summand1 = {z}")
        z = PhiPPsummand2(y, oscp, osc, jcap)
        print(f"summand2 = {z}")
        z = PhiPPsummand3(y, oscp, osc, jcap)
        print(f"summand3 = {z}")
        z = PhiPPsummand4(y, oscp, osc, jcap)
        print(f"summand4 = {z}")

# PhiPJ and PhiTPJ
def PhiPJX(y, oscp, osc, jcap):
    a = PhiPPsummand1(y, oscp, osc, jcap) + PhiPPsummand2(y, oscp, osc, jcap)
    a *= -QNorm(jcap+1) * math.sqrt(jcap)
    b = 0
    if jcap != 0:
        b = PhiPPsummand3(y, oscp, osc, jcap) + PhiPPsummand4(y, oscp, osc, jcap)
        b *= QNorm(jcap-1) * math.sqrt(jcap + 1)
    return PhiPPoverall(y, oscp, osc, jcap) * (a + b)

def testPhiPJX():
    y = 0.25
    Oscp = ( 1, 3, 5./2 )
    Osc = ( 2, 5, 9./2 )
    jcap = 4
    z = PhiPJX(y, Oscp, Osc, jcap)
    zc = -0.011606533690899203
    check(abs(z - zc) < 1e-8, f"PhiPJX({y},{Oscp}, {Osc}, {jcap}) = {z}, expecting {zc}")
    

# write in terms of Principle quantum number
def PhiPJE(y, Oscp, Osc, jcap):
    return PhiPJX(y, Osc2n(Oscp), Osc2n(Osc), jcap)

# now impose physical conditions
def PhiPJ(y, Oscp, Osc, jcap):
    if PhysicalConditionsNormal(Oscp, Osc, jcap):
        return PhiPJE(y, Oscp, Osc, jcap)
    else:
        return 0.0

# Now evaluate PhiPT to have simple turn-around properties
def PhiTPJ(y, Oscp, Osc, jcap):
    return PhiPJ(y, Oscp, Osc, jcap) + SigmaJ(y, Oscp, Osc, jcap) * 0.5

def testPhiTPJ():
    y = 0.25
    Oscp = ( 3, 5./2 )
    Osc = ( 7, 9./2 )
    jcap = 4
    z = PhiTPJ(y, Oscp, Osc, jcap)
    zc = -0.010405432375836143
    check(abs(z - zc) < 1e-8, f"PhiTPJ({y},{Oscp}, {Osc}, {jcap}) = {z}, expecting {zc}")

#
# Operators:  Abnormal Parity
#

#
# DeltaJ
#
def Deltaop(y, oscp, osc, jcap):
    return MJLDivQ(y, oscp, osc, jcap, jcap)

# write in terms of Principle quantum number
def DeltaJE(y, Oscp, Osc, jcap):
    return Deltaop(y, Osc2n(Oscp), Osc2n(Osc), jcap)

# impose physical conditions
def DeltaJ(y, Oscp, Osc, jcap):
    if PhysicalConditionsAbnormal(Oscp, Osc, jcap):
        return DeltaJE(y, Oscp, Osc, jcap)
    else:
        return 0.0

def testDeltaJ():
    y = .375
    Oscp = [2, 5./2]
    Osc = [3, 7./2]
    jcap = 4
    z = DeltaJ(y, Oscp, Osc, jcap)
    zc = -0.016876012112710687
    check(abs(z - zc) < 1e-8, f"DeltaJ({y},{Oscp}, {Osc}, {jcap}) = {z}, expecting {zc}")

#
# SigmaPJ
# 
def SigmaPrime(y, oscp, osc, jcap):
    a = -math.sqrt(jcap) / QNorm(jcap)
    a *= MJLSigma(y, oscp, osc, jcap, jcap + 1)
    b = math.sqrt(jcap + 1)/ QNorm(jcap)
    b *= MJLSigma(y, oscp, osc, jcap, jcap - 1)
    return a + b

# write in terms of principal quantum number
def SigmaPrimeJE(y, Oscp, Osc, jcap):
    return SigmaPrime(y, Osc2n(Oscp), Osc2n(Osc), jcap)

# impose physical conditions
def SigmaPJ(y, Oscp, Osc, jcap):
    if PhysicalConditionsAbnormal(Oscp, Osc, jcap):
        return SigmaPrimeJE(y, Oscp, Osc, jcap)
    else:
        return 0.0

#
# SigmaPPJ
#
def SigmaSecond(y, oscp, osc, jcap):
    a = math.sqrt(jcap + 1) * MJLSigma(y, oscp, osc, jcap, jcap + 1)
    b = math.sqrt(jcap) * MJLSigma(y, oscp, osc, jcap, jcap - 1)
    return (a + b) / QNorm(jcap)

# write in terms of principal quantum number
def SigmaSecondJE(y, Oscp, Osc, jcap):
    return SigmaSecond(y, Osc2n(Oscp), Osc2n(Osc), jcap)

# impose physical conditions
def SigmaPPJ(y, Oscp, Osc, jcap):
    if PhysicalConditionsAbnormal(Oscp, Osc, jcap):
        return SigmaSecondJE(y, Oscp, Osc, jcap)
    else:
        return 0

def testSigmaPPJ():
    y = .375
    Oscp = [2, 5./2]
    Osc = [3, 7./2]
    jcap = 4
    z = SigmaPPJ(y, Oscp, Osc, jcap)
    zc = 0.0703260145999258
    check(abs(z - zc) < 1e-8, f"SigmaPPJ({y},{Oscp}, {Osc}, {jcap}) = {z}, expecting {zc}")

#
# OmegaJ and OmegaTJ
#
def MJLDivSigoverall(y, oscp, osc, jcap):
    np,lp,jp=oscp
    n,l,j=osc
    r = (-1)**lp
    r *= QNorm(lp) * QNorm(jp) * QNorm(j) * QNorm(2*j-l) * QNorm(jcap)
    r *= SixJ( (lp, jp, 0.5), (j, 2*j - l, jcap) )
    r *= ThreeJ( (lp, 0), (jcap, 0), (2*j - l, 0)) / math.sqrt(4 * math.pi)
    return r

def MJLDivSigsummand1(y, oscp, osc, jcap):
    np,lp,jp=oscp
    n,l,j=osc
    if j == l + 0.5:  # relying on exact rep of int + 0.5 in floating point
        return -BesselElementMinus(y, (np, lp), (n, l), jcap)
    else:
        return 0

def MJLDivSigsummand2(y, oscp, osc, jcap):
    np,lp,jp=oscp
    n,l,j=osc
    if j == l - 0.5:  # relying on exact rep of int + 0.5 in floating point
        return BesselElementPlus(y, (np, lp), (n, l), jcap)
    else:
        return 0

def MJLDivSig(y, oscp, osc, jcap):
    r = MJLDivSigoverall(y, oscp, osc, jcap)
    r *= (MJLDivSigsummand1(y, oscp, osc, jcap) + MJLDivSigsummand2(y, oscp, osc, jcap))
    return r

#
# Write OmegaJ in terms of principle quantum number
#
def OmegaJE(y, Oscp, Osc, jcap):
    return MJLDivSig(y, Osc2n(Oscp), Osc2n(Osc), jcap)

# Impose physical constraints
def OmegaJ(y, Oscp, Osc, jcap):
    if PhysicalConditionsAbnormal(Oscp, Osc, jcap):
        return OmegaJE(y, Oscp, Osc, jcap)
    else:
        return 0

def OmegaTJ(y, Oscp, Osc, jcap):
    return OmegaJ(y, Oscp, Osc, jcap) + SigmaPPJ(y, Oscp, Osc, jcap) * 0.5

def testOmegaTJ():
    y=0.125
    Oscp = (2, 5./2)
    Osc = (2, 1/2)
    jcap = 3
    z = OmegaTJ(y, Oscp, Osc, jcap)
    zc = -0.016069510098686166
    check(abs(z - zc) < 1e-8, f"OmegaTJ({y},{Oscp}, {Osc}, {jcap}) = {z}, expecting {zc}")


#
# PhiJ and PhiTJ
#
def Phioverall(y, oscp, osc, jcap):
    np,lp,jp=oscp
    n,l,j=osc
    r = (-1)**(lp + 1) * 6
    r *= QNorm(lp) * QNorm(jp) * QNorm(j) * JNorm(jcap) / math.sqrt(4 * math.pi)
    return r

def Phisummand1(y, oscp, osc, jcap):
    np,lp,jp=oscp
    n,l,j=osc
    r = QNorm(l + 1) * math.sqrt(l + 1.0)
    r *= ThreeJ( (lp, 0), (jcap, 0), (l + 1, 0) )
    r *= BesselElementMinus(y, (np, lp), (n, l), jcap)
    tot = 0.0
    for lcap in range( max(0, jcap-1), jcap + 2): # including up to jcap - 1
        t = (-1)**(jcap - lcap + 1) * JNorm(lcap)
        t *= SixJ( (jcap, 1, lcap), (1, jcap, 1) )
        t *= SixJ( (jcap, 1, lcap), (l, lp, l + 1) )
        t *= NineJ( (lp, l, lcap), (0.5, 0.5, 1), (jp, j, jcap) )
        tot += t
    return r * tot

def Phisummand2(y, oscp, osc, jcap):
    np,lp,jp=oscp
    n,l,j=osc
    if l == 0:
        return 0
    r = QNorm(l - 1) * math.sqrt(l)
    r *= ThreeJ( (lp, 0), (jcap, 0), (l - 1, 0))
    r *= BesselElementPlus(y, (np, lp), (n, l), jcap)
    tot = 0.0
    for lcap in range(max(0, jcap-1), jcap+2): # including up to jcap - 1
        t = (-1)**(jcap - lcap) * JNorm(lcap)
        t *= SixJ( (jcap, 1, lcap), (1, jcap, 1) )
        t *= SixJ( (jcap, 1, lcap), (l, lp, l - 1) )
        t *= NineJ( (lp, l, lcap), (0.5, 0.5, 1), (jp, j, jcap) )
        tot += t
    return r * tot

def PhiJX(y, oscp, osc, jcap):
    r = Phioverall(y, oscp, osc, jcap)
    r *= (Phisummand1(y, oscp, osc, jcap) + Phisummand2(y, oscp, osc, jcap))
    return r

# Now write in terms of principle quantum number Np, N
def PhiJE(y, Oscp, Osc, jcap):
    return PhiJX(y, Osc2n(Oscp), Osc2n(Osc), jcap)

# impose physical constraints
def PhiJ(y, Oscp, Osc, jcap):
    if PhysicalConditionsAbnormal(Oscp, Osc, jcap):
        return PhiJE(y, Oscp, Osc, jcap)
    else:
        return 0

# Evaluate PhiT to have simple turn-around properties
def PhiTJ(y, Oscp, Osc, jcap):
    return PhiJ(y, Oscp, Osc, jcap) - SigmaPJ(y, Oscp, Osc, jcap) * 0.5

# REL
# Lower Component Operators: normal parity
# M_J^{(1)} : (Eq. D20 in thesis)
# Note complex value return
#
def m1jelement(y, oscp, osc, jcap):
    np,lp,jp=oscp
    n,l,j=osc
    a = math.pow(-1.0, 0.5 + j + jcap) * math.sqrt(jcap * (jcap + 1.0))
    b = math.sqrt( JNorm(j) * JNorm(jp) * JNorm(l) * JNorm(lp) * JNorm(jcap) / (4*math.pi))
    c = ThreeJ( (lp, 0), (jcap, 0), (l, 0) ) * SixJ( (lp, jp, 0.5), (j, l, jcap) )
    d = BesselElementLower1(y, (np,lp), (n, l), jcap)
    return a * b * c * d

# REL
#  njp: (ncapp, jp),    nj: (ncap, j)
def M1JE(y, njp, nj, jcap):
    ncapp,jp = njp
    ncap,j = nj
    return m1jelement(y, (Nodal(ncapp, jp), Lnumber(ncapp, jp), jp), (Nodal(ncap, j), Lnumber(ncap, j), j), jcap)

# REL
def M1J(y, njp, nj, jcap):
    ncapp,jp = njp
    ncap,j = nj
    if PhysicalConditionsNormal(ncapp, jp, ncap, j, jcap):
        return M1JE(y, njp, nj, jcap)
    else:
        return 0.0

# REL
# M_J^{(2)}: (Eq. D21 in thesis)
def m2jelement(y, oscp, osc, jcap):
    np,lp,jp=oscp
    n,l,j=osc
    a = math.pow(-1.0, 0.5 + j + jcap)
    b = math.sqrt( JNorm(j) * JNorm(jp) * JNorm(l) * JNorm(lp) * JNorm(jcap) / (4*math.pi))
    c = ThreeJ( (lp, 0), (jcap, 0), (l, 0)) * SixJ((lp, jp, 0.5), (j, l, jcap))
    d = BesselElementLower1(y, (np,lp), (n, l), jcap)
    return a * b * c * d

# REL
def M2JE(y, njp, nj, jcap):
    ncapp,jp = njp
    ncap,j = nj
    return m2jelement(y, (Nodal(ncapp, jp), Lnumber(ncapp, jp), jp), (Nodal(ncap, j), Lnumber(ncap, j), j), jcap)

# REL
def M2J(y, njp, nj, jcap):
    ncapp,jp = njp
    ncap,j = nj
    if PhysicalConditionsNormal(ncapp, jp, ncap, j, jcap):
        return M2JE(y, njp, nj, jcap)
    else:
        return 0.0

# Lower Component Operators: abnormal parity
# \Sigma_J^{\prime (0)} : (Eq. D.23 in thesis)
# REL
def MJLSigma0(y, Oscp, Osc, jcap, lcap):
    np,lp,jp=Oscp
    n,l,j=Osc
    a = math.pow(-1, lp)
    b = math.sqrt(JNorm(j)*JNorm(jp)*JNorm(l)*JNorm(lp)*JNorm(jcap)*JNorm(lcap)/ (4 * math.pi))
    c = ThreeJ((lp, 0), (lcap, 0), (l, 0)) * NineJSymbol( (lp, l, lcap), (0.5, 0.5, 1), (jp, j, jcap) )
    d = BesselElement(y, (np, lp), (n, l), jcap)
    return a * b * c * d

# REL
def SigmaPrime0(y, Oscp, Osc, jcap):
    a = math.sqrt(jcap*1.0)/QNorm(jcap)
    a *= MJLSigma0(y, Oscp, Osc, jcap, jcap + 1)
    b = math.sqrt(jcap + 1.0)/QNorm(jcap)
    b *= MJLSigma0(y, Oscp, Osc, jcap, jcap - 1)
    return a + b

# REL
def SigmaPrime0JE(y, njp, nj, jcap):
    ncapp,jp = njp
    ncap,j = nj
    Oscp = (Nodal(ncapp, jp), Lnumber(ncapp, jp), jp)
    Osc = (Nodal(ncap, j), Lnumber(ncap, j), j)
    return SigmaPrime0(y, Oscp, Osc, jcap)

# REL
def SigmaP0J(y, njp, nj, jcap):
    ncapp,jp = njp
    ncap,j = nj
    if PhysicalConditionsAbnormal(ncapp, jp, ncap, j, jcap):
        return SigmaPrime0JE(y, njp, nj, jcap)
    else:
        return 0.0

# \Sigma_J^{\prime\prime (0)}:
# REL
def SigmaPrimePrime0(y, Oscp, Osc, jcap):
    np,lp,jp=Oscp
    n,l,j=Osc
    a = - math.sqrt(jcap + 1.0) / QNorm(jcap)
    a *= MJLSigma0(y, Oscp, Osc, jcap, jcap + 1)
    b = math.sqrt(jcap*1.0) / QNorm(jcap)
    b *= MJLSigma0(y, Oscp, Osc, jcap, jcap - 1)
    return a + b

#REL
def SigmaPrimePrime0JE(y, njp, nj, jcap):
    ncapp,jp = njp
    ncap,j = nj
    Oscp = (Nodal(ncapp, jp), Lnumber(ncapp, jp), jp)
    Osc = (Nodal(ncap, j), Lnumber(ncap, j), j)
    return SigmaPrimePrime0(y, Oscp, Osc, jcap)

#REL
def SigmaPP0J(y, njp, nj, jcap):
    ncapp,jp = njp
    ncap,j = nj
    if PhysicalConditionsAbnormal(ncapp, jp, ncap, j, jcap):
        return SigmaPrimePrime0JE(y, njp, nj, jcap)
    else:
        return 0.0
    
# \Sigma_J^{\prime(2)}: (Eq. D.25 in thesis)
# REL
def MJLSigma2(y, Oscp, Osc, jcap, lcap):
    np,lp,jp=Oscp
    n,l,j=Osc
    a = math.pow(-1.0, lp)
    a *= math.sqrt(JNorm(j) * JNorm(jp) * JNorm(l) * JNorm(lp) * JNorm(jcap) * JNorm(lcap) / (4*math.pi))
    a *= math.sqrt(6.0)
    b = ThreeJ( (lp, 0), (lcap, 0), (l, 0) ) * NineJSymbol( (lp, l, lcap), (0.5, 0.5, 1), (jp, j, jcap) )
    c = BesselElementLower2(y, (np, lp), (n, l), lcap)
    return a * b * c

# REL
def SigmaPrime2(y, Oscp, Osc, jcap):
    np,lp,jp=Oscp
    n,l,j=Osc
    a = -math.sqrt(jcap*1.0)/QNorm(jcap)
    a *= MJLSigma2(y, Oscp, Osc, jcap, jcap + 1)
    b = math.sqrt(jcap+1.0)/QNorm(jcap)
    b *= MJLSigma2(y, Oscp, Osc, jcap, jcap - 1)
    return a + b

# REL
def SigmaPrime2JE(y, njp, nj, jcap):
    ncapp,jp = njp
    ncap,j = nj
    Oscp = (Nodal(ncapp, jp), Lnumber(ncapp, jp), jp)
    Osc = (Nodal(ncap, j), Lnumber(ncap, j), j)
    return SigmaPrime2(y, Oscp, Osc, jcap)

# REL
def SigmaP2J(y, njp, nj, jcap):
    ncapp,jp = njp
    ncap,j = nj
    if PhysicalConditionsAbnormal(ncapp, jp, ncap, j, jcap):
        return SigmaPrime2JE(y, njp, nj, jcap)
    else:
        return 0.0

# \Sigma_J^{\prime\prime(2)}: (Eq. ?? in thesis)
# REL
def SigmaPrimePrime2(y, Oscp, Osc, jcap):
    np,lp,jp=Oscp
    n,l,j=Osc
    a = math.sqrt(jcap+1.0)/QNorm(jcap)
    a *= MJLSigma2(y, Oscp, Osc, jcap, jcap + 1)
    b = math.sqrt(jcap*1.0)/QNorm(jcap)
    b *= MJLSigma2(y, Oscp, Osc, jcap, jcap - 1)
    return a + b

# REL
def SigmaPrimePrime2JE(y, njp, nj, jcap):
    ncapp,jp = njp
    ncap,j = nj
    Oscp = (Nodal(ncapp, jp), Lnumber(ncapp, jp), jp)
    Osc = (Nodal(ncap, j), Lnumber(ncap, j), j)
    return SigmaPrimePrime2(y, Oscp, Osc, jcap)

# REL
def SigmaPP2J(y, njp, nj, jcap):
    ncapp,jp = njp
    ncap,j = nj
    if PhysicalConditionsAbnormal(ncapp, jp, ncap, j, jcap):
        return SigmaPrimePrime2JE(y, njp, nj, jcap)
    else:
        return 0.0

def testPhiTJ():
    y=0.125
    Oscp = (2, 5./2)
    Osc = (2, 1/2)
    jcap = 3
    z = PhiTJ(y, Oscp, Osc, jcap)
    zc = 0.018555471962443736
    check(abs(z - zc) < 1e-8, f"PhiTJ({y},{Oscp}, {Osc}, {jcap}) = {z}, expecting {zc}")
    

def testOperatorSupport():
    print("")
    print("### Testing Operator Support")
    testParityConsNormal()
    testmjelement()
    testMJ()
    testMJLSigma()
    testSigmaJ()
    testDeltaPrime()
    testDeltaPJ()
    testDeltaTPPJ()
    testPhiPPJ()
    testPhiPJX()
    testPhiTPJ()
    testDeltaJ()
    testSigmaPPJ()
    testOmegaTJ()
    testPhiTJ()
    print("")

########################################
# Testing Code
########################################

def testThreeJ():
    print("")
    print("#### Testing ThreeJ")
    z = ThreeJ([0.5, 0.5], [0.5, -0.5], [1.0, 0.0])
    zc = 1.0/math.sqrt(6.0)
    check(abs(z - zc) < 1e-8, f"ThreeJ([1/2,1/2], [1/2,-1/2], [1,0]) = {z}, expecting {zc}")

    z = ThreeJ([1,1], [2,1], [3, -2])
    zc = - math.sqrt(2.0/21.0)
    check(abs(z - zc) < 1e-8, f"ThreeJ([1,1], [2,1], [3,-2]) = {z}, expecting {zc}")
    


def testSixJConditionTriad():
    x = SixJConditionTriad(0.5, 1.5, 1.0)
    check(x, f"6jtriad(1/2, 3/2, 1) = {x}")

    x = SixJConditionTriad(0.5, 1.5, 1.5)
    check(not x, f"6jtriad(1/2, 3/2, 3/2) = {x}")

    x = SixJConditionTriad(0.5, 1.5, 3.0)
    check(not x, f"6jtriad(1/2, 3/2, 3) = {x}")

def testSixJ():
    print("")
    print("#### Testing SixJ")
    testSixJConditionTriad()
    print()
    a = (1, 2, 3)
    b = (2, 1, 2)
    z = SixJ(a, b)
    zc = 1.0 / (5.0 * math.sqrt(21.0))
    check(abs(z - zc) < 1e-8, f"SixJ({a}, {b}) = {z}, expecting {zc}")

def testNineJ():
    print("")
    print("#### Testing NineJ")
    # all arguments are doubled.
    # z = py3nj.wigner9j( 2,  4,  2, 4,  6,  2, 2,  4,  2)
    z = NineJ([1, 2, 1], [2, 3, 1], [1, 2, 1])
    zc = 1.0 / (5.0 * math.sqrt(30.0))
    x = abs(z - zc)
    check(x < 1e-8, "NineJ([1,2,1],[2,3,1],[1,2,1]) = " + f"{z}, expecting {zc}")

    # check that triangle rule violating values
    # yield 0
    # z = py3nj.wigner9j( 2,  4,  2, 4,  6,  2, 2,  4,  20)
    z = NineJ([1,2,1], [2, 3, 1], [1, 2, 10])
    zc = 0.0
    check(z == 0.0, f"NineJ with triangle violation should be 0.0, got {z}")

def loadElements():
    global elements 
    with open("elements.yaml", "r") as f:
        elements = yaml.safe_load(f)
    for name, eht in elements.items():
        Z = eht['Z']
        isotopeht = eht['isotopes']
        idx = 0
        for A, pht in isotopeht.items(): 
            pht['N'] = A - Z
            pht['mass'] = auMeV * A + pht['delta']
            pht['idx'] = idx
            idx = idx + 1

# For debugging print out the elements table
def printElements():
    print("12C mass = ", 12 * auMeV)
    for key, eht in elements.items():
        print(key, eht['symbol'], f" Z={eht['Z']}, Ebind={eht['Ebind']}, Zeff={eht['Zeff']}, Qeff={eht['Qeff']}, fgAvg={eht['fgAvg']}")
        isotopeht = eht['isotopes']
        for A, pht in isotopeht.items(): 
            jzstr=halfstr(pht['JS']) + ','
            print(f"  {A:<2}: N={pht['N']:<2}, Jz={jzstr:<5} mass={pht['mass']:.3f} MeV, ints={pht['interactions']} ")

#
# Figure out location of elastic dir.  We used to 
# look for a subdir of "StructureResults"
# Now we use a direct specification of the Elastic directory
# that contains the one body density matrices for GS to GS
#
def setElasticDir(args):
    global elasticDir
    edir = pathlib.Path(args.edir).expanduser()
    if not edir.is_dir():
        raise ValueError(f"Supplied Elastic density matrix directory {args.edir} is not a directory")
    elasticDir = edir

# Make sure all the Elastic results files exist for all isotopes
def verifyElastic(args):
    global errcnt, elasticDir
    for key, eht in elements.items():
        symbol = eht['symbol']
        isotopeht = eht['isotopes']
        for A, pht in isotopeht.items(): 
            for x in pht['interactions']:
                f = elasticDir /  (symbol + str(A) + "_" + x + ".txt")
                if f.is_file():
                    print(f"Elastic results {f.as_posix()} exists")
                else:
                    f = elasticDir /  (symbol + str(A) + "_" + x.lower() + ".txt")
                    if f.is_file():
                        print(f"Elastic results {f.as_posix()} exists")
                    else:
                        f = elasticDir /  (symbol + str(A) + "_" + x.upper() + ".txt")
                        if f.is_file():
                            print(f"Elastic results {f.as_posix()} exists")
                        else:
                            print(f"Elastic results {f.as_posix()} is missing - fail")
                            errcnt = errcnt + 1

def run_tests():
    testDFact()
    testThreeJ()
    testSixJ()
    testNineJ()
    testBessel()
    testOperatorSupport()
    if errcnt:
        print(f"There were {errcnt} errors!")
    else:
        print(f"All tests passed.")

# key is lower case for lookup.
# first entry is used for file name construction
# second entry is used for reporting.
Interactions = {
    'ck': ['ck', "Cohen and Kurath"],
    'bw': ['bw', "Brown-Wildenthal"],
    'usda': ['usda', "USDA"],
    'usdb': ['usdb', "USDB"],
    '4hw':  ['4hw', "4hw(16O)"],
    '2hw':  ['2hw', "2hw(18O)"],
    'kbp':  ['kbp', "KBP"],
    'gx1a': ['GX1A', "GX1A"],
    'kb3g': ['KB3G', "KB3G"],
    'gcn2850': ['GCN2850', "GCN2850"],
    'jj44b':   ['jj44b', "jj44b"],
    'jun45':   ['JUN45', "JUN45"],
}

#
# Print and verify elements structure
#
def run_pelements():
    printElements()
    verifyElastic(args)
    if errcnt:
        print(f"There were {errcnt} errors!")
    else:
        print(f"All tests passed.")

def getelement(data):
    if 'Element' in data:
        ename = data['Element']
    else:
        raise ValueError("Missing Element tag in input")
    if not ename in elements:
        # might be short name, scan
        for e, edata in elements.items():
            if edata['symbol'] == ename or (('alias' in edata) and (edata['alias'] == ename)):
                ename = e # convert
                break
        if not ename in elements:
            raise ValueError(f"unknown element name or symbol {ename}")
        data['Element'] = ename
    return ename

def doabundance(data, print_details = False):
    ename = data['Element']
    edata = elements[ename]
    # get number of keys in isotope dict - the number of isotopes
    icnt = len(edata['isotopes'].keys()) # isotope dict
    AbNorm = np.zeros(icnt) # will be fractional data
    x = data['Isotope']
    itab = edata['isotopes']
    if x <= 0:
        # average over all isotopes
        atot = 0.0
        for a, idata in itab.items():
            atot += idata['abundance']
        for a, idata in itab.items():
            idx = idata['idx'] # get index into AbNorm
            AbNorm[idx] = idata['abundance'] / atot
            if print_details:
                print(f"Target is {a}{ename} at normalized abundance {AbNorm[idx]:.6f}")
    else:
        if not x in itab:
            ValueError("Asking for unknown isotope {x} of {ename}")
        idata = itab[x]
        idx = idata['idx'] # get index into AbNorm
        AbNorm[idx] = 1.0
        if print_details:
            print(f"Target is {x}{ename} at abundance 1.0")
    data['AbNorm'] = AbNorm
    # Also calculate Abar the average A, weighted by abundance AbNorm
    Abar = 0.0
    for a, idata in itab.items():
        idx = idata['idx']
        Abar += AbNorm[idx] * a
    data['Abar'] = Abar
    if print_details:
        print(f"  Average A = {Abar:0.6f}")

def getinteraction(data):
    ename = data['Element']
    if not 'Interaction' in data:
        ValueError("Interaction missing from data")
    interaction = data['Interaction'].lower()
    if not interaction in Interactions:
        ValueError(f"interaction {interaction} is an unknown interaction")
    AbNorm = data['AbNorm']
    icnt = len(AbNorm)
    edata = elements[ename]
    itab = edata['isotopes'] # isotope table
    for a, idata in itab.items():
        idx = idata['idx']
        if AbNorm[idx] > 0.0:
            # check that interaction is in the interactions list
            found = False
            for istr in idata['interactions']:
                if istr.lower() == interaction:
                    found = True
            if not found:
                raise ValueError(f"Interaction {interaction} is not valid for {a}{ename}")
    # Note:  case appears to be mixed for files on disk
    data['Interaction'] = Interactions[interaction][0] # Use table to pick case to match files.
    data['InteractionReport'] = Interactions[interaction][1] # shift to lower case

def getoscb(data, print_details = False):
    ename = data['Element']
    Abar = data['Abar']
    enamel = ename.lower()
    if enamel == 'carbon':
        print("  Carbon oscillator length scale b fixed at 1.77 fm")
        b = 1.77
    elif enamel == 'fluorine':
        print("  Fluorine oscillator length scale b fixed at 1.833 fm")
        b = 1.833
    else:
        if ('oscb' in data) and (data['oscb'] > 0.0):
            b = data['oscb']
            if print_details:
                print(f"  Using user supplied oscillator length scale b={b} fm")
        else:
            b = math.sqrt(41.467 / (45.0 * Abar**(-1.0/3.0) - 25.0 * Abar**(-2.0/3.0)))
            if print_details:
                print(f"  {ename} oscillator length scale set to {b:.6f} fm")
                print("     *Override by setting oscb value in input data to value > 0")
    data['oscb'] = b

#
# We computed the mass for each isotope on load of elements.yaml
# Here we compute the average mass according to the fractional
# abundances AbNorm.
#
def getmasses(data, print_details = False):
    ename = data['Element']
    AbNorm = data['AbNorm']
    edata = elements[ename]
    itab = edata['isotopes'] # isotope table
    Mbar = 0.0
    for a, idata in itab.items():
        idx = idata['idx']
        Mbar += AbNorm[idx] * idata['mass']
    data['Mbar'] = Mbar
    if print_details:
        print(f"  Average Nuclear Mass = {Mbar:0.6f}")

#
# Computing random energy, momentum, charge, radii used
# in formulas and saving in dict data
#
def computeThings(data, print_details = False):
    ename = data['Element']
    Mbar = data['Mbar']
    edata = elements[ename]
    Ebind = edata['Ebind']
    data['Ebind'] = Ebind
    qval = math.sqrt((Mbar / (mumass+Mbar))*((mumass - Ebind)**2 - emass**2))
    data['qval'] = qval
    data['qeff'] = edata['Qeff']
    mu = mumass * Mbar / (mumass + Mbar)
    data['mu'] = mu
    data['fgAvg'] = edata['fgAvg'] # f/g muon component average over isotopes
    Z = edata['Z']
    data['Z'] = Z
    Zeff = edata['Zeff']
    data['Zeff'] = Zeff
    data['RZ2'] = (Zeff * alpha * mu)**3 / math.pi
    if print_details:
        print(f"  Z = {Z}, Zeff = {Zeff}, RZ2 = {data['RZ2']}")
        print(f"  Momentum transfer to electron = {qval} MeV   qeff = {data['qeff']} MeV")
        print(f"  muon binding energy = {Ebind}")
        print(f"  Ratio of average muon Dirac components <f>/<g> = {data['fgAvg']}")
    
# Matrix element files (one body) are grouped into sets
# organized by Final state angular momentum JF and isospin TF,  
#              Initial state angular momentum JI and isospin TI
#       and    Operator  J0 and T0   (not sure why mathematica code uses 0 instead of O)     
#
meset = namedtuple('meset', ['JF', 'TF', 'JI', 'TI', 'J0', 'T0', 'rows'])
# Within a group we have rows which identify a bra (final) oscillator state
# by principle quantum number (quanta) and angular momentum 
# followed by the value of the matrix element.
# bra is a tuple  (N, j) and ket is a tuple (N, j)
merow = namedtuple('merow', ['bra', 'ket', 'val'])

def ISOT(TT, MT, T0):
    r = 0.0
    if T0 == 0:
        r += math.sqrt(2.0) / QNorm(TT)
    if T0 == 1 and MT != 0:
        r += MT * math.sqrt(6 / ((2.0 * TT + 1.0) * (TT + 1.0) * TT) )
    return r

# ii is the selected isotope
def FM(data, a, meset, y):
    ename = data['Element']
    edata = elements[ename]  # element data
    itab = edata['isotopes'] # table of isotopes
    idata = itab[a] # isotope data
    # print(f"a={a}, idata = {idata}")
    r = ISOT(idata['TTMS'][0], idata['TTMS'][1], meset.T0)
    tot = 0.0
    for row in meset.rows:
        tot += row.val * MJ(y, row.bra, row.ket, round(meset.J0))
    return r * tot

def testFM(data):
    if data['Element'] != "Sulfur":
        print("testFM needs Sulfer selected")
        return
    xdata = copy.deepcopy(data)
    xdata['isochar'] = [0.5, 0.5] # to match ref data
    a = 32 # an isotope of Sulpher
    z = FM(xdata, a, data['isotopeme'][a][0], 0.5)
    zc = 3.196732746006806
    check(abs(z - zc) < 1e-8, f"FM(isochar=[0.5,0.5], y=0.5) = {z}, expecting {zc}")
    a = 33 # an isotope of Sulpher
    z = FM(xdata, a, data['isotopeme'][a][2], 0.5)
    zc = 0.0
    check(abs(z - zc) < 1e-8, f"FM(isochar=[0.5,0.5], y=0.5) = {z}, expecting {zc}")

    a = 33 # an isotope of Sulpher
    z = FM(xdata, a, data['isotopeme'][a][4], 0.5)
    zc = -0.22721850919643283
    check(abs(z - zc) < 1e-8, f"FM(isochar=[0.5,0.5], y=0.5) = {z}, expecting {zc}")

    print("testFM passed")

#
# This is a check of the matrix elements we read.
# The J0=0, T0=0 result should give A
# The J0=0, T0=1 result apparently gives Z-N
def evalsumrules(data, print_details = False):
    if print_details:
        print("# Sum rules for 1 body density matrices")
        print("#  J0=0, T0=0 result should give a charge of A (nucleon count)")
        print("#  J0=0, T0=1 result should give a charge of (Z - N)")
    ename = data['Element']
    symbol = data['symbol']
    AbNorm = data['AbNorm']
    edata = elements[ename]
    itab = edata['isotopes']
    for a, isme in data['isotopeme'].items():
        idata = itab[a] # get isotope data
        idx = idata['idx']
        if AbNorm[idx] > 0.0:
            if print_details:
                print(f"  Evaluating sum rules for {a}{symbol}, Z={edata['Z']}, N={idata['N']}")
            for meset in data['isotopeme'][a]:
                if meset.J0 == 0:
                    Valu = FM(data, a, meset, 0.00000001)
                    Valu *= math.sqrt(4.0 * math.pi / (2.0 * meset.JI + 1))
                    if print_details:
                        print(f"    J0 = {meset.J0}, T0 = {meset.T0}, Charge = {Valu:0.2f}")

#
# Read interaction files
# We verify that the matrix elements are in the ground state of the isotope
def readint(data, print_details = False):
    ename = data['Element']
    edata = elements[ename]
    intstr = data['Interaction']
    itab = edata['isotopes']
    icnt = len(edata['isotopes'].keys()) # isotope dict size, range of idx
    data['isotopeme'] = dict()
    iso = data['Isotope']
    for a, idata in itab.items(): # a is nucleon count, idata is isotope data
        if iso and iso != a:
            continue
        idx = idata['idx']
        erstr = edata['symbol'] + str(a) + "_" + intstr + ".txt"
        intfile = elasticDir / erstr
        if not intfile.exists():
            intfile = elasticDir / (edata['symbol'] + str(a) + "_" + intstr.lower() + ".txt")
            if not intfile.exists():
                intfile = elasticDir / (edata['symbol'] + str(a) + "_" + intstr.upper() + ".txt")
                if not intfile.exists():
                    raise ValueError(f"Can't locate elastic results {erstr}")
            
        print(intfile)
        curmeset = None
        mesets = []
        with intfile.open() as f:
            lines = f.read().splitlines()
        for line in lines:
            tokens = line.split()
            tlen = len(tokens)
            if tlen == 0:  # empty line
                continue
            elif tlen == 6:  # meset
                jf = float(tokens[0]) * 0.5 # all half integer values
                tf = float(tokens[1]) * 0.5
                ji = float(tokens[2]) * 0.5
                ti = float(tokens[3]) * 0.5
                j0 = float(tokens[4]) * 0.5
                t0 = float(tokens[5]) * 0.5
                curmeset = meset(JF=jf, TF=tf, JI=ji, TI=ti, J0=j0, T0=t0, rows=[])
                # check GS match for isotope
                if jf == idata['JS'] and tf == idata['TTMS'][0] and ji == idata['JS'] and ti == idata['TTMS'][0]:
                    mesets.append(curmeset)
                else:
                    print(f"Ignoring: {curmeset}, doesn't match gs spec for isotope")
            elif tlen == 5: # merow
                nb = int(tokens[0])
                jb = float(tokens[1])*0.5
                nk = int(tokens[2])
                jk = float(tokens[3])*0.5
                me = float(tokens[4])
                r = merow(bra=(nb, jb), ket=(nk, jk), val=me)
                curmeset.rows.append(r)
            elif tlen == 1: # last merow
                num = int(tokens[0])
                if num >= 0:
                    ValueError(f"Expecting end of meset indicator, got {tokens[0]}")
            else:
                raise ValueError(f"Bad Format interaction file {intfile.as_posix()}")
        if print_details:
            for mes in mesets:
                print(f"meset: JF={mes.JF}, TF={mes.TF}, JI={mes.JI}, TI={mes.TI}, J0={mes.J0}, T0={mes.T0}")
                for me in mes.rows:
                    print(f"  merow: bra={me.bra}, ket={me.ket}, val={me.val}")
        data['isotopeme'][a] = mesets

#
# Note:  in the Mathematica version it uses Simplify instead of direct
# calculation.  Probably to preserve symbols like Pi and Sqrt[2] ...
# Here we just do it numerically
#
# In the Mathematica code these functions take
#   ii - the isotope index index
#   jj - the matrix elements under a specifc   JF, TF, JI, TI, J0, T0
#    y - a float  (I have to remind myself of what this is. 
#                  It flows through into the Bessel function factors)
#   
# In the python version data replaces the globals used in the Mathematica code
#  data holds the element, isotope info including matrix element data for the selected
#  interaction.
# a identifies the isotope with the nucleon count
# ms - the meset we are interested in - conceptually the same as jj.
# y - a float needed for lower level functions.  In range [0,2]
#

# A common preable to functions
def getidata(data, a):
    ename = data['Element']
    edata = elements[ename]  # element data
    itab = edata['isotopes'] # table of isotopes
    idata = itab[a] # isotope data
    return idata

def evalOp(data, a, ms, y, Op):
    idata = getidata(data, a)
    r = ISOT(idata['TTMS'][0], idata['TTMS'][1], ms.T0)
    tot = 0.0
    for row in ms.rows:
        tot += row.val * Op(y, row.bra, row.ket, round(ms.J0))
    return r * tot

#
# Spin dependent operators
#
def FPhiPP(data, a, ms, y):
    return evalOp(data, a, ms, y, PhiPPJ)

def FSigmaPP(data, a, ms, y):
    return evalOp(data, a, ms, y, SigmaPPJ)

def FDelta(data, a, ms, y):
    return evalOp(data, a, ms, y, DeltaJ)

def FSigmaP(data, a, ms, y):
    return evalOp(data, a, ms, y, SigmaPJ)

def FPhiTP(data, a, ms, y):
    return evalOp(data, a, ms, y, PhiTPJ)

#
# Lower Component Operators:
# REL
def FM1(data, a, ms, y):
    return evalOp(data, a, ms, y, M1J)

#REL
def FM2(data, a, ms, y):
    return evalOp(data, a, ms, y, M2J)

#REL
def FSigmaP0(data, a, ms, y):
    return evalOp(data, a, ms, y, SigmaP0J)

#REL
def FSigmaPP0(data, a, ms, y):
    return evalOp(data, a, ms, y, SigmaPP0J)

#REL
def FSigmaP2(data, a, ms, y):
    return evalOp(data, a, ms, y, SigmaP2J)

#REL
def FSigmaPP2(data, a, ms, y):
    return evalOp(data, a, ms, y, SigmaPP2J)

def ResOp(data, y, Op):
    AbNorm = data['AbNorm']
    ename = data['Element']
    isod = data['isochar']
    edata = elements[ename]
    intstr = data['Interaction']
    itab = edata['isotopes']
    tot = 0.0
    for a, idata in itab.items(): # going through isotopes
        idx = idata['idx'] # index to AbNorm
        ab = AbNorm[idx] * 4 * math.pi
        mslist = data['isotopeme'][a]  # list of matrix element sets
        # note: double loop because we are squaring and if the
        # first Op piece and second Op piece can reach the same state we have to
        # count them both
        for ms in mslist: # for each JF, TF, JI, TI, J0, T0 set.
            # print(f"ms JF={ms.JF}, TF={ms.TF}, JI={ms.JI}, TI={ms.TI}, J0={ms.J0}, T0={ms.T0}")
            for msj in mslist: # double loop
                # print(f"msj JF={msj.JF}, TF={msj.TF}, JI={msj.JI}, TI={msj.TI}, J0={msj.J0}, T0={msj.T0}")
                T0j = round(msj.T0)
                if ms.J0 == msj.J0:
                    T0 = round(ms.T0)
                    abms = ab / (2 * ms.JI + 1.0);
                    r = abms
                    r *= Op(data, a, ms, y)
                    r *= Op(data, a, msj, y)
                    r *= isod[T0] # note mathematica arrays start with 1
                    r *= isod[T0j]
                    tot += r
    return tot

#
# Wrapper functions for plotting
# response vs y
#
def ResM(data, y):
    return ResOp(data, y, FM)

def ResSigmaPP(data, y):
    return ResOp(data, y, FSigmaPP)

def ResSigmaP(data, y):
    return ResOp(data, y, FSigmaP)

def ResSD(data, y):
    return ResSigmaPP(data, y) + ResSigmaP(data, y)

def ResDelta(data, y):
    return ResOp(data, y, FDelta)

def ResPhiPP(data, y):
    return ResOp(data, y, FPhiPP)

def ResPhiTP(data, y):
    return ResOp(data, y, FPhiTP)


#
# Effective theory input connected to operators
# For compatability of the indices with the Mathematica
# code we will use indices starting at 1
#
def RM(i, j, cs):
    return cs[1, i] * cs[1, j].conjugate() + cs[11, i] * cs[11, j].conjugate()

def RPhiPP(i, j, cs):
    return cs[3, i] * cs[3, j].conjugate() + (cs[12,i] - cs[15,i])*(cs[12,j]-cs[15,j]).conjugate()

def RPhiPPM(i, j, cs):
    x = cs[3,i] * cs[1,j].conjugate() - (cs[12,i]-cs[15,i])*cs[11,j].conjugate()
    return x.real

def RPhiTP(i,j,cs):
    return cs[12,i]*cs[12,j].conjugate() + cs[13,i]*cs[13,j].conjugate()

def RSigmaPP(i,j,cs):
    return cs[10,i]*cs[10,j].conjugate() + (cs[4,i]-cs[6,i])*(cs[4,j]-cs[6,j]).conjugate()

def RSigmaP(i,j,cs):
    return cs[4,i]*cs[4,j].conjugate() + cs[9,i]*cs[9,j].conjugate()

def RDelta(i,j,cs):
    return cs[5,i]*cs[5,j].conjugate() + cs[8,i]*cs[8,j].conjugate()

def RDeltaSigmaP(i,j,cs):
    x = cs[5,i]*cs[4,j].conjugate() + cs[8,i]*cs[9,j].conjugate()
    return x.real

# Lower Component Operators: (Eq. B4 in long paper)
# REL
def RM1(i,j,bs):
    return bs[3,i]*bs[3,j].conjugate() + bs[7,i]*bs[7,j].conjugate()

#REL
def RM2(i,j,bs):
    return bs[2,i]*bs[2,j].conjugate() + bs[7,i]*bs[7,j].conjugate()

#REL
def RSigmaP0(i,j,bs):
    q = bs[12,i] - bs[15,i]
    t = bs[12,j] - bs[15,j]
    return q * t.conjugate() + bs[14,i]*bs[14,j].conjugate()

#REL
def RSigmaP2(i,j,bs):
    q = bs[13,i] - bs[14,i]
    t = bs[13,j] - bs[14,j]
    return q * t.conjugate() + bs[15,i]*bs[15,j].conjugate()

#REL
def RSigmaPP0(i,j,bs):
    return bs[8,i]*bs[8,j].conjugate() + bs[14,i]*bs[14,j].conjugate()

#REL
def RSigmaPP2(i,j,bs):
    q = bs[13,i] - bs[14,i]
    t = bs[13,j] - bs[14,j]
    return q * t.conjugate() + bs[16,i]*bs[16,j].conjugate()

#REL
def RSigmaP0SigmaP2(i,j,bs):
    q = bs[13,i] - bs[14,i]
    t = bs[13,j] - bs[14,j]
    x = bs[14,i]*t.conjugate() * t.conjugate() + q * bs[15,j].conjugate()
    return x.real
    
#REL
def RSigmaPP0SigmaPP2(i,j,bs):
    t = bs[13,j] - bs[14,j]
    x = bs[8,i] * bs[16,j].conjugate() + bs[13,i] * t.conjugate()
    return x.real

# Upper/Lower interference:

#REL
def RMM2(i,j,cs,bs):
    a = cs[1,i] * bs[2,j].conjugate() - cs[11,i]*bs[7,j].conjugate()
    return a.imag

#REL
def RPhiPPM2(i,j,cs,bs):
    a = (cs[12,i] - cs[15,i])*bs[7,j].conjugate() + cs[3,i] * bs[2,j].conjugate()
    return a.imag

#REL
def RSigmaP0SigmaP(i,j,cs,bs):
    a = bs[14,i]*cs[4,j].conjugate() - (bs[12,i]-bs[15,i])*cs[9,j].conjugate()
    return a.imag

# REL
def RSigmaP2SigmaP(i,j,cs,bs):
    a = -(bs[13,i]-bs[14,i])*cs[4,j].conjugate() + bs[15,i] * cs[9,j].conjugate()
    return a.imag

# REL
def RDeltaSigmaP0(i,j,cs,bs):
    t = bs[12,j] - bs[15,j]
    a = cs[5,i]*bs[14,j].conjugate() - cs[8,i] * t.conjugate()
    return a.imag

# REL
def RDeltaSigmaP2(i,j,cs,bs):
    t = bs[13,j] - bs[14,j]
    a = cs[8,i] * bs[15,j].conjugate() - cs[5,i] * t.conjugate()
    return a.imag

# REL
def RPhiTPM1(i,j,cs,bs):
    a = cs[12,i]*bs[7,j].conjugate() + cs[13,i]*bs[3,j].conjugate()
    return a.imag

# REL
def RSigmaPP0SigmaPP(i,j,cs,bs):
    t = cs[4,j] - cs[6,j]
    a = bs[13,i] * t.conjugate() - bs[8,i] * cs[10,j].conjugate()
    return a.imag

# REL
def RSigmaPP2SigmaPP(i,j,cs,bs):
    t = cs[4,j] - cs[6,j]
    a = (bs[13,i] - bs[14,i]) * t.conjugate() - bs[16,i] * cs[10,j].conjugate()
    return a.imag
#
# We now calculate the total response functions, summed over isotopes,
# weighted by the abundance,again with very small kinematic differences
# between isotopes due to masses ignored, to be later approximated by 
# an effecive average nuclear mass
#
# Note factorization in to nuclear F* and leptonic R* functions.
#
def Response(data, y, cs, FbOp, FkOp, ROp):
    AbNorm = data['AbNorm']
    ename = data['Element']
    isod = data['isochar']
    edata = elements[ename]
    intstr = data['Interaction']
    itab = edata['isotopes']
    tot = 0.0
    for a, idata in itab.items(): # going through isotopes
        idx = idata['idx'] # index to AbNorm
        ab = AbNorm[idx]
        if ab == 0.0:
            continue
        ab *= 4 * math.pi
        mslist = data['isotopeme'][a]  # list of matrix element sets
        for ms in mslist:
            T0 = round(ms.T0)
            for msj in mslist:
                if ms.J0 != msj.J0:
                    continue
                T0j = round(msj.T0)
                r = ab / (2 *ms.JI + 1)
                r *= FkOp(data, a, ms, y)
                r *= FbOp(data, a, msj, y)
                r *= ROp(T0, T0j, cs)
                tot += r
    return tot

# Version for cross from top to bottom
def Response2(data, y, cs, bs, FbOp, FkOp, ROp):
    AbNorm = data['AbNorm']
    ename = data['Element']
    isod = data['isochar']
    edata = elements[ename]
    intstr = data['Interaction']
    itab = edata['isotopes']
    tot = 0.0
    for a, idata in itab.items(): # going through isotopes
        idx = idata['idx'] # index to AbNorm
        ab = AbNorm[idx]
        if ab == 0.0:
            continue
        ab *= 4 * math.pi
        mslist = data['isotopeme'][a]  # list of matrix element sets
        for ms in mslist:
            T0 = round(ms.T0)
            for msj in mslist:
                if ms.J0 != msj.J0:
                    continue
                T0j = round(msj.T0)
                r = ab / (2 *ms.JI + 1)
                r *= FkOp(data, a, ms, y)
                r *= FbOp(data, a, msj, y)
                r *= ROp(T0, T0j, cs, bs)
                tot += r
    return tot



# Note:  no use in qm in this case
def WM(data, qm, y, cs):
    return Response(data, y, cs, FM, FM, RM)

def WPhiPP(data, qm, y, cs):
    return qm**2 * Response(data, y, cs, FPhiPP, FPhiPP, RPhiPP)

def WPhiPPM(data, qm, y, cs):
    return -2.0 * qm * Response(data, y, cs, FM, FPhiPP, RPhiPPM)

def WPhiTP(data, qm, y, cs):
    return qm**2 * Response(data, y, cs, FPhiTP, FPhiTP, RPhiTP)

def WSigmaPP(data, qm, y, cs):
    return Response(data, y, cs, FSigmaPP, FSigmaPP, RSigmaPP)

def WDelta(data, qm, y, cs):
    return qm**2 * Response(data, y, cs, FDelta, FDelta, RDelta)

def WSigmaP(data, qm, y, cs):
    return Response(data, y, cs, FSigmaP, FSigmaP, RSigmaP)

def WDeltaSigmaP(data, qm, y, cs):
    return -2.0 * qm * Response(data, y, cs, FDelta, FSigmaP, RDeltaSigmaP)

# Lower Component Operators: (Eq. B3, B5 in the long paper)
# REL
def WM1(data, qm, y, bs):
    return Response(data, y, bs, FM1, FM1, RM1)

# REL
def WM2(data, qm, y, bs):
    return Response(data, y, bs, FM2, FM2, RM2)

# REL
def WSigmaP0(data, qm, y, bs):
    return Response(data, y, bs, FSigmaP0, FSigmaP0, RSigmaP0)

# REL
def WSigmaP2(data, qm, y, bs):
    return Response(data, y, bs, FSigmaP2, FSigmaP2, RSigmaP2)

# REL
def WSigmaPP0(data, qm, y, bs):
    return Response(data, y, bs, FSigmaPP0, FSigmaPP0, RSigmaPP0)

# REL
def WSigmaPP2(data, qm, y, bs):
    return Response(data,y, bs, FSigmaPP2, FSigmaPP2, RSigmaPP2)

# Now for the top to bottom component cases
# REL
def WMM2(data, qm, y, cs, bs):
    return 2 * Response2(data, y, cs, bs, FM, FM2, RMM2)
 
# REL
def WPhiPPM2(data, qm, y, cs, bs):
    return -2 * qm * Response2(data, y, cs, bs, FPhiPP, FM2, RPhiPPM2)

# REL
def WSigmaP0SigmaP(data, qm, y, cs, bs):
    return -2 * Response2(data, y, cs, bs, FSigmaP0, FSigmaP, RSigmaP0SigmaP)

# REL
def WSigmaP2SigmaP(data, qm, y, cs, bs):
    return -2 * Response2(data, y, cs, bs, FSigmaP2, FSigmaP2, RSigmaP2SigmaP)

# REL
def WSigmaP0SigmaP2(data, qm, y, bs):
    return -2 * Response(data, y, bs, FSigmaP0, FSigmaP2, RSigmaP0SigmaP2)

# REL
def WDeltaSigmaP0(data, qm, y, cs, bs):
    return -2 * qm * Response2(data, y, cs, bs, FDelta, FSigmaP0, RDeltaSigmaP0)

# REL
def WDeltaSigmaP2(data, qm, y, cs, bs):
    return -2 * qm * Response2(data, y, cs, bs, FDelta, FSigmaP2, RDeltaSigmaP2)

# REL
def WPhiTPM1(data, qm, y, cs, bs):
    return -2 * qm * Response2(data, y, cs, bs, FPhiTP, FM1, RPhiTPM1)

# REL
def WSigmaPP0SigmaPP(data, qm, y, cs, bs):
    return 2 * Response2(data, y, cs, bs, FSigmaPP0, FSigmaPP, RSigmaPP0SigmaPP)

# REL
def WSigmaPP2SigmaPP(data, qm, y, cs, bs):
    return 2 * Response2(data, y, cs, bs, FSigmaPP2, FSigmaPP, RSigmaPP2SigmaPP)

# REL
def WSigmaPP0SigmaPP2(data, qm, y, bs):
    return 2 * Response(data, y, bs, FSigmaPP0, FSigmaPP2, RSigmaPP0SigmaPP2)


# REL
# updated for lower component contributions (if bs is set)
def Heff(data, qm, y, cs, bs = None):
    r =  WM(data, qm, y, cs)
    r += WPhiPP(data, qm, y, cs)
    r += WPhiPPM(data, qm, y, cs)
    r += WPhiTP(data, qm, y, cs)
    r += WSigmaPP(data, qm, y, cs)
    r += WDelta(data, qm, y, cs)
    r += WSigmaP(data, qm, y, cs)
    r += WDeltaSigmaP(data, qm, y, cs)
    if not (bs is None):
        fgAvg = data['fgAvg']
        # lower-lower
        t  = WM1(data, qm, y, bs)
        t += WM2(data, qm, y, bs)
        t += WSigmaP0(data, qm, y, bs)
        t += WSigmaP2(data, qm, y, bs)
        t += WSigmaPP0(data, qm, y, bs)
        t += WSigmaPP2(data, qm, y, bs)
        t += WSigmaP0SigmaP2(data, qm, y, bs)
        t += WSigmaPP0SigmaPP2(data, qm, y, bs)
        t *= fgAvg*fgAvg
        # cross
        s  = WMM2(data, qm, y, cs, bs)
        s += WPhiPPM2(data, qm, y, cs, bs)
        s += WPhiTPM1(data, qm, y, cs, bs)
        s += WSigmaP0SigmaP(data, qm, y, cs, bs)
        s += WSigmaP2SigmaP(data, qm, y, cs, bs)
        s += WDeltaSigmaP0(data, qm, y, cs, bs)
        s += WDeltaSigmaP2(data, qm, y, cs, bs)
        s += WSigmaPP0SigmaPP(data, qm, y, cs, bs)
        s += WSigmaPP2SigmaPP(data, qm, y, cs, bs)
        s *= fgAvg
        r += t + s
    return r
    
def testWM(data):
    if data['Element'] != "Copper":
        print("testWM needs Copper selected")
        return
    print("*** TestWM running on Cu isotope")
    # needs 63 Cu data loaded
    qm = 0.12087646432374868
    y = 0.3525021049904227
    cs = np.array([[0. , 0. ],
       [0.2, 0.2],
       [0. , 1. ],
       [0. , 0. ],
       [0. , 0. ],
       [0.3, 0. ],
       [0. , 0. ],
       [0. , 0. ],
       [0. , 0. ],
       [0. , 0. ],
       [0. , 0. ],
       [0. , 0. ],
       [0. , 0. ],
       [0. , 0. ],
       [0. , 0. ],
       [0. , 0. ],
       [0. , 1. ]]);
    z = WM(data, qm, y, cs)
    zc = 22.934564273016356
    check(abs(z - zc) < 1e-8, f"WM(data, {qm}, {y}, {cs}) = {z}, expecting {zc}")

# isochar:   
#   [1, 0] isoscalar(1)
#   [0,1] (isovector tau_3)
#   [0.5, 0.5] proton
#   [0.5,-0.5] neutron
def isoname(isod):
    a = round(isod[0] * 2.0)
    b = round(isod[1] * 2.0)
    if a == 2 and b == 0:
        return "isoscalar"
    if a == 0 and b == 2:
        return "isovector"
    if a == 1 and b == 1:
        return "proton"
    if a == 1 and b == -1:
        return "neutron"
    raise ValueError(f"Unknown isochar value {isod}")
#
# Plot an operator vs  y in [0, 2]
# Result is influence by isochar
def plotOp(data, OpName, OpInfo, Op):
    isod = tuple(data['isochar'])
    ison = isoname(isod)
    print(f"plot {OpName} with isochar={isod} ({ison})")
    ylist = np.arange(0.0, 2.000001, 0.05)
    flist = [ Op(data, y) for y in ylist ]
    if args.v:
        print("ylist=", ylist)
        print("flist=", flist)
    fig, ax = plt.subplots()
    A = data['Isotope']
    dn = pathlib.Path(f"./{data['Element']}")
    dn.mkdir(mode=0o770, parents=True, exist_ok=True)
    if A == 0:
        ax.set_title(f"{data['Element']} Abundance Average for {ison}")
        fn = dn.as_posix() + f"/{OpName}_{ison}_Ave.pdf"
    else:
        ax.set_title(f"{data['Element']} Isotope={A} for {ison}")
        fn = dn.as_posix() + f"/{OpName}_{ison}_{A}.pdf"
    ax.set_xlabel("y")
    ax.set_ylabel(OpInfo)
    ax.plot(ylist, flist)
    plt.savefig(fn, transparent=True)
    plt.close(fig)
    

def plotspindep(data):
    plots = data['plots']
    #
    if type(plots) is str:
        # condition on type being a string so we can use lower() method
        if plots.lower() == 'all':
            plots = ['all']
        elif plots.lower() == 'none':
            return
    print("Plotting spin dependent response functions 4 Pi / (2 Ji + 1) | <Ji|O|Ji>|^2")
    isochar = data['isochar']
    if isochar == 'proton':
        isochar = [0.5, 0.5]
    elif isochar == 'neutron':
        isochar = [0.5, -0.5]
    elif isochar == 'isoscalar':
        isochar = [1.0, 0.0]
    elif isochar == 'isoscalar':
        isochar = [0.0, 1.0]
    data['isochar'] = isochar
    if len(isochar) != 2:
        raise ValueError(f"isochar value should be a 2 entry tuple, got {isod}")
    npts = 20
    print(f"isochar={isochar}")
    if ('all' in plots) or ('vcrm' in plots):
        plotOp(data, "ResM", "Vector Charge (or S.I.) Response (M)", ResM)
    if ('all' in plots) or ('alsr' in plots):
        plotOp(data, "ResSigmaPP", "Axial Longitudinal Spin Response (Sigma'')" , ResSigmaPP)
    if ('all' in plots) or ('atsr' in plots):
        plotOp(data, "ResSigmaP", "Axial Transverse Spin Response (Sigma')" , ResSigmaP)
    if ('all' in plots) or ('ssd' in plots):
        plotOp(data, "ResSD", "Standard Spin-Dependent (or S.D.) Response" , ResSD)
    if ('all' in plots) or ('vtmr' in plots):
        plotOp(data, "ResDelta", "Vector transverse magnetic response (Delta)" , ResDelta)
    if ('all' in plots) or ('vlr' in plots):
        plotOp(data, "ResPhiPP", "Vector Longitudinal Response (Phi'')" , ResPhiPP)
    if ('all' in plots) or ('vter' in plots):
        plotOp(data, "ResPhiTP", "Vector transverse electric response (PhiT')" , ResPhiTP)


# For compatability with Mathematica
# we use  cs(1..16, 0..1) and ds(1..20,0..1) and bs(1..16,0..1)
# We oversize the array and ignore 0 indices for first dim.
# Note:  Wick didn't use real arrays, he used symbols cs[idx,0], cs[idx,1] which
# is why he could use 0,1 in Mathematica for the isoscalar/isovector index
#
def expandcsdsbs(data):
    cs = np.zeros((17, 2), dtype=np.complex128) 
    ds = np.zeros((21, 2), dtype=np.complex128)
    bs = np.zeros((17, 2), dtype=np.complex128)
    if 'cs' in data:
        for ce in data['cs']:
            idx, isoscalar, isovector = (ce)
            if idx < 1 or idx > 16:
                raise ValueError(f"cs index value out of bounds 1..16, got {idx}")
            cs[idx, 0] = isoscalar
            cs[idx, 1] = isovector
        data['origcs'] = data['cs'].copy()
        data['cs'] = cs
    else:
        print("No non-relativistic LECs (cs entry) in input data")
    if 'ds' in data:
        for de in data['ds']:
            idx, isoscalar, isovector = (de)
            if idx < 1 or idx > 20:
                raise ValueError(f"ds index value out of bounds 1..20, got {idx}")
            ds[idx, 0] = isoscalar
            ds[idx, 1] = isovector
        data['origds'] = data['ds'].copy()
        data['ds'] = ds
    else:
        print("No relativistic coefficients (ds entry) in input data")
    if 'bs' in data:
        for be in data['bs']:
            idx, isoscalar, isovector = (be)
            if idx < 1 or idx > 16:
                raise ValueError(f"bs index value out of bounds 1..16, got {idx}")
            bs[idx,0] = isoscalar
            bs[idx,1] = isovector
        data['origbs'] = data['bs'].copy()
        data['bs'] = bs
    else:
        print("No relativistic LECs (bs entry, lower comp) in input data")
        data['bs'] = None

# Correct cs coefficients for muon momentum and leptonic scale mL
def relativistic_cs(data, print_details = False):
    I = 0+1.0j
    cs = data['cs']
    ds = data['ds']
    qval = data['qval']
    mL = data['mL']
    # Now compute corrected coefficients
    cs[1,0]=cs[1,0]+ds[1,0]
    cs[1,1]=cs[1,1]+ds[1,1]
    cs[10,0]=cs[10,0]+ds[2,0] * qval/(2 * mN)
    cs[10,1]=cs[10,1]+ds[2,1] * qval/(2 * mN)
    cs[11,0]=cs[11,0]-ds[3,0]
    cs[11,1]=cs[11,1]-ds[3,1]
    cs[6,0]=cs[6,0]-ds[4,0] * qval/(2 * mN)
    cs[6,1]=cs[6,1]-ds[4,1] * qval/(2 * mN)
    cs[1,0]=cs[1,0]+ds[5,0]
    cs[1,1]=cs[1,1]+ds[5,1]
    cs[2,0]=cs[2,0]+I*ds[5,0]
    cs[2,1]=cs[2,1]+I*ds[5,1]
    cs[5,0]=cs[5,0]-ds[5,0]
    cs[5,1]=cs[5,1]-ds[5,1]
    cs[4,0]=cs[4,0]-ds[5,0] * qval/(2 * mN)
    cs[4,1]=cs[4,1]-ds[5,1] * qval/(2 * mN)
    cs[6,0]=cs[6,0]-ds[5,0] * qval/(2 * mN)
    cs[6,1]=cs[6,1]-ds[5,1] * qval/(2 * mN)
    cs[4,0]=cs[4,0]+ds[6,0] * qval/mN
    cs[4,1]=cs[4,1]+ds[6,1] * qval/mN
    cs[6,0]=cs[6,0]+ds[6,0] * qval/ mN
    cs[6,1]=cs[6,1]+ds[6,1] * qval/mN
    cs[7,0]=cs[7,0]+ds[7,0]
    cs[7,1]=cs[7,1]+ds[7,1]
    cs[10,0]=cs[10,0]+I*ds[7,0]
    cs[10,1]=cs[10,1]+I*ds[7,1] 
    cs[9,0]=cs[9,0]-ds[7,0]
    cs[9,1]=cs[9,1]-ds[7,1]
    cs[10,0]=cs[10,0]-ds[8,0] * qval/mN
    cs[10,1]=cs[10,1]-ds[8,1] * qval/mN
    cs[1,0]=cs[1,0]-ds[9,0] * qval/mL
    cs[1,1]=cs[1,1]-ds[9,1] * qval/mL
    cs[5,0]=cs[5,0]-ds[9,0] * qval/mL
    cs[5,1]=cs[5,1]-ds[9,1] * qval/mL
    cs[4,0]=cs[4,0]-ds[9,0] * qval**2/(2 * mN * mL)
    cs[4,1]=cs[4,1]-ds[9,1] * qval**2/(2 * mN * mL)
    cs[6,0]=cs[6,0]-ds[9,0] * qval**2/(2 * mN * mL)
    cs[6,1]=cs[6,1]-ds[9,1] * qval**2/(2 * mN * mL)
    cs[4,0]=cs[4,0]+ds[10,0] * qval**2/(mN * mL)
    cs[4,1]=cs[4,1]+ds[10,1] * qval**2/(mN * mL)
    cs[6,0]=cs[6,0]+ds[10,0] * qval**2/(mN * mL)
    cs[6,1]=cs[6,1]+ds[10,1] * qval**2/(mN * mL)
    cs[7,0]=cs[7,0]-ds[11,0] * qval/mL
    cs[7,1]=cs[7,1]-ds[11,1] * qval/mL
    cs[9,0]=cs[9,0]-ds[11,0] * qval/mL
    cs[9,1]=cs[9,1]-ds[11,1] * qval/mL
    cs[10,0]=cs[10,0]+ds[12,0] * qval**2/(mN * mL)
    cs[10,1]=cs[10,1]+ds[12,1] * qval**2/(mN * mL)
    cs[11,0]=cs[11,0]-I*ds[13,0]
    cs[11,1]=cs[11,1]-I*ds[13,1]
    cs[8,0]=cs[8,0]-ds[13,0]
    cs[8,1]=cs[8,1]-ds[13,1]
    cs[9,0]=cs[9,0]-ds[13,0] * qval/(2 * mN)
    cs[9,1]=cs[9,1]-ds[13,1] * qval/(2 * mN)
    cs[9,0]=cs[9,0]+ds[14,0] * qval/mN
    cs[9,1]=cs[9,1]+ds[14,1] * qval/mN
    cs[14,0]=cs[14,0]-I*ds[15,0]
    cs[14,1]=cs[14,1]-I*ds[15,1]
    cs[4,0]=cs[4,0]-ds[15,0]
    cs[4,1]=cs[4,1]-ds[15,1]
    cs[6,0]=cs[6,0]+I*ds[16,0] * qval/mN
    cs[6,1]=cs[6,1]+I*ds[16,1] * qval/mN
    cs[11,0]=cs[11,0]-ds[17,0] * qval/mL
    cs[11,1]=cs[11,1]- ds[17,1] * qval/mL
    cs[8,0]=cs[8,0]-I*ds[17,0] * qval/mL
    cs[8,1]=cs[8,1]-I*ds[17,1] * qval/mL
    cs[9,0]=cs[9,0]-I*ds[17,0] * qval**2/(2 * mN * mL)
    cs[9,1]=cs[9,1]-I*ds[17,1] * qval**2/(2 * mN * mL)
    cs[16,0]=cs[16,0]-I*ds[17,0] * qval/mL
    cs[16,1]=cs[16,1]-I*ds[17,1] * qval/mL
    cs[9,0]=cs[9,0]+I*ds[18,0] * qval**2/(mN * mL)
    cs[9,1]=cs[9,1]+I*ds[18,1] * qval**2/(mN * mL)
    cs[14,0]=cs[14,0]-ds[19,0] * qval/mL
    cs[14,1]=cs[14,1]-ds[19,1] * qval/mL
    cs[4,0]=cs[4,0]-I*ds[19,0] * qval/mL
    cs[4,1]=cs[4,1]-I*ds[19,1] * qval/mL
    cs[6,0]=cs[6,0]-I*ds[19,0] * qval/mL
    cs[6,1]=cs[6,1]-I*ds[19,1] * qval/mL
    cs[6,0]=cs[6,0]+ ds[20,0] * qval**2/(mN * mL)
    cs[6,1]=cs[6,1]+ ds[20,1] * qval**2/(mN * mL)
    if print_details:
        print("Corrected cs values")
    for i in range(1, 17):
        s = cs[i,0]
        if s.imag == 0.0:
            s = s.real
        v = cs[i,1]
        if v.imag == 0.0:
            v = v.real
        if print_details:
            print(f"   i {i}  [{s}, {v}]")

# corrections to lower components
def relativistic_bs(data, print_details = False):
    I = 0+1.0j
    bs = data['bs'] # modifies in place
    ds = data['ds']
    qval = data['qval']
    mL = data['mL']
    bs[2,0]  = bs[2,0]   + I*ds[1,0]
    bs[2,1]  = bs[2,1]   + I*ds[1,1]
    bs[3,0]  = bs[3,0]   - ds[1,0]
    bs[3,1]  = bs[3,1]   - ds[1,1]
    bs[7,0]  = bs[7,0]   + I*ds[3,0]
    bs[7,1]  = bs[7,1]   + I*ds[3,1]
    bs[2,0]  = bs[2,0]   - I*ds[5,0]
    bs[2,1]  = bs[2,1]   - I*ds[5,1]
    bs[3,0]  = bs[3,0]   + ds[5,0]
    bs[3,1]  = bs[3,1]   + ds[5,1]
    bs[8,0]  = bs[8,0]   - ds[7,0]
    bs[8,1]  = bs[8,1]   - ds[7,1]
    bs[12,0] = bs[12,0] - I*ds[7,0]
    bs[12,1] = bs[12,1] - I*ds[7,1]
    bs[2,0]  = bs[2,0]   - (I*qval/mL) * ds[9,0]
    bs[2,1]  = bs[2,1]   - (I*qval/mL) * ds[9,1]
    bs[3,0]  = bs[3,0]   + (qval/mL) * ds[9,0]
    bs[3,1]  = bs[3,1]   + (qval/mL) * ds[9,1]
    bs[8,0]  = bs[8,0]   + (qval/mL) * ds[11,0]
    bs[8,1]  = bs[8,1]   + (qval/mL) * ds[11,1]
    bs[12,0] = bs[12,0] + (I*qval/mL) * ds[11,0]
    bs[12,1] = bs[12,1] + (I*qval/mL) * ds[11,1]
    bs[15,0] = bs[15,0] + (I*qval/mL) * ds[11,0]
    bs[15,1] = bs[15,1] + (I*qval/mL) * ds[11,1]
    bs[16,0] = bs[16,0] + (qval/mL) * ds[11,0]
    bs[16,1] = bs[16,1] + (qval/mL) * ds[11,1]
    bs[7,0]  = bs[7,0]   + ds[13,0]
    bs[7,1]  = bs[7,1]   + ds[13,1]
    bs[5,0]  = bs[5,0]   + ds[15,0]
    bs[5,1]  = bs[5,1]   + ds[15,1]
    bs[13,0] = bs[13,0] + I*ds[15,0]
    bs[13,1] = bs[13,1] + I*ds[15,1]
    bs[14,0] = bs[14,0] + I*ds[15,0]
    bs[14,1] = bs[14,1] + I*ds[15,1]
    bs[7,0]  = bs[7,0]   + (I*qval/mL)*ds[17,0]
    bs[7,1]  = bs[7,1]   + (I*qval/mL)*ds[17,1]
    bs[5,0]  = bs[5,0]   - (I*qval/mL)*ds[19,0]
    bs[5,1]  = bs[5,1]   - (I*qval/mL)*ds[19,1]
    bs[13,0] = bs[13,0] + (qval/mL)*ds[19,0]
    bs[13,1] = bs[13,1] + (qval/mL)*ds[19,1]
    if print_details:
        print("Corrected bs (lower comp) values")
    for i in range(1, 17):
        s = bs[i,0]
        if s.imag == 0.0:
            s = s.real
        v = bs[i,1]
        if v.imag == 0.0:
            v = v.real
        if print_details:
            print(f"   i {i}  [{s}, {v}]")


def report_decay_rate(data):
    b = data['oscb'] / hc
    qval = data['qval']
    qeff = data['qeff']
    qm = qeff / Nmass  # also nM
    y = (qeff * b / 2.0)**2
    RZ2 = data['RZ2']
    Mbar = data['Mbar']
    cs = data['cs']
    bs = data['bs']
    #
    DecayRate = 1 / (mv**4  * 2 * math.pi)
    DecayRate *= RZ2 * qeff**2 / (1 + qval / Mbar)
    h = Heff(data, qm, y, cs, bs)
    DecayRate *= h
    DecayRate /= hbar
    if abs(DecayRate.imag) > 1e-6:
        print(f" Warning:  Decay rate has imaginary component! {DecayRate}")
    DecayRate = DecayRate.real
    print(f"  Decay rate = {DecayRate:.6e}/sec")
    # save decay rate
    data['DecayRate'] = DecayRate
    if 'ExpectedDecayRate' in data:
        ed = data['ExpectedDecayRate']
        eval = abs(DecayRate - ed) / abs(ed)
        if eval > 1e-6:
            print("Failed to match expected decay rate, got {DecayRate}, expected {ed}")
            raise ValueError(f"Failed to match expected decay rate, got {DecayRate}, expected {ed}")

#
# Process one analysis input file
#
def process(path):
    with open(path, "r") as f:
        data = yaml.safe_load(f)
    ename = getelement(data)
    data['symbol'] = elements[ename]['symbol']
    print(f"Found element {ename}, symbol {data['symbol']}")
    if not ('Isotope' in data):
        raise ValueError(f"Missing Isotope specification in input - 0 for average, A for specific isotope")
    doabundance(data)
    getinteraction(data)
    getoscb(data)
    getmasses(data)
    computeThings(data)
    pdata = copy.deepcopy(data)  # save before interaction added
    readint(data)
    testFM(data)
    testWM(data)
    if errcnt > 0:
        print("Errors found during tests")
    evalsumrules(data)
    expandcsdsbs(data)
    if ('plots' in data) and ('isochar' in data):
        plotspindep(data)
    pprint.pprint(data if args.v else pdata)
    if 'mL' in data and data['mL'] > 0.0:
        if not ('ds' in data):
            raise ValueError("mL={mL} is > 0, but no ds coeffients are specified")
        relativistic_cs(data) # relativistic corrections
        if 'bs' in data:
            relativistic_bs(data) # lower component corrections
    else:
        data['mL'] = 0.0

    report_decay_rate(data)
    return data

def main():
    global errcnt
    errcnt = 0

    if args.test:
        run_tests()
        return
    loadElements()
    setElasticDir(args)
    if args.pelements:
        run_pelements()
        return
    if len(args.path) == 0:
        print("Expecting input files to process")
        return
    for p in args.path:
        print(f"Processing {p}")
        process(p)

# For MuonConverter
def compute_decay_rate(nr_model_path, elastic_path = './Elastic/'):
    """
    Function for computing decay rate
    """
    # Open and read yaml file contianing model info
    with open(nr_model_path, "r") as f:
        data = yaml.safe_load(f)

    # Load in element list
    with open("elements.yaml", "r") as f:
        global elements
        elements = yaml.safe_load(f)
    for name, eht in elements.items():
        Z = eht['Z']
        isotopeht = eht['isotopes']
        idx = 0
        for A, pht in isotopeht.items(): 
            pht['N'] = A - Z
            pht['mass'] = auMeV * A + pht['delta']
            pht['idx'] = idx
            idx = idx + 1

    # Get the name of the element
    if 'Element' in data:
        ename = data['Element']
    else:
        raise ValueError("Missing Element tag in input")
    if not ename in elements:
        # Correct for abbreviated name
        for e, edata in elements.items():
            if edata['symbol'] == ename or (('alias' in edata) and (edata['alias'] == ename)):
                ename = e # convert
                break
        if not ename in elements:
            raise ValueError(f"unknown element name or symbol {ename}")
        data['Element'] = ename
        
    # Set symbol name
    data['symbol'] = elements[ename]['symbol']
        
    # Set elastic directory
    edir = pathlib.Path(elastic_path).expanduser()
    if not edir.is_dir():
        raise ValueError(f"Supplied Elastic density matrix directory {args.edir} is not a directory")
    global elasticDir
    elasticDir = edir

    # Compute
    doabundance(data)
    getinteraction(data)
    getoscb(data)
    getmasses(data)
    computeThings(data)
    pdata = copy.deepcopy(data)  # save before interaction added
    readint(data)
    evalsumrules(data)
    expandcsdsbs(data)
    
    if ('plots' in data) and ('isochar' in data):
        plotspindep(data)
        
    # Check for mL value
    if 'mL' in data and data['mL'] > 0.0:
        if not ('ds' in data):
            raise ValueError("mL={mL} is > 0, but no ds coeffients are specified")
        relativistic_cs(data) # relativistic corrections
        if 'bs' in data:
            relativistic_bs(data) # lower component corrections
    else:
        data['mL'] = 0.0

    b = data['oscb'] / hc
    qval = data['qval']
    qeff = data['qeff']
    qm = qeff / Nmass  # also nM
    y = (qeff * b / 2.0)**2
    RZ2 = data['RZ2']
    Mbar = data['Mbar']
    cs = data['cs']
    bs = data['bs']
    
    DecayRate = 1 / (mv**4  * 2 * math.pi)
    DecayRate *= RZ2 * qeff**2 / (1 + qval / Mbar)
    h = Heff(data, qm, y, cs, bs)
    DecayRate *= h
    DecayRate /= hbar
    if abs(DecayRate.imag) > 1e-6:
        print(f" Warning:  Decay rate has imaginary component! {DecayRate}")
    DecayRate = DecayRate.real

    return DecayRate 
    
    

if __name__ == '__main__':
    main()