#!/usr/bin/python
# -*- coding: utf-8 -*-

from math import factorial, pi
from constants import *
import sys

# A python script to calculate the analytic potential in a 
# spherical dielectric cavity with off-centre charge.

kappa = float(sys.argv[1]) #0.125713350686
q = 1.0
D1 = 80.0 # internal dielectric
D2 = 80.0 # external dielectric
dielectric_ratio = D2 / D1
#Rc = 0.90 # eccentricity
b = 1.0 # radius of the cavity
r1 = 0.0 # radius of the charge
x = kappa*b
max_n = 10000 # number of terms to add

Rc = 0.0
delta = 0.01
#while (Rc < 1.0):
for Rc in [0.0]:

    Y = Rc / b

    def calc_K(n):
        # not actually used since the recurrence relations
        # are more stable/useful

        nfac = factorial(n)
        two_nfac = factorial(2*n)
        nfac_ratio = float(nfac) / float(two_nfac)
        
        nums = []
        for s in range(0,n+1):
            
            sfac = factorial(s)
            a = 2.0**s * nfac_ratio / sfac
            b = factorial(2*n - s) / factorial(n-s)
            nums.append(a*b*x**s)
            
        return sum(nums)

    def calc_delta_K(n, Kn_minus):
        # recurrence relation for K_n+1 given K_n and K_n-1
        return x*x*Kn_minus / ((2*n+1)*(2*n-1))

    # first two values are easy to calculate by hand
    Kn_vals = [1.0, 1.0+x]

    def Kn(n):
        # returns the Kn value (and caches the result too)
        while (n > len(Kn_vals) - 1):
            last_n = len(Kn_vals) - 1
            next_n = last_n + 1
            Kn_vals.append(Kn_vals[-1] + calc_delta_K(last_n, Kn_vals[-2]))
        return Kn_vals[n]

    def K_ratio(n):
        if (x == 0): return 0
        return (-(2*n + 1)*(Kn(n+1) / Kn(n)) + (2*n + 1 + x))/x

    def q_n_x(n):
        return 1.0 - K_ratio(n)

    def Y_n(n):
        return Y**n;

    def B_n_ij(n):
        val = (Y_n(n)*Y_n(n) / b) * ((n+1)*(1.0 - dielectric_ratio) / ((n+1)*dielectric_ratio + n))
        return val

    def C_n_ij(n):
        val = -(Y_n(n)*Y_n(n) / b)*(((2*n+1)*dielectric_ratio*x*q_n_x(n))/ ( ((n+1)*dielectric_ratio + n)*((n+1)*dielectric_ratio + n + dielectric_ratio*x*q_n_x(n)) ))
        return val

    series = [B_n_ij(n) + C_n_ij(n) for n in range(0, max_n)]

    potential = q * sum(series) / (2.0*4.0*pi)
    #print series

    #potential_conversion_factor = elementary_charge * elementary_charge * kJ_in_kcal * Avogadro / (epsilon0*D1*Angstroms*1000.)
    #force_conversion_factor = 100* elementary_charge * elementary_charge * kJ_in_kcal * Avogadro/ (epsilon0*D1*Angstroms*1000.)
    potential_conversion_factor = elementary_charge * elementary_charge * Avogadro / (epsilon0*D1*Angstroms*1000.)
    #force_conversion_factor = 1e12 * elementary_charge * elementary_charge / (epsilon0*D1*Angstroms*Angstroms)
    force_conversion_factor = elementary_charge * elementary_charge * Avogadro / (epsilon0*D1*Angstroms*1000.)

    #print "energy (kJ/mol): ", potential * potential_conversion_factor
    if Y > 0:
        force_series = [n*B_C/Y for n,B_C in enumerate(series)]
        force = -q*q*sum(force_series) / (4.0*pi*b)
    else:
        force = 0
    #print "radial qE force (kJ/molA): ", force * force_conversion_factor
    print(Rc,potential*potential_conversion_factor,force*force_conversion_factor)

    Rc += delta
    
        
    
