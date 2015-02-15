# -*- coding: utf-8 -*-
"""
==========================================
Reconstructing a Coulomb occupation ladder
==========================================

Filling of the 2 degenerate orbitals of an atom in function of the chemical
potential
"""

from scipy.special import binom
from scipy.optimize import fsolve
from numpy import linspace, exp, sum, zeros
from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, tight_layout

#exact case
def spectrum(mu, orbitals, U, beta, Q):
    return binom(2*orbitals, Q) * exp(-beta *(U/2.*(Q - orbitals)**2 - mu*Q))

def expected_filling(mu, orbitals, U, beta):
    Z = sum([spectrum(mu, orbitals, U, beta, Q) for Q in range(2*orbitals+1)], axis=0)
    n_avg = sum([Q*spectrum(mu, orbitals, U, beta, Q) for Q in range(2*orbitals+1)], axis=0)
    return n_avg / Z

def fermi_dist(energy, beta):
    """ Fermi dirac distribution"""
    return 1./(exp(beta*energy) +1)

def restriction(lam, mu, orbitals, U, beta):
    """Equation that determines the restriction on lagrange multipier"""
    return 2*orbitals*fermi_dist(-(mu + lam), beta) - expected_filling(-1*lam, orbitals, U, beta)

def main(orbitals, beta, U, step):
    mu = linspace(-U*orbitals, U*orbitals, step)
    lam = fsolve(restriction, -mu, (mu, orbitals, U, beta))
    plot(mu, expected_filling(mu, orbitals, U, beta), '--', label='Exact')
    plot(mu, 2*orbitals*fermi_dist(-(mu+lam), beta), label='Slave spin approx')

    legend(loc=0)
    title('Orbitals ocupation, $\\beta = {} $, $U= {} $'.format(beta, U), fontsize=14)
    xlabel('$\mu$', fontsize=20)
    ylabel('$n$', fontsize=20)
    tight_layout()
    return mu, lam

def pressision_try(orbitals, U, beta, step):
    """perform a better initial guess of lambda
       no improvement"""
    mu, lam = main(orbitals, U, beta, step)
    mu2, lam2 = linspace(0, U*orbitals, step), zeros(step)
    for i in range(99):
        lam2[i+1] = fsolve(restriction, lam2[i], (mu2[i+1], orbitals, U, beta))
    plot(mu2, 2*orbitals*fermi_dist(-(mu2+lam2), beta), label='Test guess')
    legend(loc=0)

if __name__ == "gallery":
    mu, lam = main(2, 50, 2, 200)
