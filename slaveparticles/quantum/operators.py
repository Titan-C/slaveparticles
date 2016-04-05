# -*- coding: utf-8 -*-
"""
@author: Óscar Nájera
Created on Wed Apr  2 16:49:03 2014
"""
from __future__ import division, absolute_import, print_function
import numpy as np
import scipy.linalg as LA
from itertools import product


def fermi_dist(energy, beta):
    """ Fermi Dirac distribution"""
    exponent = np.asarray(beta*energy).clip(-600, 600)
    return 1./(np.exp(exponent) + 1)


def partition_func(beta, energies):
    """Evaluate the partition function for an array of energies"""
    return np.exp(-beta*energies).sum()


# Operators
def commutator(first, second):
    """Applies the commutation of two matrix operators"""
    return first.dot(second) - second.dot(first)


def anticommutator(first, second):
    """Applies the anti-commutation of two matrix operators"""
    return first.dot(second) + second.dot(first)


# Tools
def diagonalize(operator):
    """diagonalizes single site Spin Hamiltonian"""
    eig_values, eig_vecs = LA.eigh(operator)
    emin = np.amin(eig_values)
    eig_values -= emin

    return eig_values, eig_vecs


def gf_lehmann(eig_e, eig_states, d_dag, beta, omega, d=None):
    """Outputs the lehmann representation of the greens function
       omega has to be given, as matsubara or real frequencies"""

    ew = np.exp(-beta*eig_e)
    zet = ew.sum()
    G = np.zeros_like(omega)
    basis_create = np.dot(eig_states.T, d_dag.dot(eig_states))
    if d is None:
        tmat = np.square(basis_create)
    else:
        tmat = np.dot(eig_states.T, d.T.dot(eig_states))*basis_create

    tmat *= np.add.outer(ew, ew)
    gap = np.add.outer(-eig_e, eig_e)

    N = eig_e.size
    for i, j in product(range(N), range(N)):
        G += tmat[i, j] / (omega + gap[i, j])
    return G / zet


def expected_value(operator, eig_values, eig_states, beta):
    """Calculates the average value of an observable
       it requires that states and operators have the same base"""
    aux = np.einsum('i,ji,ji', np.exp(-beta*eig_values),
                    eig_states, operator.dot(eig_states))

    return aux / partition_func(beta, eig_values)
