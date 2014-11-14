# -*- coding: utf-8 -*-
"""
@author: Óscar Nájera
Created on Wed Apr  2 16:49:03 2014
"""
from __future__ import division, absolute_import, print_function
import numpy as np
import scipy.linalg as LA
from slaveparticles.quantum.spinmatrices import *


def fermi_dist(energy, beta):
    """ Fermi Dirac distribution"""
    return 1./(np.exp(beta*energy) + 1)


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


def gf_lehmann(eig_e, eig_states, d_dag, beta, omega):
    """Outputs the lehmann representation of the greens function
       omega has to be given, as matsubara frequencies or their analitical
       continuation"""
    zet = partition_func(beta, eig_e)
    G = 0
    for i in range(len(eig_e)):
        for j in range(len(eig_e)):
            G += np.dot(eig_states[:, j].T, d_dag.dot(eig_states[:, i]))**2 * \
                       (np.exp(-beta*eig_e[i]) + np.exp(-beta*eig_e[j])) / \
                       (omega + eig_e[i] - eig_e[j])
    return G / zet


def expected_value(operator, eig_values, eig_states, beta):
    """Calculates the average value of an observable
       it requires that states and operators have the same base"""
    aux = np.einsum('i,ij,ji', np.exp(-beta*eig_values),
                    eig_states.T, np.dot(operator, eig_states))

    return aux / partition_func(beta, eig_values)
