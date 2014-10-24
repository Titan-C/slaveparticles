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
    return 1./(np.exp(beta*energy) +1)

def partition_func(beta, energies):
    """Evaluate the partition function for an array of energies"""
    return np.exp(-beta*energies).sum()

#Operators
def commutator(first, second):
    """Applies the commutation of two matrix operators"""
    return np.dot(first, second) - np.dot(second, first)

def anticommutator(first, second):
    """Applies the anti-commutation of two matrix operators"""
    return np.dot(first, second) + np.dot(second, first)

def diagonalize(operator):
    """diagonalizes single site Spin Hamiltonian"""
    eig_values, eig_vecs = LA.eigh(operator)
    emin = np.amin(eig_values)
    eig_values -= emin

    return eig_values, eig_vecs

def expected_value(operator, eig_values, eig_states, beta):
    """Calculates the average value of an observable
       it requires that states and operators have the same base"""
    aux = np.einsum('i,ij,ji', np.exp(-beta*eig_values),
                    eig_states.T, np.dot(operator, eig_states))

    return aux / partition_func(beta, eig_values)
