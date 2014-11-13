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
    return first.dot(second) - second.dot(first)

def anticommutator(first, second):
    """Applies the anti-commutation of two matrix operators"""
    return first.dot(second) + second.dot(first)

def fermion_phase(state, index):
    """Counts the number of fermions present in the state before the indexed
       one. Then returns the fermion phase sign"""
    amount = bin(state >> index+1).count('1')
    if amount % 2 == 0:
        return 1
    else:
        return -1

from scipy.sparse import csr_matrix
def f_destruct(particles, index):
    """Fermion annihilation operator in matrix representation for a indexed
       particle in a bounded N-particles fermion fock space"""

    mat = np.zeros((2**particles, 2**particles))

    flipper = 2**index
    for i in range(2**particles):
        ispin = btest(i, index)
        if ispin == 1:
            mat[i ^ flipper, i] = fermion_phase(i, index)
    return csr_matrix(mat)

def f_creation(particles, index):
    """Fermion creation operator in matrix representation for an indexed
       particle in a bounded N-particles fermion fock space"""
    return f_destruct(particles, index).T


#Tools
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
