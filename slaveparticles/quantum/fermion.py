# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 11:33:36 2014

@author: oscar
"""
from __future__ import division, absolute_import, print_function
import numpy as np
from slaveparticles.quantum.spinmatrices import btest
from scipy.sparse import csr_matrix


def phase(state, index):
    """Counts the number of fermions present in the state before the indexed
       one. Then returns the fermion phase sign"""
    amount = bin(state >> index+1).count('1')
    if amount % 2 == 0:
        return 1
    else:
        return -1


def destruct(particles, index):
    """Fermion annihilation operator in matrix representation for a indexed
       particle in a bounded N-particles fermion fock space"""

    mat = np.zeros((2**particles, 2**particles))

    flipper = 2**index
    for i in range(2**particles):
        ispin = btest(i, index)
        if ispin == 1:
            mat[i ^ flipper, i] = phase(i, index)
    return csr_matrix(mat)