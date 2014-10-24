# -*- coding: utf-8 -*-
"""
@author: Óscar Nájera
Created on Thu Apr  3 18:47:53 2014
"""
from __future__ import division, absolute_import, print_function
import numpy as np

def btest(state, index):
    """A bit test that evaluates if 'state' in binany has a one in the
       'index' location. returns one if true"""
    return (state >> index) & 1

def spin_z(particles, index):
    """Generates the spin_z projection operator for a system of
        N=particles and for the selected spin index name. where index=0..N-1"""
    mat = np.zeros((2**particles, 2**particles))

    for i in range(2**particles):
        ispin = btest(i, index)
        if ispin == 1:
            mat[i, i] = 1
        else:
            mat[i, i] = -1
    return 1/2.*mat

def spin_gen(particles, index, gauge=1):
    """Generates the generic spin operator in z basis for a system of
        N=particles and for the selected spin index name. where index=0..N-1
        The gauge term sets the behavoir for a system away from half-filling"""
    mat = np.zeros((2**particles, 2**particles))

    flipper = 2**index
    for i in range(2**particles):
        ispin = btest(i, index)
        if ispin == 1:
            mat[i ^ flipper, i] = 1
        else:
            mat[i ^ flipper, i] = gauge
    return mat

def spin_x(particles, index):
    """Generates the spin_x projection operator in z basis for a system of
        N=particles and for the selected spin index name. where index=0..N-1"""

    return 1/2.*spin_gen(particles, index)

def spin_y(particles, index):
    """Generates the spin_y projection operator in z basis for a system of
        N=particles and for the selected spin index name. where index=0..N-1"""

    return 1j/2.*spin_gen(particles, index, -1)
