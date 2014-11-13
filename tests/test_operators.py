# -*- coding: utf-8 -*-
"""
@author: Óscar Nájera
Created on Wed Jun  4 14:37:13 2014
"""
from __future__ import division, absolute_import, print_function
from slaveparticles.quantum import operators as qo
from scipy.sparse import eye

def test_fermion_anticommut():
    """Verifies the fermion anticommutation relations"""
    particles = 3
    c_dagger = [qo.f_creation(particles, index) for index in range(particles)]
    c_destru = [qo.f_destruct(particles, index) for index in range(particles)]
    
    for i in range(particles):
        for j in range(i, particles):
            assert qo.anticommutator(c_dagger[i], c_dagger[j]).nnz == 0
            assert qo.anticommutator(c_destru[i], c_destru[j]).nnz == 0

            k_delta = qo.anticommutator(c_destru[i], c_dagger[j])
            if i == j:
                assert (k_delta - eye(2**particles)).nnz == 0
            else:
                assert k_delta.nnz == 0
    

if __name__ == "__main__":
    pass
