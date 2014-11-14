# -*- coding: utf-8 -*-
"""
@author: Óscar Nájera
Created on Wed Jun  4 14:37:13 2014
"""
from __future__ import division, absolute_import, print_function
from slaveparticles.quantum import fermion, operators
from scipy.sparse import eye


def test_fermion_anticommut():
    """Verifies the fermion anticommutation relations"""
    particles = 3
    c = [fermion.destruct(particles, index) for index in range(particles)]

    for i in range(particles):
        for j in range(i, particles):
            assert operators.anticommutator(c[i].T, c[j].T).nnz == 0
            assert operators.anticommutator(c[i], c[j]).nnz == 0

            k_delta = operators.anticommutator(c[i], c[j].T)
            if i == j:
                assert (k_delta - eye(2**particles)).nnz == 0
            else:
                assert k_delta.nnz == 0


if __name__ == "__main__":
    pass
