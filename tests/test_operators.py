# -*- coding: utf-8 -*-
"""
@author: Óscar Nájera
Created on Wed Jun  4 14:37:13 2014
"""
from __future__ import division, absolute_import, print_function
from slaveparticles.quantum import fermion, operators
from scipy.sparse import eye
import pytest
import numpy as np


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


@pytest.mark.parametrize("U, mu, beta", [(1, 0, 50), (2, -0.5, 500)])
def test_gf_lehmann(U, mu, beta):
    """Verifies the lehmann representation of the atom"""
    w = np.linspace(-1.5*U, 1.5*U, 500)
    d_up, d_dw = [fermion.destruct(2, sigma) for sigma in range(2)]
    sigma_z = d_up.T*d_up - d_dw.T*d_dw
    H = - U/2 * sigma_z * sigma_z - mu * (d_up.T*d_up + d_dw.T*d_dw)

    e, v = operators.diagonalize(H.todense())
    g_up = operators.gf_lehmann(e, v, d_up.T, beta, w)
    g_up_cross = operators.gf_lehmann(e, v, d_up.T, beta, w, d_up)

    Z = 1+2*np.exp(beta*(U/2+mu)) + np.exp(2*beta*mu)
    g_up_ref = ((1+np.exp(beta*(U/2+mu)))/(w + mu + U/2.) +
                (np.exp(beta*(U/2+mu)) + np.exp(2*beta*mu))/(w + mu - U/2.))/Z

    assert np.allclose(g_up, g_up_ref)
    assert np.allclose(g_up_cross, g_up_ref)

def test_gf_lehmann2():
    """Verifies the lehmann representation of a random Hamiltonian"""
    w = np.linspace(-1.5, 1.5, 500)
    d_up = fermion.destruct(6, 0)
    H = np.random.rand(*d_up.shape)
    H += H.T

    e, v = operators.diagonalize(H)
    g_up_diag = operators.gf_lehmann(e, v, d_up.T, 400, w)
    g_up_cross = operators.gf_lehmann(e, v, d_up.T, 400, w, d_up)
    assert np.allclose(g_up_diag, g_up_cross)

if __name__ == "__main__":
    pass
