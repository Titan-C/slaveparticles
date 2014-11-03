# -*- coding: utf-8 -*-
"""
@author: Óscar Najera
Created on Fri Apr 18 11:09:05 2014

Compares to results of the multiorbital degenerate single site half-filled
slave spin system
"""

import slaveparticles.spins as ss
import numpy as np
import scipy.linalg as LA

def test_Hamitonian_hermitian():
    """Test if the hamiltonian is hermitian, after bug of transposing operators"""
    sl = ss.Spinon()
    assert (sl.H_s - sl.H_s.T < 5e-14).all()

def test_diagonalization():
    """tests if diagonalization states give identity matrix"""
    U = 6.0
    sl = ss.Spinon(slaves=6, orbitals=3, avg_particles=3,
                   hopping=[1]*6, populations=[0.5]*6)
    assert (sl.eig_energies >= 0.).all()
    h = sl.mean_field()
    assert (abs(h) <= 1.0).all()
    if (abs(h) <= 1.0).all():
        sl.oper['Hint'] = sl.inter_spin_hamiltonian(U, 0.)
        sl.update_H(-0.5*np.ones((2, 6)), np.zeros(6))
        assert LA.norm(np.dot(sl.eig_states.T, sl.eig_states)-np.eye(64)) < 5e-14


def test_averager():
    """tests expected value against analitiacally calculated one"""

    sl = ss.Spinon(hopping=[1]*2)

    sl.update_H(-0.5*np.ones((2, 2)), np.zeros(2))
    assert abs(sl.expected(sl.H_s)+2.) < 5e-16

    for beta in [2., 60.0, 800.]:
        sl.selfconsistency(1, 0, np.zeros((2, 2)))
        assert abs(sl.expected(sl.H_s, beta) - \
                    np.exp(-beta/2.)/(2+2*np.exp(-beta/2.))) < 5e-13

        for operator in sl.oper['O']:
            assert abs(sl.expected(operator, beta)) < 5e-15
            assert abs(sl.expected(np.dot(operator, operator), beta)- 1) < 5e-13

#    sl = ss.Spinon(slaves = 8, orbitals=4, avg_particles=4, hopping=[1]*4, orbital_e=[0]*4)
#    sl.oper['Hint'] = sl.inter_spin_hamiltonian(1., 0.)
#    sl.update_H(-0.5*np.ones((2,8)), np.zeros(8))
#    assert LA.norm(sl.H_s -sl.H_s.T.conjugate()) < 5e-16


def test_equivalent_s():
    """tests if the expected value of all spin matrices is equivalent"""
    sl = ss.Spinon(slaves=6, orbitals=3, avg_particles=3,
                   hopping=[1]*6, populations=[0.5]*6)
    sl.selfconsistency(1, 0, np.zeros((2, 6)))

    spi = sl.expected(sl.oper['O'][0])
    for s in sl.oper['O']:
        assert abs(sl.expected(np.dot(s, s))-1.) < tol
        avsx = sl.expected(s)
        assert abs(avsx-spi) < tol
        spi = avsx

    for s in sl.oper['Sz']:
        assert abs(sl.expected(np.dot(s, s))-1/4.) < tol
        assert abs(sl.expected(s)) < tol

def test_quasiparticle():
    """tests for the quasiparticle weight"""
    sl = ss.Spinon(hopping=[1]*2)
    sl.selfconsistency(0, 0, np.zeros((2, 2)))
    assert (sl.quasiparticle_weight() -1. < 1e-13).all()

tol = 3e-5
