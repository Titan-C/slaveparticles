# -*- coding: utf-8 -*-
"""
@author: Óscar Najera
Created on Fri Apr 18 11:09:05 2014

Compares to results of the multiorbital degenerate single site half-filled
slave spin system
"""

import slaveparticles.spins as ss
import numpy as np
import pytest

def system_output(slsp, coulomb, results, tol = 3e-6, lone_con = False, J_coup=0.0):
    """runs a simulation on the given system setup and compares to know
        output"""

    if lone_con:
        for inter, ref in zip(coulomb, results):
            slsp.selfconsistency(inter, J_coup)
            assert (np.abs(slsp.quasiparticle_weight() - ref) < tol).all()

    #bug on solver changing initial guess
    field = -1*np.ones((2, slsp.param['slaves']))
    for inter, ref in zip(coulomb, results):
        slsp.selfconsistency(inter, J_coup, field)
        field = slsp.mean_field()
        assert (np.abs(slsp.quasiparticle_weight() - ref) < tol).all()

def test_singleorbit():
    """test the satisfactory results for the single orbital single
       site mean field hamiltonian"""
    slsp = ss.Spinon()
    coulomb = [ 0., 1.2, 2.8, 3.35, 3.4 ]
    results = [ 1., 0.87508782, 0.31992258, 0.02650960, 0.  ]

    system_output(slsp, coulomb, results, True)

test_singleorbit()

@pytest.mark.slow
def test_two_orbit():
    """test the satisfactory results for the two orbitals single
       site mean field hamiltonian"""
    coulomb = [ 0., 1.7, 3.5, 5.02, 5.1]
    slsp = ss.Spinon(slaves = 4, orbitals=2, avg_particles=2, \
                            hopping = [0.5]*4, populations=[0.5]*4)
    results = [ 1., 0.78889506, 0.44144503, 0.02329075, 0.  ]

    system_output(slsp, coulomb, results)

@pytest.mark.slow
def test_3_orbit():
    """test the satisfactory results for the three orbitals single
       site mean field hamiltonian"""
    coulomb = [ 0.  ,  2.7,  4.5,  6.74, 6.8]
    slsp = ss.Spinon(slaves=6, orbitals=3, avg_particles=3, \
                            hopping = [0.5]*6, populations=[0.5]*6)
    results = [ 1., 0.68270056, 0.43086564, 0.0110575, 0. ]

    system_output(slsp, coulomb, results)

def test_singleorbit_dopping():
    """test quasiparticle behavior on hole dopped system single orbit"""

    coulomb = [ 0., 1.2, 2.8, 3.35, 4.0 ]
    results = [[ 1., 0.87513694,  0.32415401,  0.09952643,  0.04759491],
               [ 1., 0.87528405,  0.33550111,  0.14974231,  0.08927663],
               [ 1., 0.8798186 ,  0.47724379,  0.38423417,  0.32283685],
               [ 1., 0.9422046 ,  0.8442584 ,  0.82157692,  0.80050541],
               [ 1., 0.99183007,  0.98036478,  0.97748548,  0.97462334]]

    for den, zet in  zip([0.495, 0.49, 0.45, 0.25, 0.05], results):
        slsp = ss.Spinon(populations= [den]*2)
        system_output(slsp, coulomb, zet)