# -*- coding: utf-8 -*-
"""
@author: Óscar Nájera
Created on Fri Apr 25 09:44:11 2014
Tests over the energy integrals
"""
from __future__ import division, absolute_import, print_function
from scipy.integrate import quad
import numpy as np
import random
import slaveparticles.spins.fermion as ss
import slaveparticles.quantum.dos as dos

hop = random.uniform(0.5, 2)

def test_filling_integral():
    """Tests that the integration over the semicircular density of states
       of a bethe lattices gives the appropiate particle density"""

    assert abs(dos.bethe_filling_zeroT(2*hop, hop) - 1.) < 1e-15
    assert abs(dos.bethe_filling_zeroT(0, hop)   - 1/2.) < 1e-15
    assert abs(dos.bethe_filling_zeroT(-3*hop, hop))     < 1e-15
    assert abs(dos.bethe_filling_zeroT(-hop, hop) -
              quad(dos.bethe_fermi, -2*hop, -hop, (1, 0, hop, 1e5))[0]) < 1e-15
    assert (abs(ss.fermion_avg([2*hop, 0], [hop]*2, 'ocupation')
    - np.array([1, 1/2.])) < 1e-15).all()


def test_find_fermi_energy_zeroT():
    """Test over the know cases of which fermi energy comes from a certain
       particle density"""

    assert abs(dos.bethe_findfill_zeroT(1., [0],  [hop]) - 2*hop) < 5e-8
    assert abs(dos.bethe_findfill_zeroT(1/2, [0], [hop]) - 0    ) < 5e-8
    assert -2*hop < dos.bethe_findfill_zeroT(1/4., [0], [hop]) < 0

    #multiorbitals
    assert abs(dos.bethe_findfill_zeroT(3., [0, 2*hop]*2, [hop]*4) - 2*hop) < 6e-8


def test_kinetic_energy_zeroT():
    """Test over the kinetic energy integrals"""
    assert 0 <= dos.bethe_ekin_zeroT(2*hop, hop) < 5e-15
    assert abs(dos.bethe_ekin_zeroT(0, hop) + 4*hop/3/np.pi) < 5e-15

def test_dos():
    """test over dos. functions"""
    assert dos.bethe_lattice(-2, 0.5) == 0.
    assert abs(quad(dos.bethe_lattice, -2, 2, args=(1))[0] -1.0) < 5e-16
    assert abs(quad(dos.bethe_fermi_ene, -2, 2, \
                    args=(1, 0, 1, 1e5))[0] +4/3/np.pi) < 5e-15
