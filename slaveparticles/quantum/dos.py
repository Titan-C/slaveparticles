# -*- coding: utf-8 -*-
"""
@author: Óscar Nájera
Created on Wed Apr  2 16:49:03 2014
"""
from __future__ import division, absolute_import, print_function
import numpy as np
from slaveparticles.quantum.operators import fermi_dist

#Bethe Lattice
def bethe_lattice(energy, hopping):
    """Bethe lattice in inf dim density of states"""
    energy = np.asarray(energy).clip(-2*hopping, 2*hopping)
    return np.sqrt(4*hopping**2 - energy**2) / (2*np.pi*hopping**2)

def bethe_fermi(energy, quasipart, shift, hopping, beta):
    """product of the bethe lattice dos, fermi distribution"""
    return fermi_dist(quasipart * energy - shift, beta) \
            * bethe_lattice(energy, hopping)

def bethe_fermi_ene(energy, quasipart, shift, hopping, beta):
    """product of the bethe lattice dos, fermi distribution an weighted
       by energy"""
    return energy * bethe_fermi(energy, quasipart, shift, hopping, beta)


def bethe_filling_zeroT(fermi_energy, hopping):
    """Returns the particle average count given a certan fermi energy, for the
       semicircular density of states of the bethe lattice"""
    fermi_energy = np.asarray(fermi_energy).clip(-2*hopping, 2*hopping)
    return 1/2. + fermi_energy/2 * bethe_lattice(fermi_energy, hopping) \
                    + np.arcsin(fermi_energy/2/hopping)/np.pi

from scipy.optimize import fsolve
def bethe_findfill_zeroT(particles, orbital_e, hopping):
    """Return the fermi energy that correspond to the given particle quantity
       in a semicircular density of states of a bethe lattice in a multi
       orbital case that can be non-degenerate"""

    assert 0. <= particles <= len(orbital_e)

    zero = lambda e: np.sum([bethe_filling_zeroT(e-e_m, t) \
                        for t, e_m in zip(hopping, orbital_e)]) - particles
    return fsolve(zero, 0)

def bethe_find_crystalfield(populations, hopping):
    """Return the orbital energies to have the system populates as
       desired by the given individual populations"""

    zero = lambda orb: [bethe_filling_zeroT(-em, tz) - pop \
                        for em, tz, pop in zip(orb, hopping, populations)]

    return fsolve(zero, np.zeros(len(populations)))

def bethe_ekin_zeroT(fermi_energy, hopping):
    """Returns the kinetic energy of the system up to the fermi energy
       at zero temperature but when the quasiparticle weight is not zera"""
    return bethe_lattice(fermi_energy, hopping)*(fermi_energy**2 -4*hopping**2)/3.
