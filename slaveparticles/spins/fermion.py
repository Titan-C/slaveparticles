# -*- coding: utf-8 -*-
"""
@author: Óscar Nájera
Created on Wed Jun 18 09:19:33 2014
"""
from __future__ import division, absolute_import, print_function
from slaveparticles.quantum.dos import bethe_ekin_zeroT, \
    bethe_filling_zeroT
import numpy as np

def orbital_energies(param, quasiparticle):
    """calculates the output orbital energies for the system after the spin
       problem has converged for the lagrange multipliers and knows the state
       of the system for the given populations"""

    return param['lambda'] + quasiparticle * param['orbital_e_free']

def fermion_avg(efermi, norm_hopping, func):
    """calcules for every slave it's average over the desired observable"""
    if func == 'ekin':
        func = bethe_ekin_zeroT
    elif func == 'ocupation':
        func = bethe_filling_zeroT

    return np.asarray([func(ef, tz) for ef, tz in \
            zip(efermi, norm_hopping)])
