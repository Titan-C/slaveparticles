# -*- coding: utf-8 -*-
"""
===============================================
Effect of crystalfield field degeneracy lifting
===============================================

@author: Óscar Nájera
Created on Wed Sep 24 14:55:15 2014
"""
from __future__ import division, absolute_import, print_function
import numpy as np
from slavespins.spinon import Spinon
import slavespins.plotter as ssplt

def crystalfield(interaction=np.linspace(0, 20, 201), \
                 j_hund=np.linspace(0, 0.35, 71)):
    """Aimed at reproducing the figure in paper
       L. de'Medici, PRB 83,205112 (2011)
       showing the phase diagram of a 3 band hubbard with one lifted band
       fixed population 1:1.5,1.5"""
    slsp = Spinon(slaves=6, orbitals=3, hopping=[0.5]*6, \
               populations=[1, 1, 1.5, 1.5, 1.5, 1.5])

    zet = []
    for hund_cu in j_hund:
        zet.append(ssplt.solve_loop(slsp, interaction, hund_cu)[0][0])
    np.savez('PRB_83_205112', zeta=zet, u_int=interaction, j_hund=j_hund)

if __name__ == "__main__":
    crystalfield()
