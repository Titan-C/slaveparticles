# -*- coding: utf-8 -*-
"""
======================================================
Drop of quasiparticle weight by increasing interaction
======================================================

The quasiparticle weight of the electronic system drops as the local interaction
is increased. When studying the two band case, the half-filled scenario and
other commensurate fillings
undergoe a transition into the Mott insulator. All other doping configurations
retain the strongly correlated metal.
"""

from __future__ import division, absolute_import, print_function
import slaveparticles.utils.plotter as ssplt
import numpy as np
import matplotlib.pyplot as plt

#band dop_phasediag
def plot_dop_phase(bands, int_max, hund_cu):
    """Phase plot of Quasiparticle weight for N degenerate bands
    under doping shows transition only at integer filling
    the rest are metallic states"""
    name = 'Z_dop_phase_'+str(bands)+'bands_U'+str(int_max)+'J'+str(hund_cu)
    dop = np.sort(np.hstack((np.linspace(0.01,0.99,50),
                    np.arange(1./2./bands, 1, 1/2/bands))))
    data = ssplt.calc_z(bands, dop, np.arange(0, int_max, 0.1), hund_cu, name)

    ssplt.surf_z(data, name)

plot_dop_phase(2, 6, 0.)

