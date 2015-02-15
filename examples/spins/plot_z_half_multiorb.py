# -*- coding: utf-8 -*-
"""
======================================================
Drop of quasiparticle weight by increasing interaction
======================================================

The quasiparticle weight of the electronic system drops as the local interaction
is increased. Multi orbital degenerate systems are studied
"""

from __future__ import division, absolute_import, print_function
import slaveparticles.utils.plotter as ssplt
import numpy as np
import matplotlib.pyplot as plt

#Degenerate bands
def plot_degbandshalffill():
    """Plot of Quasiparticle weight for degenerate
    half-filled bands, showing the Mott transition"""
    ulim = [3.45, 5.15, 6.85, 8.55]
    bands = range(1, 5)
    for band, u_int in zip(bands, ulim):
        name = 'Z_half_'+str(band)+'band'
        dop = [0.5]
        data = ssplt.calc_z(band, dop, np.arange(0, u_int, 0.1),0., name)
        plt.plot(data['u_int'], data['zeta'][0, :, 0], label='N={}'.format(str(band)))

    ssplt.label_saves('Z_half_multiorb.png')

plot_degbandshalffill()

