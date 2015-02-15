# -*- coding: utf-8 -*-
"""
======================================================
Drop of quasiparticle weight by increasing interaction
======================================================

The quasiparticle weight of the electronic system drops as the local interaction
is increased. When studying the single band case, only the half-filled scenario
undergoes a transition into the Mott insulator. All other doping configurations
retain the strongly correlated metal.
"""

from __future__ import division, absolute_import, print_function
import slaveparticles.utils.plotter as ssplt
import numpy as np
import matplotlib.pyplot as plt

#band dop
def plot_dop(bands, int_max, dop, hund_cu, name):
    """Plot of Quasiparticle weight for N degenerate bands
    under selected doping shows transition only at half-fill
    the rest are metallic states"""
    data = ssplt.calc_z(bands, dop, np.arange(0, int_max, 0.1), hund_cu, name)
    ssplt.plot_curves_z(data, name)

plot_dop(1, 4.6, [0.5, 0.499, 0.495, 0.49, 0.45, 0.4, 0.2, 0.1], 0., 'Z_dop_1orb')

