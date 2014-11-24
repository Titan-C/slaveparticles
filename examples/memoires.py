# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 22:37:48 2014

@author: oscar
plots of report
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
        plt.plot(data['u_int'], data['zeta'][0, :, 0], label='$N={}$'.format(str(band)))

    ssplt.label_saves('Z_half_multiorb.png')

#band dop
def plot_dop(bands, int_max, dop, hund_cu, name):
    """Plot of Quasiparticle weight for N degenerate bands
       under selected doping shows transition only at half-fill
       the rest are metallic states"""
    data = ssplt.calc_z(bands, dop, np.arange(0, int_max, 0.1), hund_cu, name)
    ssplt.plot_curves_z(data, name)

#band dop_phasediag
def plot_dop_phase(bands, int_max, hund_cu):
    """Phase plot of Quasiparticle weight for N degenerate bands
       under doping shows transition only at interger filling
       the rest are metallic states"""
    name = 'Z_dop_phase_'+str(bands)+'bands_U'+str(int_max)+'J'+str(hund_cu)
    dop = np.sort(np.hstack((np.linspace(0.01,0.99,50),
                    np.arange(1./2./bands, 1, 1/2/bands))))
    data = ssplt.calc_z(bands, dop, np.arange(0, int_max, 0.1), hund_cu, name)

    ssplt.imshow_z(data, name)
    ssplt.surf_z(data, name)

if __name__ == "__main__":
    plot_degbandshalffill()
    plot_dop(1, 4.6, [0.5, 0.499, 0.495, 0.49, 0.45, 0.4, 0.2, 0.1], 0., 'Z_dop_1orb')
    plot_dop(3, 7, np.asarray([3, 4, 5])/6, 0., 'Z_3b_int_fil')
    for U, N in zip([4.6, 6, 8], range(1, 4)):
        plot_dop_phase(N, U, 0.)
