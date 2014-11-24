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

# Hund Coupling as L. de'Medici PRB 83, 205112 (2011)
def hund_coup(bands, dop, u_lim):
    hund_cu = [0, 0.05, 0.1, 0.15, 0.2, 0.25]
    plt.figure()
    for int_max, hund in zip(u_lim, hund_cu):
        name = 'Z_'+str(bands)+'bands_U'+str(int_max)+'J'+str(hund)+'n'+str(dop)
        data = ssplt.calc_z(2, [dop], np.arange(0, int_max, 0.1), hund, name)
        plt.plot(data['u_int'], data['zeta'][0, :, 0], label='$J/U={}$'.format(str(hund)))
    plt.title('Quasiparticle weight for {} electron(s) in {} bands'.format(str(dop*2*bands), str(bands)))
    ssplt.label_saves('Z_{}bands_{}n_Hund.png'.format(str(bands), str(dop)))

if __name__ == "__main__":
    hund_coup(2, 0.5, [5.2, 3.8, 3.4, 3.1, 3.0, 3.0])
    hund_coup(2, 0.25, [5.0, 5.5, 6.0, 7.5, 10.0, 15.0])
