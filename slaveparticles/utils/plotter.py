# -*- coding: utf-8 -*-
"""
@author: Óscar Nájera
Created on Tue Apr  8 14:09:37 2014
"""

import numpy as np
from slaveparticles.spins import Spinon, orbital_energies
import matplotlib.pyplot as plt

def solve_loop(slsp, u_span, j_coup):
    """Calculates the quasiparticle for the input loop of:
        @param slsp: Slave spin Object
        @param Uspan: local Couloumb interation
        @param J_coup: Fraction of Uspan of Hund coupling strength"""

    zet, lam, eps, hlog, mean_f = [], [], [], [], [None]

    for u in u_span:
        print(u, j_coup)
        hlog.append(slsp.selfconsistency(u, j_coup, mean_f[-1]))
        mean_f.append(slsp.mean_field())
        zet.append(slsp.quasiparticle_weight())
        lam.append(slsp.param['lambda'])
        eps.append(orbital_energies(slsp.param, zet[-1]))

    return np.asarray([zet, lam, eps]), hlog, mean_f

def calc_z(bands, filling, interaction, hund_cu, name):
    """Calculates the quasiparticle weight of degenerate system of N-bands
       at a given filling within an interaction range and saves the file"""

    while True:
        try:
            data = np.load(name+'.npz')
            break
        except IOError:
            dopout = []
            for dop in filling:
                slsp = Spinon(slaves=2*bands, orbitals=bands, \
                             hopping=[0.5]*2*bands, populations=[dop]*2*bands)
                dopout.append(solve_loop(slsp, interaction, hund_cu)[0][0])
            np.savez(name, zeta=dopout, u_int=interaction, doping=filling, hund=hund_cu)

    return data

def label_saves(name):
    """Labels plots and saves file"""
    plt.legend(loc=0)
    plt.ylim([0, 1.025])
    plt.xlabel('$U/D$', fontsize=20)
    plt.ylabel('$Z$', fontsize=20)
    plt.savefig(name, dpi=300, format='png',
            transparent=False, bbox_inches='tight', pad_inches=0.05)

def plot_curves_z(data, name, title=None):
    """Generates a simple plot of the quasiparticle weight decay curves given
       data object with doping setup"""

    plt.figure()
    for zet, c in zip(data['zeta'], data['doping']):
        plt.plot(data['u_int'], zet[:, 0], label='$n={}$'.format(str(c)))
    if title != None:
        plt.title(title)
    label_saves(name+'.png')

def pick_flat_z(data):
    """Generate a 2D array of the quasiparticle weight by only selecting the
    first particle data"""
    zmes = []
    for i in data['zeta']:
        zmes.append(i[:, 0])
    return np.asarray(zmes)

def imshow_z(data, name):
    """2D color plot of the quasiparticle weight as a function of interaction
       and doping"""

    zmes = pick_flat_z(data)

    plt.figure()
    plt.imshow(zmes.T, origin='lower', \
            extent=[data['doping'].min(), data['doping'].max(), \
            0, data['u_int'].max()], aspect=.16)
    plt.colorbar()
    plt.xlabel('$n$', fontsize=20)
    plt.ylabel('$U/D$', fontsize=20)
    plt.savefig(name+'_imshow.png', dpi=300, format='png',
            transparent=False, bbox_inches='tight', pad_inches=0.05)

def surf_z(data, name):
    """Surface plot of the quasiparticle weight as fuction of U/D and dop"""
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    dop, u_int = np.meshgrid(data['doping'], data['u_int'])
    zmes = pick_flat_z(data)
    ax.plot_surface(u_int, dop, zmes.T, rstride=1, cstride=1, cmap=cm.coolwarm,
            linewidth=0, antialiased=False)
    ax.contourf(u_int, dop, zmes.T, zdir='z', offset=0, cmap=cm.coolwarm)

    ax.set_zlim(0, 1)
    ax.view_init(elev=11, azim=-34)


    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    ax.set_xlabel('$U/D$', fontsize=20)
    ax.set_ylabel('$n$', fontsize=20)
    ax.set_zlabel('$Z$', fontsize=20)

    plt.savefig(name+'_surf.png', dpi=300, format='png',
            transparent=False, bbox_inches='tight', pad_inches=0.05)

def plot_mean_field_conv(N=1, n=0.5, Uspan=np.arange(0, 3.6, 0.5)):
    """Generates the plot on the convergenge of the mean field in single
       site spin hamiltonian under with N degenerate half-filled orbitals """

    sl = Spinon(slaves=2*N, orbitals=N, avg_particles=2*n,
                       hopping=[0.5]*2*N, orbital_e=[0]*2*N)
    hlog = solve_loop(sl, Uspan, [0.])[1]

    f, (ax1, ax2) = plt.subplots(2, sharex=True)
    for field in hlog:
        field = np.asarray(field)
        ax1.semilogy(abs(field[1:]-field[:-1]))
        ax2.plot(field)#, label = 'h, U = {}'.format(Uint))

    plt.title('Convergence of selfconsintent mean field')
    ax1.set_ylabel('$\\Delta h$')
    ax2.set_ylabel('mean field $h$')
    plt.xlabel('iterations')

    return hlog
