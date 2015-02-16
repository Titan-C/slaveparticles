# -*- coding: utf-8 -*-
"""
===
oeu
===

#@author: Óscar Nájera
#Created on Wed Sep 24 14:55:15 2014
eu
"""

from __future__ import division, absolute_import, print_function

import slavespins.storage as st
def follow_cf(save, Uspan, target_cf, nup, n_tot=5.0, slsp=None):
    """Calculates the quasiparticle weight in single
       site spin hamiltonian under with N degenerate half-filled orbitals """

    if slsp == None:
        slsp = Spinon(slaves=6, orbitals=3, avg_particles=n_tot,
                       hopping=[0.5]*6, populations = np.asarray([n_tot]*6)/6)

    zet, lam, mu, mean_f = [], [], [], []

    for co in Uspan:
        print('U=', co, 'del=', target_cf)
        res=root(targetpop, nup[-1],(co,target_cf,slsp, n_tot))
        print(res.x)
        if res.x>nup[-1]: break

        nup.append(res.x)
        slsp.param['populations']=population_distri(nup[-1])
        mean_f.append(slsp.mean_field())
        zet.append(slsp.quasiparticle_weight())
        lam.append(slsp.param['lambda'])
        mu.append(orbital_energies(slsp.param, zet[-1]))

#    plt.plot(np.asarray(zet)[:,0], label='d={}, zl'.format(str(target_cf)))
#    plt.plot(np.asarray(zet)[:,5], label='d={}, zh'.format(str(target_cf)))

    case = save.createGroup('cf={}'.format(target_cf))
    varis = st.setgroup(case)
    st.storegroup(varis, Uspan[:len(zet)], zet, lam, mu, nup[1:],target_cf,mean_f)
#    save.close()

def targetpop(upper_density, coul, target_cf, slsp, n_tot):
    """restriction on finding the right populations that leave the crystal
    field same"""
    if upper_density < 0.503: return 0.
    trypops=population_distri(upper_density, n_tot)
    slsp.set_filling(trypops)
    slsp.selfconsistency(coul,0)
    efm_free = dos_bethe_find_crystalfield(trypops, slsp.param['hopping'])
    orb_ener = slsp.param['lambda']+ slsp.quasiparticle_weight()*efm_free
    obtained_cf = orb_ener[5] - orb_ener[0]
    return target_cf - obtained_cf

def population_distri(nup, n_tot=5.0):
    nup=float(nup)
    ndw=(n_tot-2*nup)/4
    return np.asarray([ndw]*4+[nup]*2)

#### Trying on lucas statement
#But it doesn't work, fixing the fermion orbital energy from the free case
#solution is bad when raising the interaction. The lagrange multiplier has
#to high values compared to this energy and the band populations change to dras
#tically even at low interaction
#def targetpop(density,coul,hund, slsp):
#    """restriction on finding the right populations that leave the crystal
#    field same"""
#    ndw,nup = density
#
#    trypops=np.asarray([ndw]*4+[nup]*2)
#    slsp.set_filling(trypops)
#    slsp.selfconsistency(coul,hund)
##    orb_ener = orbital_energies(slsp.param, slsp.quasiparticle_weight())
##    obtained_cf = orb_ener[5] - orb_ener[0]
#    efermi = slsp.param['lambda'] - slsp.param['orbital_e_shift']
#    res_pops = fermion_avg(efermi, slsp.param['hopping'], 'ocupation')
##    return [(trypops-res_pops)[-1],ndw-(5-2*nup)/4.]
#    return [ndw-res_pops[0], nup-res_pops[5]]

if __name__ == "__main__":
    save = st.Dataset('dop3b5_02e.nc', 'w', format='NETCDF4')
    deltacf = np.arange(0., 2.51, 0.1)
    coul = -1.84*deltacf+3.5
    for cf, ul in zip(deltacf, coul):
#        Uspan = np.concatenate((np.arange(0, ul, 0.1), \
#                                np.arange(ul, 4.651, 0.02)))
        Uspan = np.arange(0, 6, 0.1)
        nglo = 5.02
        nhi = [nglo/6.]
        follow_cf(save, Uspan, cf, nhi, nglo)

    save.close()
