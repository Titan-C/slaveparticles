# -*- coding: utf-8 -*-
"""
@author: Óscar Nájera
Created on Wed Jun 18 13:05:04 2014
"""
from __future__ import division, absolute_import, print_function
import numpy as np
from scipy.optimize import root

from slaveparticles.quantum.dos import bethe_findfill_zeroT, \
     bethe_find_crystalfield, bethe_ekin_zeroT, bethe_filling_zeroT
from slaveparticles.quantum.operators import diagonalize, expected_value, \
    spin_gen, spin_z

# Fermionic part of the slave spin formulation
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

# Spin form of the Slave spin formulation
def estimate_gauge(density):
    """Calculates the gauge term for the generic spin matrices in the
       case of single site"""
    return 1/np.sqrt(density*(1-density)) - 1

from scipy.sparse import csr_matrix
def spinflipandhop(slaves):
    """Calculates the interaction term of a spin flip and pair hopping"""

    Sdw = [csr_matrix(spin_gen(slaves, i, 0)) for i in range(slaves)]
    Sup = [mat.T for mat in Sdw]

    sfh = 0.
    orbitals = slaves//2
    for n in range(orbitals):
        for m in range(n+1, orbitals):
            sfh = Sup[2*n  ] * Sdw[2*n+ 1] * Sup[2*m + 1] * Sdw[2*m  ] + sfh
            sfh = Sup[2*n+1] * Sdw[2*n   ] * Sup[2*m    ] * Sdw[2*m+1] + sfh

            sfh = Sup[2*n] * Sup[2*n+ 1] * Sdw[2*m] * Sdw[2*m+1] + sfh
            sfh = Sup[2*m] * Sup[2*m+ 1] * Sdw[2*n] * Sdw[2*n+1] + sfh

    return sfh

def spin_z_op(param, oper):
    """Generates the required Sz operators, given the system parameter setup
       and the operator dictionary"""
    slaves = param['slaves']
    oper['Sz'] = np.array([spin_z(slaves, spin) for spin in range(slaves)])
    oper['Sz+1/2'] = oper['Sz'] + 0.5*np.eye(2**slaves)
    oper['sumSz2'] = oper['Sz'].sum(axis=0)**2 #because Sz is diagonal
    Sz_mat_shape = oper['Sz'].reshape(param['orbitals'], 2, 2**slaves, 2**slaves)
    oper['sumSz-sp2'] = (Sz_mat_shape.sum(axis=1)**2).sum(axis=0)
    oper['sumSz-or2'] = (Sz_mat_shape.sum(axis=0)**2).sum(axis=0)

def spin_gen_op(oper, gauge):
    """Generates the generic spin matrices for the system"""
    slaves = len(gauge)
    oper['O'] = np.array([spin_gen(slaves, i, c) for i, c in enumerate(gauge)])
    oper['O_d']  = np.transpose(oper['O'], (0, 2, 1))
    oper['O_dO'] = np.einsum('...ij,...jk->...ik', oper['O_d'], oper['O'])
    oper['Sfliphop'] = spinflipandhop(slaves)

class Spinon(object):
    """Holds the matrix operators for a single site slave-spin system"""

    def __init__(self, **kwargs):
        """
        Generates general matrices

        Parameters
        ----------
        slaves : scalar
            number of slave particles
        orbitals : scalar
            number of orbital to alocate particles
        avg_particles : scalar
            real particle ocupacy
        hopping: list
            the hooping amplitude between each orbital of the same type

        """
        self.param = {'slaves': 2,
                      'orbitals': 1,
#                      'avg_particles': 1, MAKE SUM OF POPULATIONS
                      'populations': 0.5*np.ones(2),
                      'hopping': 0.5*np.ones(2),
                      'orbital_e': np.zeros(2),
                      'tol': 1e-8}
        for key, item in kwargs.items():
            self.param[key] = item

        self.param['lambda'] = np.zeros(self.param['slaves'])

        self.H_s = None
        self.eig_energies = None
        self.eig_states = None
        self.oper = {}

        self.set_filling(self.param['populations'])
        spin_z_op(self.param, self.oper)

        self.selfconsistency(0, 0)
        self.param['orbital_e_free'] = bethe_find_crystalfield(\
                            self.param['populations'], self.param['hopping'])
        self.param['orbital_e_shift']=orbital_energies(self.param,self.quasiparticle_weight())

    def set_filling(self, populations):
        """Sets the orbital enenergies for on the reference of the free case.
          By setting the desired local populations on every orbital.
          Then generate the necesary operators to respect such configuraion"""
        populations=np.asarray(populations)
#
#        self.param['orbital_e'] -= bethe_findfill_zeroT( \
#                                        self.param['avg_particles'],
#                                        self.param['orbital_e'],
#                                        self.param['hopping'])
        efermi = - bethe_find_crystalfield(\
                            populations, self.param['hopping'])
        self.param['populations'] = populations
#        fermion_avg(efermi, self.param['hopping'], 'ocupation')
        self.param['ekin'] = fermion_avg(efermi, self.param['hopping'], 'ekin')

        spin_gen_op(self.oper, estimate_gauge(populations))

    def reset(self, populations, lag, mu, u_int, j_coup, mean_f):
        """Resets the system into the last known state as given by the input
           values"""

        self.set_filling(populations)
        self.param['lambda'] = lag
        self.param['orbital_e'] = mu
        self.selfconsistency(u_int, j_coup, mean_f)

    def update_H(self, mean_field, l):
        """Updates the spin hamiltonian and recalculates its eigenbasis"""

        self.H_s = self.spin_hamiltonian(mean_field, l)
        try:
            self.eig_energies, self.eig_states = diagonalize(self.H_s)
        except np.linalg.linalg.LinAlgError:
            np.savez('errorhamil', H=self.H_s, fiel=mean_field, lamb=l)
            raise
        except ValueError:
            np.savez('errorhamil', H=self.H_s, fiel=mean_field, lamb=l)
            print(mean_field, l)
            raise

    def spin_hamiltonian(self, h, l):
        """Constructs the single site spin Hamiltonian"""
        H_s  = np.einsum('i,ijk', h[1], self.oper['O'])
        H_s += np.einsum('i,ijk', h[0], self.oper['O_d'])
        H_s += np.einsum('i,ijk', l, self.oper['Sz+1/2'])

        H_s += self.oper['Hint']
        return H_s

    def inter_spin_hamiltonian(self, U_inter, J_coup):
        """Calculates the interaction Hamiltonian. The Hund coupling is a
           fraction of the coulom interaction"""
        J_coup *= U_inter
        Hint  = (U_inter - 2*J_coup)/2.*self.oper['sumSz2']
        Hint += J_coup*self.oper['sumSz-sp2']
        Hint -= J_coup/2.*self.oper['sumSz-or2']
        Hint -= J_coup*self.oper['Sfliphop']

        return Hint


    def expected(self, observable, beta=1e5):
        """Wrapper to the expected_value function to fix the eigenbasis"""
        return expected_value(observable,
                         self.eig_energies,
                         self.eig_states,
                         beta)#this helps to account for zero temperature

    def quasiparticle_weight(self):
        """Calculates quasiparticle weight"""
        return np.array([self.expected(op)**2 for op in self.oper['O']])


    def mean_field(self):
        """Calculates mean field"""
        mean_field = []
        for sp_oper in [self.oper['O'], self.oper['O_d']]:
            avgO = np.array([self.expected(op) for op in sp_oper])
            avgO[abs(avgO) < 1e-10] = 0.
            mean_field.append(avgO*self.param['ekin'])

        return np.array(mean_field)

    def selfconsistency(self, U_inter, J_coup, mean_field_prev=None):
        """Iterates over the hamiltonian to get the stable selfcosistent one"""
        if mean_field_prev == None:
            mean_field_prev = np.array([self.param['ekin']]*2)
        hlog = [mean_field_prev]

        self.oper['Hint'] = self.inter_spin_hamiltonian(U_inter, J_coup)
        converging = True
        half_fill = (self.param['populations'] == 0.5).all()
        while converging:
            if half_fill:
                self.update_H(hlog[-1], self.param['lambda'])
            else:
                res = root(self.restriction, self.param['lambda'], (hlog[-1]))#, method='lm')
                if not res.success:
                    res.x = res.x *0.5 + 0.5*self.param['lambda']
                    self.update_H(self.mean_field()*0.5 +0.5*hlog[-1], res.x)
                    print('fail', self.param['populations'][3:5])
                if (self.quasiparticle_weight() <0.001 ).all():
                    return hlog
                self.param['lambda'] = res.x

            hlog.append(self.mean_field())
            converging = (abs(hlog[-1] - hlog[-2]) > self.param['tol']).all() \
                    or (abs(self.restriction(self.param['lambda'], hlog[-1])) > self.param['tol']).all()

        return hlog

    def restriction(self, lam, mean_field):
        """Lagrange multiplier in lattice slave spin"""
        self.update_H(mean_field, lam)
        restric = np.array([self.expected(op) - n for op, n in zip(self.oper['Sz+1/2'], self.param['populations'])])
        return restric
