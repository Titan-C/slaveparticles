# -*- coding: utf-8 -*-
"""
@author: Óscar Nájera
Created on Wed Jun  4 14:37:13 2014
"""
from __future__ import division, absolute_import, print_function
import numpy as np
from scipy import linalg as LA
from slaveparticles.quantum import operators as qo
from slaveparticles.spins.spinon import spinflipandhop

def test_spinmatrix():
    """Verifies spin matrices commute"""
    orbitals = 3
    Sx = [qo.spin_x(orbitals*2, spin) for spin in range(orbitals*2)]
    Sy = [qo.spin_y(orbitals*2, spin) for spin in range(orbitals*2)]
    Sz = [qo.spin_z(orbitals*2, spin) for spin in range(orbitals*2)]

    for S in Sx + Sy + Sz:
        assert (np.dot(S, S) - 1/4.*np.eye(4**orbitals) < 5e-16).all()
        assert (abs(qo.anticommutator(S, S) - 1/2.*np.eye(4**orbitals)) < 5e-16).all()
        assert LA.norm(S - S.T.conjugate()) < 5e-16

    #same slave particle
    for sx, sy, sz in zip(Sx, Sy, Sz):
        assert (qo.commutator(sx, sz) == -1j*sy).all()
        assert (qo.commutator(sx, sy) ==  1j*sz).all()
        assert (qo.commutator(sy, sz) ==  1j*sx).all()
    for sx, sy, sz in zip(Sx, Sy, Sz):
        assert (abs(qo.anticommutator(sx, sz)) < 5e-16).all()
        assert (abs(qo.anticommutator(sx, sy)) < 5e-16).all()
        assert (abs(qo.anticommutator(sy, sz)) < 5e-16).all()

    #to other slaves
    for S in [Sx, Sy, Sz]:
        for i in range(len(Sx)):
            for j in range(i+1, len(Sx)):
                assert (abs(qo.commutator(S[i], Sx[j])) < 5e-16).all()
                assert (abs(qo.commutator(S[i], Sy[j])) < 5e-16).all()
                assert (abs(qo.commutator(S[i], Sz[j])) < 5e-16).all()

    for op in sum(Sx) + sum(Sz):
        assert isinstance(op, np.ndarray)

def test_spinflipandhop():
    """Test the operator has the required amount of non-zero entries"""
    for orbitals in range(2, 6):
        sfh = spinflipandhop(orbitals*2)
        assert sfh.nnz == orbitals*(orbitals-1)/2 * 4**(orbitals-1)


if __name__ == "__main__":
    pass
