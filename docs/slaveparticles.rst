Slave Particles
===============

Representation of physical electron operators. The fermion is expressed in
terms of constrained(slave) auxiliary fields that enlarge the Hilbert space
but are subject to local restrictions that eliminate unphysical states.


Slave Spin
----------

Start from the observation of the possible occupancies of a spinless fermion on
a given site, $n_d=0$ and $n_d=1$, they can be represented as two possible
states of a spin-$1/2$ variable, $S^z=-1/2$ and $S^z=+1/2$. To ensure, in this
fermionic context, the anti-commutation properties and auxiliary fermion $f$ is
introduced, with the additional local constrain:

.. math::  f^\dagger f = S^z + \frac{1}{2}
   :label: slavespinrestriction

In this manner vacuum state and occupied states are represented as

.. math::
   \ket{\emptyset} & = \ket{n_f=0,S^z = -1/2} \\
   \ket{1} & \equiv d^\dagger\ket{\emptyset} = \ket{n_f=1,S^z = +1/2}

And the constrain eliminates the unphysical states

.. math::
   \ket{n_f=1,S^z = -1/2} \\
   \ket{n_f=0,S^z = +1/2}

It is important to note that the spin-$1/2$ variable has nothing to do with the
actual spin of the particle. Spins are treated independently and the spin-$1/2$
variable is just to give the presence of the particle. To exemplify this, the
orbital base $\{\ket{\emptyset}, \ket{\uparrow}, \ket{\downarrow}, \ket{\uparrow
\downarrow}\}$ is written in this extended Hilbert space as:

.. math::
   \ket{\emptyset} &= \ket{n_f = 0; S^z_\uparrow = -1/2}\ket{n_f = 0; S^z_\downarrow = -1/2} \\
   \ket{\uparrow} &= \ket{n_f = 1; S^z_\uparrow = 1/2}\ket{n_f = 0; S^z_\downarrow = -1/2} \\
   \ket{\downarrow} &= \ket{n_f = 0; S^z_\uparrow = -1/2}\ket{n_f = 1; S^z_\downarrow = 1/2} \\
   \ket{\uparrow\downarrow} &= \ket{n_f = 1; S^z_\uparrow = 1/2}\ket{n_f = 1; S^z_\downarrow = 1/2}

The isolated atom
'''''''''''''''''

The Hamiltonian for the isolated atom, with degenerated $2N$ fermions in spin
and orbital reads:

.. math::
   \mathcal{H} = \frac{U}{2} \left( \sum_{n} d_n^\dagger d_n - N \right)^2
    -\mu \sum_{n} d_n^\dagger d_n

When introducing the new operators in this Hamiltonian the Lagrange multipliers
are also included to account for the restriction :eq:`slavespinrestriction` to
avoid the unphysical states.

.. math::
   \mathcal{H} = \frac{U}{2} \left( \sum_{n} S_n^z \right)^2
     + \lambda \sum_{n} \left( S_n^z +\frac{1}{2} - f_n^\dagger f_n \right)
    -\mu \sum_{n} f_n^\dagger f_n

This separates in 2 Hamiltonians:

.. math::
   \mathcal{H}_f &= -(\lambda + \mu) \sum_{n} f_n^\dagger f_n \\
   \mathcal{H}_s &= \frac{U}{2} \left( \sum_{n} S_n^z \right)^2
                    +\lambda \sum_{n} (S_n^z + \frac{1}{2})

The spin Hamiltonian has a spectrum :math:`E_{Q \leftarrow S+1/2} = U/2(Q-N)^2 +\lambda Q`,
where $Q$ is the present charge. To find the Lagrange multiplier $\lambda$ take
$0=\frac{\partial <\mathcal{H}>}{\partial \lambda}$ because restrictions are treated always in average.

.. math::
   0 &=-\sum_{n} < f_n^\dagger f_n>_f + \sum_{n} <S_n^z>_s \\
   2Nn_F(\lambda +\mu) &= 2N<S_n^z + \frac{1}{2}> = <S^z> + N = <Q> \\
   2Nn_F(\lambda +\mu) &=
   \mathcal{Z}^{-1} \sum_{Q=0}^{2N} \binom{2N}{Q} Q e^{\beta(U/2(Q-N)^2 +\lambda Q)}


Where every term is averaged with its corresponding Hamiltonian, $n_F$ is the Fermi
distribution and the Grand Canonical Partition function is
$\mathcal{Z}^{-1} \sum_{Q=0}^{2N} \binom{2N}{Q} Q e^{\beta(U/2(Q-N)^2 +\lambda Q)}$.
Then by numerical root finding the
multiplier $\bar{\lambda}(\mu,\beta)$ allows for to describe the mean fermion occupation,
which is $2Nn_F(\mu + \bar{\lambda}(\mu,\beta))$ and can recuperate the complete Coulomb ladder. It
then has to be compared to the exact solution:

.. math::
   2N<d^\dagger d> =  \mathcal{Z}^{-1} \sum_{Q=0}^{2N} \binom{2N}{Q} Q e^{\beta(U/2(Q-N)^2 -\mu Q)}

.. plot:: Luca.py

The lattice model
'''''''''''''''''

The Hamiltonian in this case when there is no Hund coupling($J=0$) and in a lattice reads:

.. math::
   \mathcal{H} = -\sum_m t_m \sum_{<i,j>, \sigma} (d^\dagger_{im\sigma}d_{jm\sigma} +h.c.)
    + \sum_{im\sigma}(\epsilon_m - \mu)d^\dagger_{im\sigma}d_{im\sigma}
    + \frac{U}{2} \sum_i \left( \sum_{m\sigma} d_{im\sigma}^\dagger d_{im\sigma} - N \right)^2

The focus now for simplicity is the case of zero crystal-field splitting
$\epsilon_m=0$ and half-filling of each band one electron per site in each
orbital $\mu=0$.

Then when dealing with a multi orbital system, $2N$  new spin-$1/2$ variables,
$S^z_{m\sigma}$ and $2N$ auxiliary fermions $f_{m\sigma}$ are introduced, where
$m=1, \cdots, N$ is the number of orbitals. And the local constrain is applied
on each lattice site($i$):

.. math::  f_{im\sigma}^\dagger f_{im\sigma} = S_{im\sigma}^z + \frac{1}{2}
   :label: slavespinrestriction_multiorbitalsite

using the Lagrange multiplier $\lambda_{im\sigma}$.

When rewriting the Hamiltonian in terms of the auxiliary fermions and the slave
spins the interaction term turn easily into:

.. math:: \mathcal{H}_{int} = \sum_i \left( \sum_{m\sigma} S^z_{im\sigma} \right)^2

For the non interacting part, an appropriate representation of the creation
operator has to be chosen. The direct possibility $d^\dagger \rightarrow S^+ f^\dagger$,
although correct leads to problems with the spectral weight conservation because
$S^+$ and $S^-$ don't commute. Instead the representation $d^\dagger \rightarrow
2S^xf^\dagger$, $d \rightarrow 2S^xf$, which is identical on the physical Hilbert
space and involves commuting slave spin operators in chosen. Then the non interacting
Hamiltonian reads:

.. math::
   \mathcal{H}_0 = -\sum_m t_m \sum_{<i,j>, \sigma} 4S^x_{im\sigma}S^x_{jm\sigma}(f^\dagger_{im\sigma}f_{jm\sigma} +h.c.)
    + \sum_{im\sigma}(\epsilon_m - \mu)f^\dagger_{im\sigma}f_{im\sigma}
