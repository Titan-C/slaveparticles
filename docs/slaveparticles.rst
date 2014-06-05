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
actual spin of the particle, it is only a label. The slave spin-$1/2$
variable is introduced for every fermion species to give the presence of the particle,
as such slave spins are also labeled with a spin index, $S^z_{\sigma}$. To exemplify this, the
orbital base $\{\ket{\emptyset}, \ket{\uparrow}, \ket{\downarrow}, \ket{\uparrow
\downarrow}\}$ is written in this extended Hilbert space as:

.. math::
   \ket{\emptyset} &= \ket{n_f = 0; S^z_\uparrow = -1/2}\ket{n_f = 0; S^z_\downarrow = -1/2} \\
   \ket{\uparrow} &= \ket{n_f = 1; S^z_\uparrow = 1/2}\ket{n_f = 0; S^z_\downarrow = -1/2} \\
   \ket{\downarrow} &= \ket{n_f = 0; S^z_\uparrow = -1/2}\ket{n_f = 1; S^z_\downarrow = 1/2} \\
   \ket{\uparrow\downarrow} &= \ket{n_f = 1; S^z_\uparrow = 1/2}\ket{n_f = 1; S^z_\downarrow = 1/2}

Then when dealing with a multi orbital system, $2N$  new spin-$1/2$ variables
$S^z_{m\sigma}$ and $2N$ auxiliary fermions $f_{m\sigma}$ are introduced, where
$m=1, \cdots, N$ is the number of orbitals.

Case: The isolated atom
'''''''''''''''''''''''

The Hamiltonian for the isolated atom, with $2N$ degenerated fermions in spin
and orbital in particle-hole symmetry formulation, reads:

.. math::
   \mathcal{H} = \frac{U}{2} \left( \sum_{n} d_n^\dagger d_n - N \right)^2
    -\mu \sum_{n} d_n^\dagger d_n

When introducing the new operators in this Hamiltonian a Lagrange multiplier
$\lambda$, considering all particles are indistingushable in spin and orbital,
has also to included to account for the restriction :eq:`slavespinrestriction`
avoid the unphysical states.

.. math::
   \mathcal{H} = \frac{U}{2} \left( \sum_{n} S_n^z \right)^2
     + \lambda \sum_{n} \left( S_n^z +\frac{1}{2} - f_n^\dagger f_n \right)
    -\mu \sum_{n} f_n^\dagger f_n

Then is possible to separate it in 2 Hamiltonians:

.. math::
   \mathcal{H}_f &= -(\lambda + \mu) \sum_{n} f_n^\dagger f_n \\
   \mathcal{H}_s &= \frac{U}{2} \left( \sum_{n} S_n^z \right)^2
                    +\lambda \sum_{n} (S_n^z + \frac{1}{2})

The spin Hamiltonian has a spectrum :math:`E_{Q \leftarrow S+1/2} = U/2(Q-N)^2 +\lambda Q`,
where $Q$ is the charge present. To find the Lagrange multiplier $\lambda$ take
$0=\frac{\partial <\mathcal{H}>}{\partial \lambda}$ because restrictions are treated always in average.

.. math::
   0 &=-\sum_{n} < f_n^\dagger f_n>_f + \sum_{n} <S_n^z>_s \\
   2Nn_F(-\lambda -\mu) &= 2N<S_n^z + \frac{1}{2}> = <S^z> + N = <Q> \\
   2Nn_F(-\lambda -\mu) &=
   \mathcal{Z}^{-1} \sum_{Q=0}^{2N} \binom{2N}{Q} Q \exp({\beta(U/2(Q-N)^2 +\lambda Q)})


Where every term is averaged with its corresponding Hamiltonian, $n_F$ is the Fermi
distribution and the Grand Canonical Partition function from the spin
hamiltonian is
$\mathcal{Z}= \sum_{Q=0}^{2N} \binom{2N}{Q} \exp({\beta(U/2(Q-N)^2 +\lambda Q)})$.
Then by numerical root finding the
multiplier $\bar{\lambda}(\mu,\beta)$ can be estimated and allows to describe the mean fermion occupation,
which is $2Nn_F(-\mu - \bar{\lambda}(\mu,\beta))$ and can recuperate the complete Coulomb ladder.
Comparing to the exact solution:

.. math::
   2N<d^\dagger d> =  \mathcal{Z}^{-1} \sum_{Q=0}^{2N} \binom{2N}{Q} Q e^{\beta(U/2(Q-N)^2 -\mu Q)}

As shown in the next plot, the slave spin approximation is capable of
recovering the coulomb ocupation ladder, for the isolated atom with degenerate
fermions in spin and orbital. The approximation works best around half-filling.

.. plot:: Luca.py


Case: The lattice model
'''''''''''''''''''''''

When in a lattice, atoms have overlapping orbitals and electrons are capable to
move along this lattice. Then for the hamiltonian this term needs to be
included appearing in the Tight-Binding formulation. Then as simple extension
of the previous isolated atom case and in a multiorbital scenario, the
Hamiltonian reads. The focus now for simplicity is the case of zero crystal-field splitting
$\epsilon_m=0$ and half-filling of each band one electron per site in each
orbital $\mu=0$.

.. math::
   \mathcal{H} = -\sum_m t_m \sum_{<i,j>, \sigma} (d^\dagger_{im\sigma}d_{jm\sigma} +h.c.)
    + \sum_{im\sigma}(\epsilon_m - \mu)d^\dagger_{im\sigma}d_{im\sigma}
    + \frac{U}{2} \sum_i \left( \sum_{m\sigma} d_{im\sigma}^\dagger d_{im\sigma} - N \right)^2

Here it is needed to enforce the restriction:

.. math::  f_{im\sigma}^\dagger f_{im\sigma} = S_{im\sigma}^z + \frac{1}{2}
   :label: slavespinrestriction_multiorbitalsite

using the Lagrange multiplier $\lambda_{im\sigma}$, which can be used declaring
specific contrains to lattice site, orbital, and spin.

When rewriting the Hamiltonian in terms of the auxiliary fermions and the slave
spins the interaction term turns easily into:

.. math:: \mathcal{H}_{int} = \frac{U}{2} \sum_i \left( \sum_{m\sigma} S^z_{im\sigma} \right)^2

For the non interacting part, an appropriate representation of the creation
operator has to be chosen. The direct possibility $d^\dagger \rightarrow S^+ f^\dagger$,
although correct leads to problems with the spectral weight conservation because
$S^+$ and $S^-$ don't commute. Instead the representation $d^\dagger \rightarrow
2S^xf^\dagger$ and $d \rightarrow 2S^xf$ is chosen, which is identical on the physical Hilbert
space and involves commuting slave spin operators. Then the non interacting
Hamiltonian reads:

.. math::
   \mathcal{H}_0 = -\sum_m t_m \sum_{<i,j>, \sigma} 4S^x_{im\sigma}S^x_{jm\sigma}(f^\dagger_{im\sigma}f_{jm\sigma} +h.c.)
   + \sum_{im\sigma}(\epsilon_m - \mu)f^\dagger_{im\sigma}f_{im\sigma}


The constrain is treated on average using a static and
site, orbital and particle independent Lagrange multiplier $\lambda$.
Then the Hamiltonian reads:


.. math:: \mathcal{H} = &\frac{U}{2} \sum_i \left( \sum_{m\sigma} S^z_{im\sigma} \right)^2 \\
   &-\sum_m t_m \sum_{<i,j>, \sigma} 4S^x_{im\sigma}S^x_{jm\sigma}(f^\dagger_{im\sigma}f_{jm\sigma} +h.c.) \\
   &+\lambda\sum_{im\sigma} \left( S_{im\sigma}^z + \frac{1}{2} - f_{im\sigma}^\dagger f_{im\sigma} \right)

Using a Hartree-Fock approximation for the operators $S$ and $f$:

.. math::
   S^x_{im\sigma}S^x_{jm\sigma}(f^\dagger_{im\sigma}f_{jm\sigma} +h.c.)
   \approx <S^x_{im\sigma}S^x_{jm\sigma}>(f^\dagger_{im\sigma}f_{jm\sigma}
   +h.c.) \\
   +S^x_{im\sigma}S^x_{jm\sigma}<f^\dagger_{im\sigma}f_{jm\sigma} +h.c.> \\
   -<S^x_{im\sigma}S^x_{jm\sigma}(f^\dagger_{im\sigma}f_{jm\sigma} +h.c.)>

it is then possible to decouple the Hamiltonian into two effective ones:

.. math:: \mathcal{H}^f_{eff} = &-\sum_m t_m^{eff} \sum_{<i,j>, \sigma} (f^\dagger_{im\sigma}f_{jm\sigma} +h.c.) \\
   &-\sum_{im\sigma} \lambda_m f_{im\sigma}^\dagger f_{im\sigma}
   :label: hamileff_fermion

.. math:: \mathcal{H}^S_{eff} = &-\sum_m 4J^{eff}_m \sum_{<i,j>, \sigma} S^x_{im\sigma}S^x_{jm\sigma} \\
   &+\sum_{im\sigma} \lambda_m \left( S_{im\sigma}^z + \frac{1}{2} \right)
   +\frac{U}{2} \sum_i \left( \sum_{m\sigma} S^z_{im\sigma} \right)^2
   :label: hamileff_spin

Where the effective hopping and the effective exchange constants are
determined self consistently from:

.. math::
   t^{eff}_m &= 4t_m<S^x_{im\sigma}S^x_{jm\sigma}>
   :label:
           eff_hopping
.. math::
   J^{eff}_m &= t_m<f^\dagger_{im\sigma}f_{jm\sigma} +h.c.>
   :label: eff_exchange

The fermion field hamiltonian is a non-interacting one, and it's analytical
solution is well known. For the slave spin hamiltonian, it can be treated
in a single-site using the Weiss mean field approximation.

.. math:: \mathcal{H}_s = &\sum_{m\sigma} 2h_mS^x_{m\sigma}
   +\sum_{m\sigma} \lambda_m \left( S_{m\sigma}^z + \frac{1}{2} \right)
   +\frac{U}{2} \left( \sum_{m\sigma} S^z_{m\sigma} \right)^2
   :label: hamil_spin_meanfield


Here the mean field $h_m$ has to be determined self-consistently from:

.. math:: h
   _m = -2zJ^{eff}_m<S^x_{m\sigma}> = 4<S^x_{m\sigma}>\frac{1}{N_s}\sum_k \epsilon_{km}<f^\dagger_{km\sigma}f_{km\sigma}>

where $z$ is the coordination number, $\epsilon_{km}=-t_m\sum_{\{\vec{a}\}}e^{-i\vec{k}\cdot\vec{a}}$
with $\{\vec{a}\}$ the set of vectors to the nearest neighbors

The effective fermion hamiltonian is

.. math:: \mathcal{H}^f_{eff} = &\sum_{km\sigma} (-t_m^{eff} \sum_{\{\vec{a}\}} e^{-i\vec{k}\cdot\vec{a}} - \lambda_m) f^\dagger_{km\sigma}f_{km\sigma} \\
   &=\sum_{km\sigma} (Z_m\epsilon_{mk} - \lambda_m) f^\dagger_{km\sigma}f_{km\sigma}

where $Z_m=4<S^x_{im\sigma}>^2$ is the quasiparticle weight.

In the ordered spin basis $\{\ket{\uparrow\uparrow}, \ket{\uparrow\downarrow}, \ket{\downarrow\uparrow}, \ket{\uparrow\downarrow}\}$, where the spin labeling the operators are then

.. math::
   S^z_{\uparrow} = \frac{1}{2} \left[\begin{smallmatrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & -1 & 0\\0 & 0 & 0 & -1\end{smallmatrix}\right]
.. math::
   S^z_{\downarrow} = \frac{1}{2} \left[\begin{smallmatrix}1 & 0 & 0 & 0\\0 & -1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & -1\end{smallmatrix}\right]
.. math::
   S^x_{\uparrow} = \frac{1}{2} \left[\begin{smallmatrix}0 & 0 & 1 & 0\\0 & 0 & 0 & 1\\1 & 0 & 0 & 0\\0 & 1 & 0 & 0\end{smallmatrix}\right]
.. math::
   S^x_{\downarrow} = \frac{1}{2} \left[\begin{smallmatrix}0 & 1 & 0 & 0\\1 & 0 & 0 & 0\\0 & 0 & 0 & 1\\0 & 0 & 1 & 0\end{smallmatrix}\right]

.. plot::

   from slavespins.plotter import plot_quasiparticle_weight
   plot_quasiparticle_weight()
   plot_quasiparticle_weight(2)
   import matplotlib.pyplot as plt
   plt.show()

