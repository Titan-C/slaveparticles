=============================
Slave Spins Mean Field Theory
=============================

Construction
------------

Start from the observation that the possible occupancies of a real fermion($d$)
on a given site, $n_d=0$ and $n_d=1$, can be
represented as two possible states of a spin-$1/2$ variable, $S^z=-1/2$ and
$S^z=+1/2$. The idea then is to duplicate the original particle into a spin-like
degree of freedom and a fermion charge degree of freedom by introducing a spin-field and
an associated auxiliary fermion($f$). In this manner vacuum state and occupied
states are represented as:

.. math::
   \ket{\emptyset} & = \ket{n_f=0,S^z = -1/2} \\
   \ket{1} & \equiv d^\dagger\ket{\emptyset} = \ket{n_f=1,S^z = +1/2}

In this context, the anti-commutation properties of the original fermion
operator $d$ are insured by the introduction of the auxiliary fermion
operator $f$. Then, one formulates the constrain:

.. math::  f^\dagger f = S^z + \frac{1}{2}
   :label: slavespinrestriction

which insures that the number of auxiliary fermions and the spin's states are
the same and, above all, it eliminates the non-physical states

.. math::
   \ket{n_f=1,S^z = -1/2} \\
   \ket{n_f=0,S^z = +1/2}

It is important to note that the spin-$1/2$ variable has nothing to do with the
actual spin of the particle. The slave spin-$1/2$ variable is introduced for
every fermion species to give the presence of the particle. Therefore, the
original particle $d_\sigma$ is replicated into the slave variables with the
same quantum numbers: $n_{f\sigma}$ and $S^z_\sigma$. To exemplify this, the
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

In order to test the slave-spin method, we first consider the case of a
degenerate multi-orbital atom, whose analytical solution is well known. The
Hamiltonian in particle-hole symmetry formulation reads:

.. math::
   \mathcal{H} = \frac{U}{2} \left( \sum_{n} d_n^\dagger d_n - N \right)^2
    -\mu \sum_{n} d_n^\dagger d_n

The first step is to recast the Hamiltonian in terms of the new slave-spin
operators. The first choice is to use the auxiliary
fermions for the non-interacting term($d_n^\dagger d_n\rightarrow f_n^\dagger
f_n$) and the slave spins for the interacting case($d_n^\dagger d_n \rightarrow
S^z_n + \frac{1}{2}$). We seek a paramagnetic solution, then the constraint
:eq:`slavespinrestriction`  which avoids non-physical states
is included through a unique Lagrange multiplier $\lambda$, since all particles
are indistinguishable in spin and orbital quantum numbers. After this treatment
the Hamiltonian reads:

.. math::
   \mathcal{H} = \frac{U}{2} \left( \sum_{n} S_n^z \right)^2
     + \lambda \sum_{n} \left( S_n^z +\frac{1}{2} - f_n^\dagger f_n \right)
    -\mu \sum_{n} f_n^\dagger f_n

which is possible to separate into a fermion and a spin Hamiltonians:

.. math::
   \mathcal{H}_f &= -(\lambda + \mu) \sum_{n} f_n^\dagger f_n \\
   \mathcal{H}_s &= \frac{U}{2} \left( \sum_{n} S_n^z \right)^2
                    +\lambda \sum_{n} (S_n^z + \frac{1}{2})

The spin Hamiltonian has a spectrum :math:`E_{Q \leftarrow S+1/2} = U/2(Q-N)^2 +\lambda Q`,
where $Q$ is the charge present.The Lagrange multiplier $\lambda$ is fixed at the mean-field level by
determining the stationary point of the mean-field averaged Hamiltonian:
$0=\frac{\partial \langle\mathcal{H}\rangle}{\partial \lambda}$. Within this
approach, the restriction :eq:`slavespinrestriction` is therefore
respected at the mean field level too.

.. math::
   0 &=-\sum_{n} < f_n^\dagger f_n>_f + \sum_{n} <S_n^z>_s \\
   2Nn_F(-\lambda -\mu) &= 2N<S_n^z + \frac{1}{2}> = <S^z> + N = <Q> \\
   2Nn_F(-\lambda -\mu) &=
   \mathcal{Z}^{-1} \sum_{Q=0}^{2N} \binom{2N}{Q} Q \exp({\beta(U/2(Q-N)^2 +\lambda Q)})


Where every term is averaged with its corresponding Hamiltonian, $n_F$ is the Fermi
distribution and the Grand Canonical Partition function from the spin
Hamiltonian is
$\mathcal{Z}= \sum_{Q=0}^{2N} \binom{2N}{Q} \exp({\beta(U/2(Q-N)^2 +\lambda Q)})$.
Then by numerical root finding the
multiplier $\bar{\lambda}(\mu,\beta)$ can be estimated and allows to describe
the mean fermion occupation,
which is $2Nn_F(-\mu - \bar{\lambda}(\mu,\beta))$ . The
solution that we find can be conveniently represented by the total fermion
occupation as a function of the chemical potential $\mu$. The resulting curve
displays the well known Coulomb ladder-shape,
where the system acquires only integer fillings. The change from an integer
occupation to an adjacent one takes place abruptly as a function of the
chemical potential $\mu$.

It is necessary to compare the slave spins solution to the exact solution,
given by:

.. math::
   2N<d^\dagger d> =  \mathcal{Z}^{-1} \sum_{Q=0}^{2N} \binom{2N}{Q} Q e^{\beta(U/2(Q-N)^2 -\mu Q)}

As shown in the next plot, the slave spin approximation is capable of
recovering the coulomb occupation ladder, for the isolated atom with degenerate
fermions in spin and orbital. The approximation works best around half-filling.

.. plot::  degenerate_2orb_filling.py


Case: The lattice model - The Hubbard Model
'''''''''''''''''''''''''''''''''''''''''''

When in a lattice, atoms have overlapping orbitals and electrons are capable to
move along this lattice. Then for the Hamiltonian this term needs to be
included appearing in the Tight-Binding formulation. Then as simple extension
of the previous isolated atom case and in a multi-orbital scenario, the
Hamiltonian reads. The focus now for simplicity is the case of zero crystal-field splitting
$\epsilon_m=0$ and half-filling of each band one electron per site in each
orbital $\mu=0$.

.. math::
   \mathcal{H} = -\sum_m t_m \sum_{<i,j>, \sigma} (d^\dagger_{im\sigma}d_{jm\sigma} +h.c.)
    + \sum_{im\sigma}(\epsilon_m - \mu)d^\dagger_{im\sigma}d_{im\sigma}
    + \frac{U}{2} \sum_i \left( \sum_{m\sigma} d_{im\sigma}^\dagger d_{im\sigma} - N \right)^2
   :label: Hubbard_multiorb

Here it is needed to enforce the restriction:

.. math::  f_{im\sigma}^\dagger f_{im\sigma} = S_{im\sigma}^z + \frac{1}{2}
   :label: slavespinrestriction_multiorbitalsite

using the Lagrange multiplier $\lambda_{im\sigma}$, which can be used declaring
specific constrains to lattice site, orbital, and spin.

When rewriting the Hamiltonian in terms of the auxiliary fermions and the slave
spins the interaction term turns easily into:

.. math:: \mathcal{H}_{int} = \frac{U}{2} \sum_i \left( \sum_{m\sigma} S^z_{im\sigma} \right)^2

For the non interacting part, an appropriate representation of the creation
operator has to be chosen. The direct possibility $d^\dagger \rightarrow S^+ f^\dagger$,
although correct leads to problems with the spectral weight conservation because
$S^+$ and $S^-$ don't commute. Instead the representation $d^\dagger \rightarrow
2S^xf^\dagger$ and $d \rightarrow 2S^xf$ is chosen, which is identical on the physical Hilbert
space and involves commuting slave spin operators.

The constrain is treated on average using a static and
site, orbital and particle independent Lagrange multiplier $\lambda_{im\sigma}$.
Then the Hamiltonian reads:


.. math:: \mathcal{H} = &\frac{U}{2} \sum_i \left( \sum_{m\sigma} S^z_{im\sigma} \right)^2 \\
   &-\sum_m t_m \sum_{<i,j>, \sigma} 4S^x_{im\sigma}S^x_{jm\sigma}(f^\dagger_{im\sigma}f_{jm\sigma} +h.c.) \\
   &+ \sum_{im\sigma}(\epsilon_m - \mu)f^\dagger_{im\sigma}f_{im\sigma} \\
   &+\sum_{im\sigma} \lambda_{im\sigma}\left( S_{im\sigma}^z + \frac{1}{2} - f_{im\sigma}^\dagger f_{im\sigma} \right)

Using a Hartree-Fock approximation for the operators $S$ and $f$:

.. math::
   S^x_{im\sigma}S^x_{jm\sigma}(f^\dagger_{im\sigma}f_{jm\sigma} +h.c.)
   \approx <S^x_{im\sigma}S^x_{jm\sigma}>(f^\dagger_{im\sigma}f_{jm\sigma}
   +h.c.) \\
   +S^x_{im\sigma}S^x_{jm\sigma}<f^\dagger_{im\sigma}f_{jm\sigma} +h.c.> \\
   -<S^x_{im\sigma}S^x_{jm\sigma}(f^\dagger_{im\sigma}f_{jm\sigma} +h.c.)>

it is then possible to decouple the Hamiltonian into two effective ones:

.. math:: \mathcal{H}^f_{eff} = &-\sum_m t_m^{eff} \sum_{<i,j>, \sigma} (f^\dagger_{im\sigma}f_{jm\sigma} +h.c.) \\
   &+\sum_{im\sigma} (\epsilon_m - \mu - \lambda_{im\sigma}) f_{im\sigma}^\dagger f_{im\sigma}
   :label: hamileff_fermion

.. math:: \mathcal{H}^S_{eff} = &-\sum_m 4J^{eff}_m \sum_{<i,j>, \sigma} S^x_{im\sigma}S^x_{jm\sigma} \\
   &+\sum_{im\sigma} \lambda_{im\sigma} \left( S_{im\sigma}^z + \frac{1}{2} \right)
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

The fermion field Hamiltonian is a non-interacting one, and it's analytical
solution is well known. For the slave spin Hamiltonian, it can be treated
in a single-site using the Weiss mean field approximation.

.. math:: \mathcal{H}_s = &\sum_{m\sigma} 2h_mS^x_{m\sigma}
   +\sum_{m\sigma} \lambda_{m\sigma} \left( S_{m\sigma}^z + \frac{1}{2} \right)
   +\frac{U}{2} \left( \sum_{m\sigma} S^z_{m\sigma} \right)^2
   :label: hamil_spin_meanfield


Here the mean field $h_m$ has to be determined self-consistently from:

.. math::
    h_m \equiv -2zJ^{eff}_m<S^x_{m\sigma}> = 4<S^x_{m\sigma}>\frac{1}{N_s}\sum_k \epsilon_{km}<f^\dagger_{km\sigma}f_{km\sigma}>
    :label: mean_field

$z$ is the coordination number, $\epsilon_{km}=-t_m\sum_{\{\vec{a}\}}e^{-i\vec{k}\cdot\vec{a}}$
with $\{\vec{a}\}$ the set of vectors to the nearest neighbours

The effective fermion Hamiltonian is

.. math:: \mathcal{H}^f_{eff} = &\sum_{km\sigma} (-t_m^{eff} \sum_{\{\vec{a}\}} e^{-i\vec{k}\cdot\vec{a}} - \lambda_m) f^\dagger_{km\sigma}f_{km\sigma} \\
   &=\sum_{km\sigma} (Z_m\epsilon_{mk} + \epsilon_m - \mu - \lambda_{m\sigma}) f^\dagger_{km\sigma}f_{km\sigma}
   :label: 1site_quasipartH

where $Z_m=4<S^x_{im\sigma}>^2$ is the quasiparticle weight.

In the ordered spin basis $\{\ket{\uparrow\uparrow}, \ket{\uparrow\downarrow}, \ket{\downarrow\uparrow}, \ket{\uparrow\downarrow}\}$, where the spin labelling the operators are then

.. math::
   S^z_{\uparrow} = \frac{1}{2} \left[\begin{smallmatrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & -1 & 0\\0 & 0 & 0 & -1\end{smallmatrix}\right]
.. math::
   S^z_{\downarrow} = \frac{1}{2} \left[\begin{smallmatrix}1 & 0 & 0 & 0\\0 & -1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & -1\end{smallmatrix}\right]
.. math::
   S^x_{\uparrow} = \frac{1}{2} \left[\begin{smallmatrix}0 & 0 & 1 & 0\\0 & 0 & 0 & 1\\1 & 0 & 0 & 0\\0 & 1 & 0 & 0\end{smallmatrix}\right]
.. math::
   S^x_{\downarrow} = \frac{1}{2} \left[\begin{smallmatrix}0 & 1 & 0 & 0\\1 & 0 & 0 & 0\\0 & 0 & 0 & 1\\0 & 0 & 1 & 0\end{smallmatrix}\right]

As seen in equations :eq:`Hubbard_multiorb`
there is no hybridization between bands, hopping preserves then the orbital
quantum number.
When treating the system within a local mean field, in absence of hybridization
the $\vec{k}$ dependence enters the problem only through each band dispersion
as seen in equations :eq:`mean_field`, :eq:`1site_quasipartH`. Sums
over momenta can thus be replaced by integrals over the energy weighted by the
density of states $D(\epsilon)$, which is specific to the
lattice geometry and dimension. For this work, as commonly employed in the
literature, the Bethe lattice will be used. It has a very simple semi-circular
form of the density of states:

.. math::
    D(\epsilon) = \frac{1}{2 \pi t^{2}} \sqrt{4 t^{2} - \epsilon^{2}}
    :label: bethe_dos


and allows to simplify calculations in a great amount. Here $t$ is the
hopping amplitude and the half-bandwidth is $D=2t$, which is set as
the energy unit throughout this work. It is known moreover that the Bethe
lattice
well portrays the salient physical properties of the Mott-Hubbard transition and
it has immediate connection with the dynamical mean field
theory [Georges1996]_, which is exact in the infinite dimensions limit
and which we intend to implement in future
work.

The mean field in equation :eq:`mean_field` is then simplified into:

.. math::
    h_{m\sigma} = \langle O_{m\sigma} \rangle \int_{-\infty}^\infty
    \epsilon D(\epsilon)  n_F(Z_{m\sigma}\epsilon + \epsilon_m - \mu -
    \lambda_{m\sigma}) d\epsilon
    :label: DOS_meanfield

where $n_F$ is the Fermi distribution function. In the same fashion to estimate
the average particle number per site, orbital and spin, one easily uses the
relation:

.. math::
    \langle n_{im\sigma}\rangle  = \int_{-\infty}^\infty
    D(\epsilon)  n_F(Z_{m\sigma}\epsilon + \epsilon_m - \mu -
    \lambda_{m\sigma}) d\epsilon
    :label: DOS_avgparticles


In this work all calculations are done at zero temperature,
where the Fermi distribution can be approximated into a step function. That
implies for equations :eq:`DOS_avgparticles` and :eq:`DOS_meanfield`
that:

.. math::
    Z_{m\sigma}\epsilon_{F_0}(n) = - \epsilon_m + \mu + \lambda_{m\sigma}
    :label: fermi_energy

in which $\epsilon_{F_0}$ is the Fermi energy at zero temperature for the
non-interacting system such that

.. math::
    \int_{-\infty}^{\epsilon_{F_0}} D(\epsilon)  d\epsilon=n
    :label: fermi_energy_cut

This procedure of defining a zero temperature
non-interacting Fermi energy($U=0$ and thus $Z=1$) allows to keep the particle
population fixed when correlations are included into the
problem [Yu2011]_ [Florens2004]_.

.. plot::

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

Introducing doping
-----------------

The previous introduction to the slave spins treats the Hubbard Hamiltonian
and deals with the quartic term that deals with the correlations. But it is
only valid in the simplified case of degenerate orbitals at
half-filling. When aiming to introduce doping, a more elaborate formulation is
required. The fact is that the Spin Hamiltonian is unable to manage doping
cases as previously formulated, since it is always populated by 2 spin which
flip by the action of the operator $S^x$ in a balanced manner. This is
equivalent to the real electron hopping from one lattice site to the next.

In the case of electron or hole doping, the fermion Hamiltonian deals with it
through the chemical potential, which doesn't appear in the Spin Hamiltonian.
Then a new generic Spin operator needs to be introduced, to treat this spin
flip unbalance in the Spin system and deal with the doping cases. It is an
mixture of the spin raising and lowering operators.

.. math:: O^\dagger = S^+ + cS^-
   :label: generic_spin_op

The inclusion of the gauge parameter $c$ originates from the study of the
action of the real creation and annihilation operators $(d^\dagger, d)$ into
the
real states.  And then how the new operators act into the extended Hilbert
space. The reference case is such

.. math::
   :nowrap:

   \begin{align*}
     d\ket{\emptyset} &= 0 &   d\ket{1} &= \ket{0} \\
     d^\dagger\ket{1} &= 0 &   d^\dagger\ket{\emptyset} &= \ket{1}
   \end{align*}

The first conditions on the left can be assured by the fermion operator.

.. math::
   f O \ket{n^f=0, S^z = -\frac{1}{2}} = 0 \\
   f^\dagger O^\dagger \ket{n^f = 1, S^z = +\frac{1}{2}} = 0

the $O$ operator does not play any role. For second set of conditions

.. math::
   f O \ket{n^f=1, S^z = +\frac{1}{2}} = \ket{n^f = 0, S^z = -\frac{1}{2}} \\
   f^\dagger O^\dagger \ket{n^f = 0, S^z = -\frac{1}{2}} = \ket{n^f = 1, S^z =
   +\frac{1}{2}}

only three out of four matrix elements of the $O$ operator can be determined,
which implies that

.. math:: O = \left( \begin{matrix} 0 & c \\ 1 & 0 \end{matrix} \right)
   :label: generic_spin_op_matrix

Which is in correspondence  with the previous formulation of equation
:eq:`generic_spin_op`. $c$ can be an arbitrary complex number and it is tuned
in order to give rise to the most physical approximation scheme, by imposing
that it correctly reproduces solvable limits of the problem such as the
non-interacting limit.

The non-interacting limit
'''''''''''''''''''''''''

Starting with a Tight Binding Hamiltonian that only includes the kinetic energy term and
a contribution from the chemical potential to control the doping, the new
operators are introduced. Treating a single orbital case of atoms in a lattice.

.. math::
   \mathcal{H}_0 &= -t\sum_{<ij>, \sigma} (d^\dagger_{i\sigma} d_{j\sigma}
   +h.c.)
   - \mu\sum_{i\sigma} d^\dagger_{i\sigma}d_{i\sigma} \\
   \rightarrow \mathcal{H}_0 &= -t\sum_{<ij>, \sigma}
   (O^\dagger_{i\sigma}O_{i\sigma} f^\dagger_{i\sigma}f_{i\sigma} +h.c.) -
   \mu\sum_{i\sigma} f^\dagger_{i\sigma} f_{i\sigma} \\
   &+ \sum_{i\sigma} \lambda_i(S^z_{i\sigma} + \frac{1}{2} -
     f^\dagger_{i\sigma}f_{i\sigma})

Here the Lagrange multiplier is treated as individual to every site. Using
first a Hartree-Fock approximation in the tight binding term and respecting the
restriction set by the Lagrange multiplier, it is possible to separate the
Hamiltonian into 2 coupled effective ones.

.. math::
   \mathcal{H}_f &= -t\sum_{<ij>,\sigma}( Q_{ij}f^\dagger_{i\sigma}f_{i\sigma}
   +h.c.) - \sum_{i\sigma}(\mu + \lambda_i) f^\dagger_{i\sigma}f_{i\sigma} \\
   \mathcal{H}_s &= -\sum_{<ij>,\sigma} ( J_{ij}O^\dagger_{i\sigma}O_{j\sigma}
   +h.c.) + \sum_{i\sigma} \lambda_i(S^z_{i\sigma} + \frac{1}{2})

The parameters $Q_{ij}$ hopping renormalization factor, $J_{ij}$ slave-spin
exchange constant and $\lambda_i$ in these expressions are determined from the
coupled self-consistency equations:

.. math::
   Q_{ij} = \langle O^\dagger_{i\sigma}O_{i\sigma} \rangle
   :label: hopping_renorm

.. math::
   J_{ij} = t \langle f^\dagger_{i\sigma}f_{j\sigma} \rangle
   :label: ss_exchange

.. math::
   \langle n^f_{i\sigma} \rangle_f = \langle S^z_{i\sigma} \rangle_s +
   \frac{1}{2}
   :label: restriction

One further approximation needs to be applied, and it is to treat the spin
Hamiltonian within a Weiss mean field, in which a single site is embedded in an
effective field of its surroundings. The spin Hamiltonian becomes:

.. math::
   \mathcal{H}_S = \sum_\sigma (h_\sigma O^\dagger_\sigma + h.c.)
   +\sum_{\sigma} \lambda(S^z_{\sigma} + \frac{1}{2})

Where the mean field $h_\sigma$

.. math:: h_\sigma = \langle O_\sigma \rangle \frac{1}{N_s} \sum_k \epsilon_k
   \langle f^\dagger_{k\sigma}f_{k\sigma} \rangle


Choosing the gauge c
""""""""""""""""""""

The condition set, to reproduce the non-interacting case is $Q_{ij} = Z = 1$,
where the quasiparticle residum $Z=\langle O_\sigma \rangle$. In the single
site approximation $Q_{ij} = Z$ by construction, and only it's unitary value
remains to be enforced $Z=1$. The non-interacting spin Hamiltonian is treated
suppressing the spin index $\sigma$ since in this case up-spin and down-spin
fermions are decoupled.

.. math::
   \mathcal{H}_S &= hO^\dagger+\overline{h}O+\lambda(S^z_\sigma + \frac{1}{2})\\
   &= \begin{pmatrix} \lambda & c \overline{h} + h\\h \overline{c} +
   \overline{h} & 0 \end{pmatrix}

It is possible to diagonalize the Hamiltonian for one slave spin in the $S^z =
\pm 1/2$ basis. The ground state eigenvalue $E_{GS}$ and the corresponding
eigenstate are.

.. math::
   E_{GS} &= \frac{\lambda}{2} - \sqrt{\frac{\lambda^2}{4} + |a| ^2} \equiv
   \frac{\lambda}{2} - R \\
   \ket{GS} &= \begin{pmatrix} - \frac{a}{N} \\ \frac{\lambda /2+R}{N}
   \end{pmatrix}

Where $N=\sqrt{2R(\lambda /2 +R)}$ and $a=c \overline{h} + h$. The expected
values of $S^z$ and $O$ are:

.. math::

   \langle S^z \rangle = -\frac{\lambda}{4R} \\
   \langle O   \rangle = - \frac{a +\overline{a}c}{2R}

It is clearly seen that the Lagrange multiplier $\lambda$ depends on the
density $n$ and it is adjusted to satisfy the constraint equation:

.. math:: n-\frac{1}{2} = \langle S^z \rangle = -\frac{\lambda}{4R}
   :label: Sz_expected_constrain

and $c$ needs to be tuned to match the condition $Z=1$

.. math:: Z = \langle O \rangle^2 = \frac{|a +\overline{a}c| ^2 }{4R^2}=1
   :label: Zadjust

It is possible to eliminate $\lambda$ from the conditions by squaring
:eq:`Sz_expected_constrain` and using the relation :eq:`Zadjust`, following the next derivation:

.. math::
   \frac{|a| ^2}{4R^2} +(n-\frac{1}{2})^2 &= -\frac{\lambda^2}{16R^2} +
   \frac{|a| ^2}{4R^2} \\
   \frac{|a| ^2}{|a+\overline{a}c| ^2} &= n - n^2

Then it is possible to choose $c$ real, making also $h$ and $a$ real. The
expression for $c$ is found to be independent of the mean field $h$:

.. math:: c = \frac{1}{\sqrt{n(1-n)}} -1

.. plot::

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

    #band dop_phasediag
    def plot_dop_phase(bands, int_max, hund_cu):
        """Phase plot of Quasiparticle weight for N degenerate bands
        under doping shows transition only at integer filling
        the rest are metallic states"""
        name = 'Z_dop_phase_'+str(bands)+'bands_U'+str(int_max)+'J'+str(hund_cu)
        dop = np.sort(np.hstack((np.linspace(0.01,0.99,50),
                        np.arange(1./2./bands, 1, 1/2/bands))))
        data = ssplt.calc_z(bands, dop, np.arange(0, int_max, 0.1), hund_cu, name)

        ssplt.surf_z(data, name)

    plot_dop(1, 4.6, [0.5, 0.499, 0.495, 0.49, 0.45, 0.4, 0.2, 0.1], 0., 'Z_dop_1orb')
    plot_dop_phase(2, 6, 0.)

