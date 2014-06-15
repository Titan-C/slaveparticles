Slave Spins Mean Field Theory
=============================

The previous introduction to the slave spins, slave particle formulation, to
treat the Hubbard Hamiltonian and deal with the quartic term that deals with
the correlations is only valid in the simplified case of degenerate orbitals at
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

The Hubbard Model
-----------------

Through this work all developments are performed around the Hubbard
Hamiltonian.

The non-interacting limit
'''''''''''''''''''''''''

Starting with a free Hamiltonian that only includes the kinetic energy term and
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



