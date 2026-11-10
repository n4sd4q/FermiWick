<p align="center">
  <img src="Logo/FermiWick_logo.svg" width="120" alt="FermiWick logo: two Wick contractions of c c-dagger c c-dagger, drawn as brackets with arrows from the creator (open circle) to the annihilator (filled circle)">
</p>

<h1 align="center">FermiWick</h1>

<p align="center">
  <b>Wick's theorem, mean field, Gutzwiller projectors and exact diagonalisation for fermions</b><br>
  a Wolfram Language package for ladder operators, as symbols and as matrices
</p>

<p align="center">
  <a href="https://github.com/n4sd4q/FermiWick/releases/latest"><img src="https://img.shields.io/github/v/release/n4sd4q/FermiWick?label=release&color=8A2880" alt="latest release"></a>
  <img src="https://img.shields.io/badge/Mathematica-13.0%2B-8A2880" alt="Mathematica 13.0 or later">
  <a href="Tests"><img src="https://img.shields.io/badge/tests-85%20passed-2E7D32" alt="85 tests passed"></a>
  <img src="https://img.shields.io/badge/license-GPL--3.0--or--later-8A2880" alt="GNU GPL v3 or later">
</p>

<p align="center">
  <a href="https://n4sd4q.github.io/FermiWick/"><b>Website</b></a>&nbsp;·
  <a href="Showcase/FermiWick_Showcase.pdf"><b>Showcase</b></a>&nbsp;·
  <a href="Documentation/FermiWick_Manual.pdf"><b>Manual</b></a>&nbsp;·
  <a href="Examples.wl"><b>Guided&nbsp;tour</b></a>&nbsp;·
  <a href="#installation"><b>Installation</b></a>&nbsp;·
  <a href="CHANGELOG.md"><b>Changelog</b></a>
</p>

FermiWick does by computer the algebra one does on paper with fermion ladder operators, in the order of a calculation:

1. **Operators.** It brings products to normal order site by site, keeping the local operators — densities, spins,
   double occupancy, doublons, projectors — compact. It decomposes products of ladder operators into these local
   operators and recognizes the operators that act on several sites: $\mathbf S_i\cdot\mathbf S_j$, the scalar
   chirality $\mathbf S_i\cdot(\mathbf S_j\times\mathbf S_k)$ and operators of your own.
2. **Contractions.** It evaluates expectation values by Wick's theorem in a Gaussian state that you declare by its
   symmetries — spin rotations, rotations about z, charge, time reversal — or by rules for its contractions, with normal
   and anomalous contractions; `ContractionTable` shows which contractions survive.
3. **Mean field.** It gives the Hartree–Fock–Bogoliubov Hamiltonian of an operator of any order, with its constant.
4. **Gutzwiller projectors.** It handles projected operators and their expectation values in the Gutzwiller
   approximation.
5. **Matrices.** It writes any operator as a sparse matrix on the Fock space of a few sites, or on a sector of fixed
   particle number and spin, and finds ground states, expectation values and Green's functions: a small cluster is solved
   exactly, and a symbolic result checked.
6. **LaTeX.** It writes any result in LaTeX, and writes Wick's theorem out with the contractions drawn.

Sites and spins are symbols, so the same calculation holds on any lattice and for any filling; with numbers for the sites,
the same expressions become matrices.

## Installation

Install the paclet of the latest [release](https://github.com/n4sd4q/FermiWick/releases) once:

```mathematica
PacletInstall["https://github.com/n4sd4q/FermiWick/releases/download/v0.2.0/FermiWick-0.2.0.paclet"]
```

and load it in every session:

```mathematica
Needs["FermiWick`"]
Needs["FermiWick`Gutzwiller`"]   (* Gutzwiller projectors, when needed *)
Needs["FermiWick`Fock`"]         (* matrices and exact diagonalisation, when needed *)
```

There are no dependencies; the package needs Mathematica 13.0 or later and is tested with 15.0. The quickest start is the
guided tour, a notebook to run cell by cell, which comes with the paclet:

```mathematica
NotebookOpen[PacletObject["FermiWick"]["AssetLocation", "Examples"]]
```

(its source is [`Examples.wl`](Examples.wl), which also runs as a script). ``?FermiWick`*`` lists the functions and
`?Name` explains each one. To use a clone of the repository instead, without installing,
`PacletDirectoryLoad["/path/to/FermiWick/FermiWick"]` before `Needs`.

**Writing operators.** `c[i, s]` and `cDag[i, s]`, with `s` = `SpinUp`, `SpinDown` (or `"↑"`, `"↓"`, or a symbol) —
or any labels, `c[i]`, `c[i, orbital, s]`, `c[k]`, with `DeclareModeLabels` saying what they are;
`**` multiplies in the given order. The spin sum is `Sum[..., {s, $Spins}]`:

```mathematica
hopping = Sum[cDag[1, s] ** c[2, s] + cDag[2, s] ** c[1, s], {s, $Spins}];
```

## A first example

**Operators.** The commutator of two exchange terms that share a site, written back in spin operators, is the scalar
chirality; and the spin of a site is $\frac12$ only when the site is singly occupied:

```mathematica
Needs["FermiWick`"]

ToSpinChargeBasis[OperatorCommutator[SdotS[1, 2], SdotS[2, 3]]]
```

$$-i \mathbf S_{1}\cdot(\mathbf S_{2}\times\mathbf S_{3})$$

```mathematica
ToSpinChargeBasis[SdotS[1, 1]]
```

$$\frac{3 n_{1}}{4}-\frac{3 D_{1}}{2}$$

**Contractions.** A state invariant under spin rotations, with hopping $\chi$ and singlet pairing $\Delta$ on every
bond; in the exchange, hopping and pairing enter with the same weight:

```mathematica
state = SpinSingletState["Pairing" -> Δ];
WickVEV[SdotS[1, 2], state]
```

$$-\frac{3}{2} \Delta_{12} \Delta_{12}^{\ast}-\frac{3}{2} \chi_{12} \chi_{12}^{\ast}$$

**Mean field.** The mean-field Hamiltonian of $J\thinspace \mathbf S_1\cdot\mathbf S_2$ has a hopping field, a pairing field
and a constant, the coefficients of $c^\dagger_{1\uparrow}c_{2\uparrow}$, of $c^\dagger_{1\uparrow}c^\dagger_{2\downarrow}$ and of
$1$:

```mathematica
mf = MeanField[J SdotS[1, 2], state];
fields = OperatorCoefficients[mf];
{fields[cDag[1, "↑"] ** c[2, "↑"]], fields[cDag[1, "↑"] ** cDag[2, "↓"]], fields[1]}
```

$$\left\\{-\frac{3}{4} J \chi_{12}^{\ast},\frac{3}{4} J \Delta_{12},\frac{3}{2} J \Delta_{12} \Delta_{12}^{\ast}+\frac{3}{2} J \chi_{12} \chi_{12}^{\ast}\right\\}$$

The same works at any order. For the six-fermion chirality $-J_\chi\thinspace \mathbf S_1\cdot(\mathbf S_2\times\mathbf S_3)$
each bond receives a hopping from the path through the third site:

```mathematica
mfc = OperatorCoefficients[MeanField[-Jc SdotScrossS[1, 2, 3], state]];
{mfc[cDag[1, "↑"] ** c[2, "↑"]], mfc[cDag[1, "↑"] ** cDag[2, "↓"]]}
```

$$\left\\{\frac{3}{4} i J_\chi \Delta_{13} \Delta_{23}^{\ast}+\frac{3}{4} i J_\chi \chi_{23} \chi_{13}^{\ast},\frac{3}{4} i J_\chi \Delta_{13} \chi_{23}^{\ast}-\frac{3}{4} i J_\chi \Delta_{23} \chi_{13}^{\ast}\right\\}$$

**Matrices.** The same operators on the two-site sector with two electrons and $S^z=0$: the spectrum of the Hubbard
dimer, exact in $t$ and $U$:

```mathematica
Needs["FermiWick`Fock`"]
```

```mathematica
hop = Sum[cDag[1, s] ** c[2, s] + cDag[2, s] ** c[1, s], {s, $Spins}];
sector = FockSector[FockSpace[{1, 2}], "N" -> 2, "Sz" -> 0];
Eigenvalues[Normal[FockMatrix[-t hop + U (DOp[1] + DOp[2]), sector]]]
```

$$\left\\{0,U,\frac{1}{2} \left(U-\sqrt{16 t^2+U^2}\right),\frac{1}{2} \left(\sqrt{16 t^2+U^2}+U\right)\right\\}$$

**LaTeX.** `WickTeX` writes the expansion with the contractions drawn (LaTeX package `simpler-wick`):

```mathematica
WickTeX[cDag[1, SpinUp] ** cDag[2, SpinDown] ** c[3, SpinDown] ** c[4, SpinUp], state]
```

<p align="center"><img src="Documentation/readme_wick.svg" width="640" alt="Wick's theorem for c-dagger 1 up, c-dagger 2 down,
c 3 down, c 4 up in a pairing state: two pairings drawn as contractions, equal to Delta 34 Delta 12 star plus chi 14 chi 23"></p>

**Next.** The [showcase](Showcase/FermiWick_Showcase.pdf) goes through everything the package does, in this order:
local operators and their algebra, the spin–charge basis and the recognition of spin operators, operators of your own,
states declared by rules, magnetic and general states, the selection of contractions, the mean field of the Hubbard
interaction and of spin models, Gutzwiller projectors with the Brinkman–Rice transition and the enhancement of the
exchange, and the matrices of the operators: Hubbard clusters solved exactly, and Wick's theorem checked on them. The [manual](Documentation/FermiWick_Manual.pdf) describes every function. In a notebook, ``?FermiWick`*``
lists the functions and `?Name` explains each one.

## What it does

- **Modes with any labels.** `c[l1, l2, ...]`: sites, spins, orbitals, momenta, or labels of your own, declared with
  `DeclareModeLabels`; generic labels give Kronecker deltas that are kept and substituted.
- **Normal order, site by site.** Ladder operators with symbolic sites and spins. A product is separated site by site —
  every factor acts on one site and is even or odd, so the sites are reordered with the fermion sign only — and each site
  is reduced once. Local operators ($n$, $S^z$, $S^\pm$, $D$, $B$, $B^\dagger$, projectors) stay compact.
- **Local regrouping and the spin–charge basis.** Products of ladder operators are decomposed into the local operators of
  each site ($c^\dagger_{1\uparrow}c_{1\downarrow}c^\dagger_{2\downarrow}c_{2\uparrow}\to S^+_1S^-_2$), and the operators on several
  sites are recognized: $\mathbf S_i\cdot\mathbf S_j$ (also from its components, also with an XXZ remainder), the chirality
  with its antisymmetry, and operators you define.
- **The Hermitian conjugate.** `Dagger` reverses products, exchanges creators and annihilators and conjugates the
  coefficients (amplitudes complex, couplings real); `PlusHC` adds the "+ h.c." and `HermitianOperatorQ` checks a
  Hamiltonian.
- **States declared by their symmetries.** The general Gaussian state is written in spin channels — hopping χ, pairing
  Δ and their spin parts χ<sup>a</sup>, Δ<sup>a</sup> (magnetization, triplet pairing) — always Hermitian. Each symmetry
  sets channels to zero or makes them real: `GaussianState[{"SpinRotation", "TimeReversal"}]` is the state with real,
  spin-independent hopping. `ContractionTable` shows the contractions that survive. States can also be declared by rules,
  and contractions selected by any criterion (only local, no anomalous, ...) or truncated by their number between sites.
- **Mean field of any order.** The Hartree–Fock–Bogoliubov Hamiltonian, with its constant, of a product of 4, 6, ...
  fermions, in any of these states: two operators out in every way, the rest contracted. It has the same expectation value
  and the same first variation as the operator.
- **Gutzwiller projectors.** The variational Gutzwiller operator, projected operators, and expectation values to first
  order in the inter-site contractions (the Gutzwiller approximation).
- **Matrices and exact diagonalisation** (``FermiWick`Fock` ``). Any expression of the package becomes a sparse matrix
  on a Fock space, or on a sector of fixed particle number, $S^z$ or occupations, built directly in the sector; then
  ground states (dense, or Arnoldi with complete degenerate multiplets), expectation values, the Green's function at zero
  temperature in the Lehmann representation, $\mathbf S^2$, particle–hole and time reversal. A Hubbard ring of twelve
  sites at half filling (a sector of 853 776 states) is written as a matrix in about four seconds.
- **LaTeX export.** Any expression in LaTeX; Wick's theorem written out, pairing by pairing, with the contractions drawn
  at different heights when they overlap, their values and the result.

<details>
<summary>Everything it does, with the names of the functions</summary>

**Operators**
- `c[i, s]`, `cDag[i, s]`, `SpinUp`, `SpinDown`, `$Spins`; local operators `nOp`, `nTotOp`, `SzOp`, `SpOp`, `SmOp`, `DOp`, `B`, `Bdag`, `SpinOperator`;
  `SdotS`, `SdotScrossS`, `EmptyProjector`, `SingleOccupationProjector`; operators of your own with `DefineOperator`
- products `**` and `·` (normal-ordered), `NormalOrder`, `OperatorCommutator`, `OperatorAnticommutator`,
  `ToSpinChargeBasis`, `ToLadderBasis`, `GroupLocalOperators`, `OperatorCoefficients`, `CleanDeltas`
- `Dagger`, `PlusHC`, `HermitianOperatorQ`

**States, Wick and mean field**
- `GaussianState` (from symmetries `"SpinRotation"`, `"SpinRotationZ"`, `"Charge"`, `"TimeReversal"`, or from rules),
  `ContractionTable`, `SpinSingletState`, `GeneralState`, `StateSelect` with `LocalPairQ` and `AnomalousPairQ`,
  `Contraction`, `$State`
- `WickVEV` (options `"MaxOrder"`, `"NormalOrder"`), `VerifyNormalOrder`, `MeanField`

**LaTeX**
- `ToTeX`, `WickPairings`, `WickTeX`

**Gutzwiller** (``FermiWick`Gutzwiller` ``)
- `GutzwillerProjector`, `LocalSandwich`, `GutzwillerVEV`, `GutzwillerSimplify`

**Lattices**
- `DefinePeierlsHopping`, `DefineFlux`

**Matrices** (``FermiWick`Fock` ``)
- `FockSpace`, `FockSector`, `ModeIndex`, `FockMatrix`, `FockExpectation`, `OccupationBasis`, `BasisState`, `StateToKet`
- `LowestStates`, `LehmannGreen`, `GreenMatrix`, `SpinSquared`, `ParticleHoleOperator`, `TimeReversal`, `SectorBlock`
- `ClearFockCache`, `$FockFullLimit`, `$FockDenseLimit`

</details>

**Scope.** FermiWick computes exact symbolic expressions for spin-½ fermions, and their matrices on clusters of up to
about a million states. Bosons and other flavours are outside its scope, and so are the lattice symmetries (momentum,
point group) that dedicated codes use to diagonalise larger clusters.

## Documentation

| Document | Contents |
|---|---|
| [Website](https://n4sd4q.github.io/FermiWick/) | the package at a glance: what it does, the first example, installation, the documents |
| [Showcase](Showcase/FermiWick_Showcase.pdf) (PDF) | the package step by step, from the anticommutation relations to the parton mean field of a chiral spin model, the Brinkman–Rice transition and Hubbard clusters solved exactly |
| [Manual](Documentation/FermiWick_Manual.pdf) (PDF) | every function, option and convention, with examples |
| [Guided tour](Examples.wl) | a notebook to run cell by cell, in the paclet (`FermiWick/Examples.nb`) |
| [Changelog](CHANGELOG.md) | what changed in each version |

Every output printed in the showcase and in the manual is computed by the code shown before it
(`Development/run_examples.wl`).

<details>
<summary>Conventions in brief</summary>

* `c[i, s]` annihilates a fermion on site `i` with spin `s` (`"↑"`, `"↓"` or a symbol). Different symbolic sites are
  different sites.
* `a ** b` is the product in the given order; `a · b` (`\[CenterDot]`) is the product in normal order, with every
  `cDag` to the left of every `c`.
* Contractions in spin channels: $\langle c^\dagger_{i\sigma}c_{j\sigma'}\rangle=\chi_{ij}\delta_{\sigma\sigma'}+\sum_a\chi^a_{ij}\sigma^a_{\sigma'\sigma}$,
  $\langle c_{i\sigma}c_{j\sigma'}\rangle=[(\Delta_{ij}+\sum_a\Delta^a_{ij}\sigma^a)\epsilon]_{\sigma\sigma'}$, $\epsilon_{\uparrow\downarrow}=1$;
  $\chi_{ji}=\chi_{ij}^{\ast}$, $\chi^a_{ii}=\langle S^a_i\rangle$, $\Delta_{ji}=\Delta_{ij}$, $\Delta^a_{ji}=-\Delta^a_{ij}$.

</details>

## Checks

The test suite (`math -script Tests/RunTests.wls`, 85 tests) checks the algebra, the mean-field formulas of the two-
and three-spin terms and of the Hubbard interaction, and the states. It checks the symbolic engine against the matrices:
the matrix of an expression equals that of its normal order and of its spin–charge form; `WickVEV` and `MeanField` give
the exact values in random Gibbs states of quadratic Hamiltonians with pairing, on three sites and with two orbitals per
site, and in the free ground state of a ring; the symmetries of the states hold on exact Gibbs states of symmetric
Hamiltonians. The matrices themselves are checked against the Jordan–Wigner construction, and exact diagonalisation on the
Hubbard dimer, a magnetic plaquette and Hubbard rings.
`sh Development/check_docs.sh` checks the documents: every example of the showcase and of the manual is run again in a
fresh kernel and must give the output printed, the guided tour must run without messages, and the code of the README must
give the formulas shown after it.

## Related software

[SNEG](https://github.com/rokzitko/sneg) (Žitko, Comput. Phys. Commun. **182**, 2259 (2011)) is the closest relative: a
Mathematica package for second quantization, with normal ordering, vacuum expectation values in the empty state or the
Fermi sea, composite spin and charge operators, and much more: bosons, Majorana fermions, arbitrary spin and exact
diagonalization. Its Wick contractions pair a creator with an annihilator.
[QuantumAlgebra.jl](https://github.com/jfeist/QuantumAlgebra.jl) (Julia) handles bosons and fermions, equations of motion
and cumulant expansions. [wick&d](https://github.com/fevangelista/wickd), [SeQuant](https://github.com/ValeevGroup/SeQuant),
[pdaggerq](https://github.com/edeprince3/pdaggerq), [drudge](https://github.com/tschijnmo/drudge) and SymPy's
[secondquant](https://docs.sympy.org/latest/modules/physics/secondquant.html) are built for quantum chemistry: index
and tensor algebra and the generation of coupled-cluster equations.

FermiWick adds four things to these:
- Wick's theorem in a Gaussian state that you declare, with normal and anomalous contractions and a selection of them;
- the symbolic mean-field Hamiltonian of an operator of any order;
- Gutzwiller-projected expectation values;
- the same operators as matrices on symmetry sectors, so that every symbolic result can be checked exactly on a cluster.

It keeps local operators compact and recognizes multi-site ones in the results.

## Citing

If you use FermiWick in published work, please cite it as described in [`CITATION.cff`](CITATION.cff).

## License

Released under the [GNU General Public License, version 3 or later](LICENSE), with two additional terms
([`NOTICE`](NOTICE)): the attribution to the author must be kept, and modified versions must say that they are modified,
by whom and when. Author: Nazzareno Africani (SISSA).
