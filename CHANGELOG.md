# Changelog

## Unreleased

- Installation in one line, from the GitHub release:
  `PacletInstall["https://github.com/n4sd4q/FermiWick/releases/download/v0.2.0/FermiWick-0.2.0.paclet"]`.
- Documents: a website, https://n4sd4q.github.io/FermiWick/ (`index.html`, written by `Development/make_site.py` from the
  README); the showcase and the manual in the common style of the author's packages (`Development/n4sd4q.sty`: the title
  band in the gradient of the package); a new logo in the common language (the same contractions on a badge, in the
  colours of the package), with a white and a small version.

## 0.2.0 (2026-11-10)

- ``FermiWick`Fock` ``: the operators of FermiWick as matrices, and exact diagonalisation.
  `FockSpace` (modes `c[i, s]`, `c[i, o, s]` with `"Orbitals"`, spinless `c[i]` with `"Spins" -> None`), `FockSector` (fixed
  `"N"`, `"Sz"`, `"Counts"`), `ModeIndex`; `FockMatrix[expr, target]` reads any expression of the package, products in the
  order written, symbolic coefficients, with no normal ordering; `FockExpectation`, `OccupationBasis`, `BasisState`,
  `StateToKet`; `LowestStates` (dense or Arnoldi, complete degenerate multiplets), `LehmannGreen`, `GreenMatrix`;
  `SpinSquared`, `ParticleHoleOperator`, `TimeReversal`, `SectorBlock`; `ClearFockCache`, `$FockFullLimit`,
  `$FockDenseLimit`.
- Operators are built directly in a sector, which is indexed through the combinatorial number system, read 12 bits at a
  time from tables; diagonal strings need no look-up of their images. A Hubbard ring of twelve sites at half filling (a
  sector of 853 776 states) is written as a matrix in 3.8 s.
- Fixed: a product with a vanishing local factor (`c[2, "↑"] ** c[2, "↑"]` inside it) gave `2*0**0` instead of `0` in
  `ToSpinChargeBasis` when its coefficient was not 1 (found by comparing the matrices of expressions with those of their
  spin–charge forms).
- Tests: `Tests/Fock.wlt` (30 tests, 85 in all): the Jordan–Wigner matrices; the matrix of an expression equals that of its
  normal order and of its spin–charge form; Wick's theorem exact in a free ground state; sectors; the
  Hubbard dimer, a magnetic plaquette, Hubbard rings with Arnoldi. The exact checks of `WickVEV` and `MeanField` use
  `FockMatrix` instead of matrices written by hand.
- Documents: a section on matrices and exact diagonalisation in the showcase and in the manual, part 7 of the guided tour,
  the README.

## 0.1.0 (2026-11-09)

First version.

- Modes with any labels, `c[l1, l2, ...]`, declared with `DeclareModeLabels` (site, spin, generic labels); the general
  state `GeneralState` for any labels.
- Gaussian states (`GaussianState`) declared by their symmetries or by rules.
- Guided tour `Examples.wl` (`FermiWick/Examples.nb` in the paclet); `Development/check_docs.sh` checks the documents.
- License: GNU GPL v3 or later, with the additional terms of `NOTICE`.

- Ladder operators with symbolic sites and spins; normal ordering done site by site, with compact local operators
  (`nOp`, `nTotOp`, `SzOp`, `SpOp`, `SmOp`, `DOp`, `B`, `Bdag`, projectors); spin–charge basis with recognition of
  `SdotS`, `SdotScrossS` and of operators defined with `DefineOperator`; `OperatorCoefficients`.
- Gaussian states as objects, declared by their symmetries (`GaussianState[{"SpinRotation", "SpinRotationZ", "Charge",
  "TimeReversal"}]`, written in spin channels χ, χ^a, Δ, Δ^a, always Hermitian) or by rules; `ContractionTable`;
  `SpinSingletState`, `GeneralState`, `StateSelect`.
- `Dagger`, `PlusHC`, `HermitianOperatorQ`.
- The recognition of operators keeps the fermion sign when the sites are out of order.
- `WickVEV` (normal and anomalous contractions, `"MaxOrder"`), `MeanField` of any order.
- ``FermiWick`Gutzwiller` ``: `GutzwillerProjector`, `LocalSandwich`, `GutzwillerVEV`.
- `DefinePeierlsHopping`, `DefineFlux`.
- LaTeX: `ToTeX`, `WickPairings`, `WickTeX` (Wick's theorem written out, contractions drawn with simpler-wick).
- `SpinUp`, `SpinDown`, `$Spins`; `OperatorCoefficients[expr, monomial]`.
- Showcase, manual, guided tour (`Examples.nb`), 39 tests.
