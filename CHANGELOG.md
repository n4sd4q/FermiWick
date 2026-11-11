# Changelog

## Unreleased

- Gaussian states written by spin components, one amplitude for each contraction: $\langle c^\dagger_{i\sigma}c_{j\sigma'}\rangle$ =
  `χ[i, j, s, t]` and $\langle c_{i\sigma}c_{j\sigma'}\rangle$ = `Δ[i, j, s, t]` in the general state; `χ[i, j, s]`, and `Δ[i, j, s, t]`
  with `s ≠ t`, when $S^z$ is conserved; `χ[i, j]` and `Δ[i, j]` with SU(2): every symmetry removes indices, and the table of
  the general state reads at sight. The spin channels `χ[i, j, a]`, `Δ[i, j, a]` are gone (the magnetization is half the
  difference of the densities, the singlet and the triplet are combinations of the components). The results of states
  with less than SU(2) symmetry are written differently; `SpinSingletState` and every state invariant under spin rotations
  give the same expressions as before. `"Filling" -> ν` alone now sets $\langle n_{i\sigma}\rangle=\nu$ with no magnetization
  (before, with `"SpinRotationZ"`, $\chi^z_{ii}$ stayed free); `"Magnetization"` adds it.
- `GeneralState` (ρ, κ) is gone: `GaussianState[]` is the general state, also for modes with more labels, whose
  amplitudes group the labels by type with the spins last (`χ[i, j, o, o', s, t]`), and for spinless modes (`χ[i, j]`).
- Documents: the Gaussian states introduced from the general one, every contraction and nothing fixed, and its
  symmetries one at a time (the charge, $S^z$, spin rotations, time reversal), each removing contractions, in the wiki, the
  showcase, the manual and the guided tour; `SpinSingletState` presented as the shortcut for the state invariant under spin
  rotations at half filling, where that symmetry is explained; the on-site options after.
- Website: the link of the icon carries a hash of the file, so that browsers, which keep a favicon for long, fetch
  it again when the logo changes.

## 0.4.0 (2026-11-11)

- The adjoint of an operator defined with `DefineOperator` stays compact: `Dagger` gives the operator itself when it is
  Hermitian, the operator with the sites permuted when that is its adjoint (a hopping), otherwise
  `AdjointOperator[head][sites]`, registered automatically, displayed with a dagger, recognized by `ToSpinChargeBasis` and
  written by `ToTeX`. Operators of the user are written by `ToTeX` with their display.
- Fixed: `?GaussianState` (and `?` of other functions) in a notebook printed messages, because the display rules of
  $|x|^2$ and of a state evaluated the held code of the definitions; they now apply only to evaluated expressions.
- The logo of the banner shown at load draws both contractions, as the full logo does (the outer one in white from the
  creator on the right to the annihilator on the left, the inner one in the tint the other way), with the heavier lines of
  the small badge.
- Documents: exact diagonalisation done analytically in the showcase, the wiki and the video: the Hubbard dimer's
  eigenstates as kets, its free energy and entropy at finite temperature, and its Green's function as a Lehmann sum.
- Documents: a demo video (`Media/demo.mp4`, 1920×1080 at 60 fps, made by `Development/make_demo_video.py` from the first
  example of the README), its looping preview at the top of the README and the video under the hero of the website.
- In a notebook, the conjugate of an amplitude or of a symbol is shown as $x^\ast$ and the product $x^\ast x$ as
  $|x|^2$, as `ToTeX` writes them; only the display changes (copied, the output is the same expression).
- Documents: the bar of the website has a button, All packages, with the mark of the hub
  (https://n4sd4q.github.io/), beside GitHub.
- Documents: the README links the showcase and the manual on GitHub Pages, so that they open as PDFs in the browser,
  not on GitHub's page of the file (as the website did already, and now the hub too).
- Logo: the small badge (16 to 32 px: the icon of the website and of the hub) draws both contractions, as the full one,
  with heavier lines.
- Documents: the README and the website say what the package does in six lines, from the simplest to the heaviest, and
  build a first example up from the ladder operators to the Hubbard dimer solved exactly, the state you assume and the
  mean field. The guided tour goes the same way, one piece at a time, and ends with Wick's theorem checked against the
  exact ground state of free fermions. A wiki (`wiki/`, published by `.github/workflows/wiki.yml`): installation, writing
  operators, declaring a state, mean field with a self-consistent loop, exact diagonalisation, Gutzwiller projectors,
  LaTeX, scripts, conventions and a FAQ, every block checked by `Development/check_wiki.wl` (run by `check_docs.sh`).
- Documents: the picture of Wick's theorem follows the theme of the page (transparent, the ink in the colour of the
  text: `Development/make_readme_picture.py`); the wiki's page on Gutzwiller projectors shows projectors of your own and
  why the sandwiches are cheap.

## 0.3.0 (2026-11-10)

- The banner and the update check, as in FunctionalDiagrams: loaded in a notebook, the package prints a banner (version,
  license, how to cite it, links) and, at most once a day, in the background, asks GitHub for the latest release; if it is
  newer, the banner shows the line that installs it. `$FermiWickUpdateCheck = False` before `Needs` turns the check off,
  `FermiWickUpdateCheck[]` checks by hand; scripts print nothing and never check (`SECURITY.md`). Tested without the
  network on a saved answer of GitHub (`Tests/LatestRelease.json`).
- Installation in one line, from the GitHub release:
  `PacletInstall["https://github.com/n4sd4q/FermiWick/releases/download/v0.3.0/FermiWick-0.3.0.paclet"]`.
- Documents: a website, https://n4sd4q.github.io/FermiWick/ (`index.html`, written by `Development/make_site.py` from the
  README); the showcase and the manual in the common style of the author's packages (`Development/n4sd4q.sty`: the title
  band in the gradient of the package); a new logo in the common language (the same contractions on a badge, in the
  colours of the package), with a white and a small version.
- The first example of the README (and of the website) is one session: operators, contractions, the mean field of the
  Hubbard interaction, the exact spectrum of the Hubbard dimer, Wick's theorem written out.

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
