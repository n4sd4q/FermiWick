# Ideas and things to do

A working list, roughly by value; strike or remove an item when it is done (the CHANGELOG records what was done).

## Showcase and documents

- **A teaser page** at the start of the showcase (and of the README): three results in a few lines each, before any
  explanation: Wick's theorem drawn, the t–J model derived from the Hubbard model, a spectral function of a cluster.
- **The t–J model from the Hubbard model**, the most iconic symbolic result missing: second-order perturbation theory in
  t/U with the projectors on singly occupied sites (`SingleOccupationProjector`, `DOp`), the result recognized by
  `ToSpinChargeBasis` as J (S·S − n n/4) with J = 4t²/U; then the three-site terms, and the ring exchange at fourth order
  as a stretch goal. Possibly a function of its own (see below).
- **Pictures.** The showcase has formulas only. Candidates: the spectral function A(ω) of a ring of four or six sites as
  U grows (a waterfall, from `LehmannGreen` with a broadening), the entropy and the specific heat of the dimer against T
  (from `PartitionFunction`), mean field against exact for the energy of a ring.
- **Mean field against exact**, one page: the Hartree–Fock and the BCS mean field of a small cluster, solved
  self-consistently, compared with exact diagonalisation of the same cluster (energy, double occupancy, spin
  correlations): what the package is for, in one figure.
- **Correlation functions at a temperature**, once they exist (below): the spin susceptibility of the dimer, χ(T) with
  the Curie law at high T and the singlet gap at low T.
- A short **"recipes" page in the wiki**: the ten most common tasks as copyable blocks (a Hamiltonian, its spectrum, a
  thermal average, a Green's function, a mean-field Hamiltonian, LaTeX of a result).

## Features

- **Thermal expectation values**: `ThermalExpectation[H, space, T, μ, op]` = Tr(e^{-(H-μN)/T} op)/Z, exact with symbols
  (the levels and projectors of `LehmannGreen` already give it); today `FockExpectation` takes a single vector.
- **Two-particle correlation functions** in the Lehmann representation, `LehmannCorrelation[H, space, A, B]`: the spin
  and charge susceptibilities, the dynamical structure factor, at zero and finite temperature, with the same options as
  `LehmannGreen`.
- **Modes as linear combinations** in `LehmannGreen` (and `FockMatrix`): `(c[1, ↑] + c[2, ↑])/Sqrt[2]` as a mode, so that
  the Green's function of a bonding orbital, or of a momentum, needs no change of basis by hand.
- **Spectral functions**: `SpectralFunction[g, ω, η]` (−Im G/π with a broadening) and the Matsubara Green's function on a
  list of frequencies, ready to plot.
- **From a quadratic Hamiltonian to a Gaussian state**: diagonalise the output of `MeanField` (Bogoliubov–de Gennes
  matrix, quasiparticle energies) and return the `GaussianState` of its ground or Gibbs state, with numbers. This closes
  the loop MeanField → state → MeanField and makes a built-in self-consistent solver possible
  (`SelfConsistentMeanField`), which the wiki now writes by hand.
- **Effective Hamiltonians**: `EffectiveHamiltonian[H0, V, P, order]` by Schrieffer–Wolff (or degenerate perturbation
  theory) for a projector P, recognized in the spin–charge basis (the t–J model above as its example).
- **Lattices**: a small builder for chains, ladders and square or triangular clusters with periodic boundaries and
  fluxes (the hopping lists that the documents write by hand), and Fourier modes c_k for rings.
- More symmetries of `GaussianState`: translation (amplitudes that depend on the distance only) and inversion.

## Performance

- `FockMatrix` of a power of a local operator is much slower than of its expanded form: on a ring of eight sites
  `2 Sum[(nTotOp[i] - 1)^2, {i, 8}]` takes 0.25 s against 0.045 s for `4 Sum[DOp[i], {i, 8}]` (the same operator up to a
  constant and a number term). Reduce local polynomials site by site (`NormalOrder` of the local factor, or a cache of
  the local matrices) before building the matrix.
- `LehmannGreen` with symbols at a finite temperature returns residues that need a `FullSimplify` of tens of seconds to be
  read; group the residues by symmetry (spin, the orbitals of a point group) and simplify each level once.
- Benchmark script, as in FunctionalDiagrams (`Development/benchmark.wl`): the timings of twenty representative
  calculations, each result compared with a reference, so that a faster kernel is checked to give the same results.

## Repository

- Continuous integration: run `Tests/RunTests.wls` on GitHub with the Wolfram Engine (the free engine needs an
  activation key as a secret); today the tests run only locally.
- Publish the paclet on the Wolfram Paclet Repository, so that `PacletInstall["FermiWick"]` works without the URL.
- A DOI for each release (Zenodo's GitHub integration), quoted in `CITATION.cff`.
- `CONTRIBUTING.md` and issue templates (a bug report asks for `$FermiWickVersion`, `$Version` and a minimal input).

## Rough edges seen while writing the documents

- `FullSimplify` does not always give the same form in two sessions (its caches depend on what was computed before): an
  output shown in the showcase and in the video, which run different sessions, must use `Simplify`, `Factor` or
  `Together`.
- `LehmannGreen` at zero temperature with symbols fails when the ground state depends on their values; the message says
  so, but the assumptions needed (here |μ| below the gap) are not obvious: say in the message which levels compete.
