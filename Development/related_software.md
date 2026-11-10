# FermiWick and related software

Survey date: 2026-11-09. All claims about other tools come from the sources linked in the text, read on that date. "unverified" means that no source was found for the claim. The FermiWick column describes what the package is designed to do. Those features are claimed by the package and were not tested for this note.

## 1. Tools and status

| Tool | Language | Licence | Latest activity seen | Reference |
|---|---|---|---|---|
| **FermiWick** | Wolfram Language (no NCAlgebra) | to be decided | in development | — |
| SNEG | Wolfram Language | GPL (README: v2 or later; GitHub sidebar: GPL-3.0) | v2.1.1, commits Sep 2026, tested up to Mathematica 15 | Žitko, CPC 182, 2259 (2011), [arXiv:1105.4697](https://arxiv.org/abs/1105.4697); [repo](https://github.com/rokzitko/sneg) |
| QUANTUM (Gómez-Muñoz, Delgado) | Wolfram Language | "free"; licence not stated | paper 2016, Mathematica 8–10 | [J. Phys.: Conf. Ser. 698, 012019 (2016)](https://doi.org/10.1088/1742-6596/698/1/012019) |
| NCAlgebra | Wolfram Language | BSD | v6.0.5 (date not shown) | [ncalgebra.github.io](https://ncalgebra.github.io), [repo](https://github.com/NCAlgebra/NC) |
| Wolfram built-in `NonCommutativeAlgebra` family | Wolfram Language | proprietary | introduced 14.3 (2025), updated 15.0 | [ref page](https://reference.wolfram.com/language/ref/NonCommutativeAlgebra.html) |
| SymPy `physics.secondquant` | Python | BSD | SymPy 1.14.0 | [docs](https://docs.sympy.org/latest/modules/physics/secondquant.html) |
| QuantumAlgebra.jl | Julia | MIT | v1.7.1, Oct 2026 | [repo](https://github.com/jfeist/QuantumAlgebra.jl), [release notes](https://jfeist.github.io/QuantumAlgebra.jl/dev/release_notes) |
| QuantumCumulants.jl | Julia | MIT | active repository | [Quantum 6, 617 (2022)](https://quantum-journal.org/papers/q-2022-01-04-617/) |
| wick&d (`wickd`, formerly `wicked`) | C++/Python | MIT | v1.3.0, Sep 2026 | Evangelista, [JCP 157, 064111 (2022)](https://doi.org/10.1063/5.0097858); [repo](https://github.com/fevangelista/wickd) |
| SeQuant | C++ | LGPL-3.0 | release 26.2 (2026) | [JCP 164, 142502 (2026), arXiv:2511.09943](https://arxiv.org/abs/2511.09943) |
| pdaggerq | Python/C++ | Apache-2.0 | commits Nov 2026 | Rubin, DePrince, [Mol. Phys. 119, e1954709 (2021)](https://arxiv.org/abs/2106.06850) |
| drudge / gristmill | Python (on SymPy) | MIT | citation entry 2016–2025 | J. Zhao, PhD thesis, Rice (2018); [repo](https://github.com/tschijnmo/drudge), [docs](https://drudgecas.github.io/drudge/drudge_docs/0_9_0/api/algebras.html) |
| SpinAdaptedSecondQuantization.jl | Julia | MIT | v1.0 (2025) | [arXiv:2508.16342](https://arxiv.org/abs/2508.16342) |
| OpenFermion | Python | Apache-2.0 (licence not shown on the releases page; from memory, unverified) | v1.8.1 | [API symbols](https://quantumai.google/reference/python/openfermion/all_symbols) |

## 2. Feature matrix

Y = yes, P = partial or indirect (see the paragraph on the tool), N = no documented support, ? = could not be verified.

| Tool | Normal ordering (CAR) | Symbolic site/spin indices | Wick, normal contractions | Wick, anomalous contractions | User-declared contraction values / state | Compact composite spin/charge ops | Operator recognition (S·S etc.) | Mean-field decoupling, any order | Gutzwiller projector / approximation | Lattice fluxes | Matrices / ED |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **FermiWick** | Y | Y (spin ½ only) | Y | Y | Y | Y | Y (S·S, S·(S×S), user operators) | Y (incl. HFB, constant term) | Y | Y | Y (``FermiWick`Fock` ``: N, S_z and occupation sectors) |
| SNEG | Y | Y | Y (empty or Fermi-sea vacuum) | N | N | P | P | N | P (occupancy projectors only) | P (bond phases) | Y |
| QUANTUM | N | N | N | N | N | N | N | N | N | N | P (Dirac notation, qubits) |
| NCAlgebra | P (user must encode the CAR) | ? | N | N | N | N | N | N | N | N | N |
| WL built-in (14.3+) | P (relations or Gröbner reduction) | ? | N | N | N | N | N | N | N | N | N |
| SymPy secondquant | Y | P (no spin label) | Y (Fermi vacuum) | N | N | N | N | N | N | N | P (Fock states) |
| QuantumAlgebra.jl | Y | Y | P (cumulants) | P | P | N | N | P | N | N | N |
| QuantumCumulants.jl | ? (no fermions in README) | Y | P (cumulants) | ? | N | N | N | P (equations of motion) | N | N | N (ODE numerics) |
| wick&d | Y (general Fermi vacua) | Y (orbital spaces) | Y | N | P | N | N | N | N | N | N |
| SeQuant | Y | Y | Y | ? | N | N | N | N | N | N | numerical tensor evaluation |
| pdaggerq | Y (true or Fermi vacuum) | Y (orbital labels) | Y | N | N | N | N | N | N | N | einsum prototypes |
| drudge / gristmill | Y | Y | Y | ? | Y (`contractor`) | N | N | N | N | N | code generation |
| SpinAdaptedSQ.jl | Y | Y | Y | ? | N | P (spin-adapted operators) | N | N | N | N | N |
| OpenFermion | Y (integer modes only) | N | ? | N | N | N | N | N | N | ? | Y (sparse, Jordan–Wigner) |

## 3. Notes per tool

**SNEG.** This is the closest relative: a Mathematica package for fermionic second quantization, maintained since 2006. In its [source](https://github.com/rokzitko/sneg/blob/main/sneg.m), `wick`, `vevwick` and `vevwicknew` evaluate vacuum expectation values. The vacuum is EMPTY or SEA (`ordering`). Contractions are pairs of a creation and an annihilation operator that reduce to Kronecker deltas. No creator–creator pairs are generated, so there are no anomalous contractions. There is no mechanism for declaring symbolic contraction values such as χ_ij. SNEG has generators for many composite operators: `spinxyz`, `spinplus`, `spinminus`, `spinspin`, `isospin`, `chargecharge`, `hubbard`, `projector0/1/2`, `nambu`, `hop`, `hopphi`, `currentphi` and `anomaloushop`. These generators return fermion strings. `SnegSimplifySpin` and `SnegSimplifyNumber` then fold products back into `spinspin` and `number` using heuristic rules. No spin-chirality or user-defined pattern recognition was found. SNEG does more than FermiWick in several areas: bosons, Majorana fermions, Grassmann numbers, phonons, arbitrary spin, Dirac bras and kets, symbolic sums (`sum`, `sumSimplify`), and basis construction with matrix representations in (Q,S) or (Q,S_z) subspaces. It also underlies the NRG Ljubljana code. No mean-field or Gutzwiller routines were found.

**QUANTUM.** This add-on handles Dirac notation, operator algebra and commutators, quantum circuits, and quantum Hamilton dynamics. Its 2016 paper lists second quantization as future work ([paper](https://doi.org/10.1088/1742-6596/698/1/012019)). No later fermionic module was found (unverified).

**NCAlgebra.** A general non-commutative algebra suite: `**` multiplication, NCGB Gröbner bases and NCSDP semidefinite programming. Its applications are in control theory. It has no physics-specific rules. FermiWick has its own expansion of `**` and does not depend on it.

**Wolfram built-in (`NonCommutativeAlgebra`, `NonCommutativeExpand`, `NonCommutativeGroebnerBasis`, …).** The reference pages say these were introduced in 14.3 (2025), not 14.1. They provide generators, relations and Clifford, Grassmann and Weyl structures. The documentation contains no ladder-operator, normal-ordering or Wick examples, and symbolic-index families are not documented.

**SymPy `secondquant`.** Provides `F`/`Fd`/`B`/`Bd`, `NO`, `wicks` (with `keep_only_fully_contracted`), `contraction`, `substitute_dummies`, `AntiSymmetricTensor` and `evaluate_deltas`. Indices are dummies above or below the Fermi level. A contraction of two creators or two annihilators gives zero, and there are no thermal or user-defined states ([docs](https://docs.sympy.org/latest/modules/physics/secondquant.html)).

**QuantumAlgebra.jl.** Handles bosons, fermions, Pauli and two-level operators, and spins and so(n)/u(n) generators (marked experimental), with symbolic indices and sums. It provides `normal_form`, `expval`, `vacExpVal`, Heisenberg equations of motion and Julia code generation. `expval_as_corrs` writes ⟨…⟩ as a sum over partitions into connected correlators, which can include ⟨b b⟩-type ones. `dropcorr` and `droplen` truncate them. Dropping the connected parts beyond second order gives a Gaussian (Wick-like) factorization, but the documentation does not call it Wick's theorem. Spin operators are separate generators, not composites of fermions ([docs](https://jfeist.github.io/QuantumAlgebra.jl/dev/interface)).

**QuantumCumulants.jl.** Derives Heisenberg and Lindblad equations for bosonic and N-level systems, truncates them by a cumulant expansion of any order, and solves them numerically through ModelingToolkit. Its "mean field" means truncating an equation hierarchy, not decoupling a Hamiltonian. The README does not mention fermions.

**wick&d.** Implements Wick's theorem for general Fermi vacua, including multireference vacua, with up to eight user-defined orbital spaces. Contractions are expressed through the density and cumulant tensors of the reference. It derives CC, MBPT and MR equations with LaTeX output and handles fermions only. It does not document anomalous contractions.

**SeQuant.** Symbolic tensor algebra with a graph-based canonicalizer, a Wick engine for fermions and bosons, and lowering to numerical tensor evaluation ([arXiv:2511.09943](https://arxiv.org/abs/2511.09943)). Its focus is quantum chemistry. Anomalous contractions were not found (unverified).

**pdaggerq.** Normal-orders operator strings relative to the true vacuum or the Fermi vacuum. It covers fermions and photons (QED-CC), CC up to CCSDTQ and EOM-CC, and spin blocking. It generates C++ or Python and NumPy einsum code.

**drudge / gristmill.** A tensor algebra on SymPy with canonicalization under symmetries and symbolic summations. It has CAR, CCR and Clifford algebras and an su(2) lattice algebra (`SU2LatticeDrudge`). `WickDrudge` accepts a user-defined `contractor`, so contraction values can be declared. In `FockDrudge` only the annihilator–creator contraction can be customized (`ancr_contractor`), and the creator–creator rules are fixed. Gristmill optimizes and generates numerical code.

**SpinAdaptedSecondQuantization.jl.** Works directly with spin-adapted (spin-summed) excitation operators for coupled-cluster theory and is used in the eT program. This is the only other tool found that keeps spin-summed composites compact, but they are E_pq-type excitations, not local spin or charge operators.

**OpenFermion.** Numerical: `FermionOperator` with integer mode indices, `normal_ordered`, Hubbard-model builders, `s_squared_operator`/`sz_operator`, Jordan–Wigner, and sparse matrices. `mean_field_dwave` builds a fixed d-wave BCS mean-field model. It is not a decoupling routine. FQE has a `wick` module that was not assessed.

**Related tools with a different scope.** ADG ([CPC 240, 202 (2019)](https://arxiv.org/abs/1809.01187), GPLv3) generates and evaluates Bogoliubov MBPT diagrams, which do contain anomalous propagators, and HF-MBPT diagrams. Pymablock ([SciPost Phys. Codebases 50 (2025)](https://scipost.org/10.21468/SciPostPhysCodeb.50)) does symbolic or numerical Schrieffer–Wolff perturbation theory to any order, with second-quantized bosons, fermions and spins. H-wave ([arXiv:2308.00324](https://arxiv.org/abs/2308.00324)) and MeanFi ([PyPI](https://pypi.org/project/meanfi/)) are numerical Hartree–Fock codes and do not produce symbolic decoupled operators. PhysLean formalizes Wick's theorem in Lean 4 ([arXiv:2505.07939](https://arxiv.org/abs/2505.07939)). No software was found that computes Gutzwiller-projected expectation values symbolically. Searches returned only methods papers.

## 4. Where FermiWick adds something

These points rest on the evidence above and on FermiWick working as designed.

1. **Wick's theorem in a quasi-free state that the user declares.** FermiWick works with symbolic sites and spins, with normal (ρ, χ) and anomalous (κ, Δ) contractions, ready-made SU(2) presets, and selection or truncation of local and non-local contractions. No surveyed symbolic tool combines these. SNEG and SymPy contract only against a vacuum or Fermi sea. drudge allows a custom contractor, but its Fock algebra fixes the creator–creator rules. QuantumAlgebra.jl reaches anomalous Gaussian factorization only indirectly through cumulant truncation.
2. **Symbolic mean-field Hamiltonians.** FermiWick derives the Hartree–Fock(–Bogoliubov) Hamiltonian of an operator of any order from the variational principle, including the constant term. No surveyed package offers this as a function. Existing mean-field codes are numerical.
3. **Gutzwiller tools.** FermiWick handles P X P sandwiches, the Gutzwiller approximation to first order in inter-site contractions. No software doing this symbolically was found. SNEG has occupancy projectors but no Gutzwiller approximation.
4. **Compact site-local composites with recognition.** n, S, D, doublon and the projectors are multiplied as compact composites site by site. FermiWick then recognizes S_i·S_j, S_i·(S_j×S_k) and user-defined operators. SNEG comes closest, with heuristic `spinspin` folding but no chirality or user-defined pattern recognition.
5. **No dependencies.** FermiWick needs nothing but the Wolfram Language (no NCAlgebra).

## 5. Where others are stronger

- **Breadth of the algebra.** SNEG, QuantumAlgebra.jl, SeQuant, pdaggerq, drudge and SymPy also handle bosons. SNEG also handles Majorana fermions, Grassmann numbers and arbitrary spin S. QuantumAlgebra.jl also handles Pauli operators, qudits and Lie-algebra generators. FermiWick handles spin-½ fermions only.
- **Exact diagonalization inside one package.** SNEG builds bases and matrix representations in (Q,S) and (Q,S_z) subspaces and is used in production NRG. OpenFermion produces sparse matrices and maps to qubits. FermiWick writes its operators as sparse matrices on sectors of fixed N, S_z or occupations (``FermiWick`Fock` ``), but has no sectors of total S and no lattice symmetries.
- **Tensor, index and coupled-cluster machinery.** wick&d, SeQuant, pdaggerq, drudge/gristmill and SpinAdaptedSQ.jl canonicalize dummy-index tensor expressions and derive CC, EOM and multireference equations. They also generate optimized numerical code (einsum, C++, TiledArray-style backends) and scale to very large expressions. FermiWick does none of this.
- **Dynamics and diagrammatics.** QuantumAlgebra.jl and QuantumCumulants.jl derive equations of motion and cumulant hierarchies and connect them to ODE solvers. ADG generates diagrams of any order in Bogoliubov MBPT. Pymablock performs systematic Schrieffer–Wolff expansions to any order.
- **Maturity and access.** SNEG has had 20 years of use and a CPC paper. The Python, Julia and C++ tools are free to run, while FermiWick requires a Mathematica licence and has no publication yet.
