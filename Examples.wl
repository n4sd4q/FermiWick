(* ::Package:: *)

(* ::Title:: *)
(*FermiWick*)

(* ::Subtitle:: *)
(*A guided tour: fermion operators, Wick's theorem, mean field*)

(* ::Section:: *)
(*0. Welcome*)

(* ::Text:: *)
(*FermiWick does by computer the algebra one does on paper with fermion ladder operators: products brought to normal order, local operators (densities, spins, double occupancy) kept compact and recognized in the results, expectation values by Wick's theorem in a Gaussian state that you declare, the mean-field Hamiltonian of an operator of any order, Gutzwiller projectors, the same operators as matrices for exact diagonalisation, and LaTeX of everything.*)

(* ::Text:: *)
(*How to use this notebook: evaluate the cells in order (Shift+Enter, or Evaluation > Evaluate Notebook); the whole tour takes a few seconds. The comment at the end of each input cell says what you should get. The same file runs as a script: math -script Examples.wl. The showcase (Showcase/FermiWick_Showcase.pdf) goes through everything in detail, and the manual describes every function.*)

(* ::Subsection:: *)
(*Loading the package*)

(* the paclet next to this file (the repository), otherwise the installed one *)
With[{folder = Quiet[If[$InputFileName =!= "", DirectoryName[$InputFileName], NotebookDirectory[]]]},
  If[StringQ[folder],
   Scan[If[FileExistsQ[FileNameJoin[{#, "PacletInfo.wl"}]], PacletDirectoryLoad[#]] &, {FileNameJoin[{folder, "FermiWick"}], folder}]]];
Needs["FermiWick`"]; Needs["FermiWick`Gutzwiller`"]; Needs["FermiWick`Fock`"];

$FermiWickVersion
(* the version of the package, e.g. "0.2.0" *)

(* ::Section:: *)
(*1. Writing operators*)

(* ::Text:: *)
(*c[i, s] annihilates a fermion on site i with spin s, cDag[i, s] creates it; the spin is SpinUp, SpinDown or a symbol, and sites can be symbols. ** multiplies in the given order; NormalOrder moves every cDag to the left.*)

NormalOrder[c[1, SpinUp] ** cDag[1, SpinUp]]
(* 1 - cDag[1, up] ** c[1, up]: the anticommutation relation *)

NormalOrder[c[i, SpinUp] ** cDag[j, SpinUp]]
(* -cDag[j, up] ** c[i, up]: different symbolic sites are different sites *)

hopping = Sum[cDag[1, s] ** c[2, s], {s, $Spins}];
PlusHC[t hopping]
(* the hopping between sites 1 and 2 plus its Hermitian conjugate, four terms *)

(* ::Section:: *)
(*2. Local operators and spins*)

(* ::Text:: *)
(*Densities, spins, double occupancy (DOp), doublons (B, Bdag) and projectors stay compact; ToSpinChargeBasis writes any expression with them and recognizes the spin operators on several sites.*)

ToSpinChargeBasis[SpOp[1] \[CenterDot] SmOp[1]]
(* n/2 + Sz - D on site 1: the algebra of one site *)

ToSpinChargeBasis[SdotS[1, 1]]
(* 3/4 n - 3/2 D: a spin 1/2 only on singly occupied sites *)

ToSpinChargeBasis[OperatorCommutator[SdotS[1, 2], SdotS[2, 3]]]
(* -I S1.(S2 x S3): the scalar chirality *)

(* ::Section:: *)
(*3. Gaussian states and Wick's theorem*)

(* ::Text:: *)
(*A Gaussian state is fixed by its contractions; declare it by its symmetries, and ContractionTable shows what survives. WickVEV gives expectation values.*)

state = GaussianState[{"SpinRotation", "Charge"}, "Filling" -> 1/2];
MatrixForm[ContractionTable[state, {1, 2}, "Normal"]]
(* hopping chi_12 between the sites, 1/2 on each site, nothing between opposite spins *)

WickVEV[SdotS[1, 2], state]
(* -3/2 |chi_12|^2 *)

pairing = GaussianState["SpinRotation", "Filling" -> 1/2, "OnSitePairing" -> 0];
WickVEV[SdotS[1, 2], pairing]
(* -3/2 (|chi_12|^2 + |Delta_12|^2): hopping and singlet pairing enter with the same weight *)

(* ::Section:: *)
(*4. Any labels, and states by rules*)

(* ::Text:: *)
(*A state can be declared by rules for its contractions, with any spin structure: here a hopping that flips the spin.*)

flip = GaussianState[{{cDag[i_, s_], c[j_, t_]} :> If[i === j, KroneckerDelta[s, t]/2, If[s === t, g, h]]}];
WickVEV[SzOp[1] ** SzOp[2], flip]
(* (h^2 - g^2)/2: the spin-flip hopping h enters with the opposite sign *)

(* ::Text:: *)
(*Modes can carry any labels: DeclareModeLabels says what each label is. With orbitals, the anticommutator gives a delta of the orbitals:*)

DeclareModeLabels[{"Site", "Orbital", "Spin"}];
NormalOrder[c[1, a, SpinUp] ** cDag[1, b, SpinUp]]
(* KroneckerDelta[a, b] - cDag[1, b, up] ** c[1, a, up] *)

(* ::Section:: *)
(*5. Mean field*)

(* ::Text:: *)
(*MeanField gives the Hartree-Fock-Bogoliubov Hamiltonian of an operator of any order: two operators kept out in every way, the rest contracted. The fields are its coefficients.*)

MeanField[U DOp[1], GaussianState["SpinRotation", "Filling" -> \[Nu], "OnSitePairing" -> \[Eta]]]
(* Hartree terms U nu n and BCS terms with eta, plus the constant *)

fields = OperatorCoefficients[MeanField[-Jc SdotScrossS[i, j, k], pairing]];
fields[cDag[i, SpinUp] ** c[j, SpinUp]]
(* the hopping field of the bond i -> j: (3 I/4) Jc (chi_jk chi_ki + Delta*_jk Delta_ki), with chi_ki = Conjugate[chi_ik] *)

(* ::Section:: *)
(*6. Gutzwiller projectors*)

real = GaussianState[{"SpinRotation", "Charge", "TimeReversal"}, "Filling" -> 1/2];
hop2 = Sum[cDag[1, s] ** c[2, s] + cDag[2, s] ** c[1, s], {s, $Spins}];
Simplify[GutzwillerVEV[LocalSandwich[GutzwillerProjector, {1, 2}, hop2], real]/WickVEV[hop2, real]]
(* Sin[2 theta]^2: the Gutzwiller reduction of the hopping, 8 d (1 - 2 d) with d = Sin[theta]^2/2 *)

(* ::Section:: *)
(*7. Matrices and exact diagonalisation*)

(* ::Text:: *)
(*FermiWick`Fock` writes any operator as a sparse matrix on the Fock space of a few sites, or on a sector with fixed N and Sz: the modes are the ladder operators, numbered site by site with the spin fastest.*)

sector = FockSector[FockSpace[{1, 2}], "N" -> 2, "Sz" -> 0];
Eigenvalues[Normal[FockMatrix[-t hop2 + U (DOp[1] + DOp[2]), sector]]]
(* the Hubbard dimer: 0, U and (U -+ Sqrt[U^2 + 16 t^2])/2 *)

{e, v} = LowestStates[FockMatrix[-hop2 + 4 (DOp[1] + DOp[2]), sector], 1];
FockExpectation[SdotS[1, 2], First[v], sector]
(* about -0.64: the two spins of the ground state form mostly a singlet *)

(* ::Section:: *)
(*8. LaTeX*)

ToTeX[ToSpinChargeBasis[OperatorCommutator[SdotS[1, 2], SdotS[2, 3]]]]
(* "-i \mathbf S_{1}\cdot(\mathbf S_{2}\times\mathbf S_{3})" *)

WickTeX[cDag[1, SpinUp] ** c[2, SpinUp] ** cDag[2, SpinUp] ** c[1, SpinUp], state]
(* Wick's theorem written out, with the contractions drawn (LaTeX package simpler-wick) *)
