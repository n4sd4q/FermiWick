(* ::Package:: *)

(* ::Title:: *)
(*FermiWick*)

(* ::Subtitle:: *)
(*A guided tour: from the ladder operators to a solved model*)

(* ::Section:: *)
(*0. Welcome*)

(* ::Text:: *)
(*FermiWick lets you write fermion operators as on paper and does the algebra for you. This tour starts from the ladder operators and builds everything up, one piece at a time: local operators and spins, a Hamiltonian, its exact solution, expectation values by Wick's theorem in the state you assume, mean field, and at the end a symbolic result checked against the exact one.*)

(* ::Text:: *)
(*How to use this notebook: evaluate the cells in order (Shift+Enter, or Evaluation > Evaluate Notebook); the whole tour takes a few seconds. The comment under each input says what you should get. The same file runs as a script: math -script Examples.wl.*)

(* ::Subsection:: *)
(*Loading the package*)

(* the paclet next to this file (the repository), otherwise the installed one *)
With[{folder = Quiet[If[$InputFileName =!= "", DirectoryName[$InputFileName], NotebookDirectory[]]]},
  If[StringQ[folder],
   Scan[If[FileExistsQ[FileNameJoin[{#, "PacletInfo.wl"}]], PacletDirectoryLoad[#]] &, {FileNameJoin[{folder, "FermiWick"}], folder}]]];
Needs["FermiWick`"]; Needs["FermiWick`Fock`"]; Needs["FermiWick`Gutzwiller`"];

$FermiWickVersion
(* the version of the package, e.g. "0.4.0" *)

(* ::Section:: *)
(*1. Ladder operators*)

(* ::Text:: *)
(*c[i, s] annihilates a fermion on site i with spin s, and cDag[i, s] creates it; the spin is SpinUp, SpinDown or a symbol. ** multiplies in the order written, and NormalOrder brings every cDag to the left.*)

NormalOrder[c[1, SpinUp] ** cDag[1, SpinUp]]
(* 1 - cDag[1, up] c[1, up]: the anticommutation relation *)

NormalOrder[cDag[1, SpinUp] ** cDag[1, SpinUp]]
(* 0: the Pauli principle *)

(* ::Text:: *)
(*Sites can be symbols: different symbols are different sites, so a result holds on any lattice.*)

NormalOrder[c[i, SpinUp] ** cDag[j, SpinUp]]
(* -cDag[j, up] c[i, up] *)

PlusHC[cDag[1, SpinUp] ** c[2, SpinUp]]
(* the hopping from 2 to 1 plus its Hermitian conjugate *)

(* ::Section:: *)
(*2. Local operators*)

(* ::Text:: *)
(*Products on one site are read as local operators, which stay compact: the densities nOp and nTotOp, the spin SzOp, SpOp, SmOp, the double occupancy DOp. ToSpinChargeBasis writes an expression with them.*)

ToSpinChargeBasis[cDag[1, SpinUp] ** c[1, SpinUp] ** cDag[1, SpinDown] ** c[1, SpinDown]]
(* DOp[1]: the double occupancy *)

ToSpinChargeBasis[SdotS[1, 1]]
(* 3/4 n - 3/2 D: a spin 1/2 lives only on a singly occupied site *)

(* ::Section:: *)
(*3. Spins on several sites*)

(* ::Text:: *)
(*Products on several sites are recognized as spin operators: the exchange, written with ladder operators, is S_i.S_j plus a density term.*)

ToSpinChargeBasis[1/2 Sum[cDag[i, s] ** c[i, r] ** cDag[j, r] ** c[j, s], {s, $Spins}, {r, $Spins}]]
(* S_i.S_j + n_i n_j/4 *)

(* ::Text:: *)
(*The algebra of spins comes for free: the commutator of two exchange terms that share a site is the scalar chirality.*)

ToSpinChargeBasis[OperatorCommutator[SdotS[1, 2], SdotS[2, 3]]]
(* -I S1.(S2 x S3) *)

(* ::Section:: *)
(*4. A model, as on paper*)

(* ::Text:: *)
(*The Hubbard model on a ring of L sites: hopping t between neighbours, plus its Hermitian conjugate, and the repulsion U on every site. The dimer, two sites joined once, is written apart.*)

hubbard[L_] := -t Sum[PlusHC[cDag[i, s] ** c[Mod[i, L] + 1, s]], {i, L}, {s, $Spins}] +
   U Sum[nOp[i, SpinUp] ** nOp[i, SpinDown], {i, L}];
dimer = -t Sum[PlusHC[cDag[1, s] ** c[2, s]], {s, $Spins}] + U Sum[nOp[i, SpinUp] ** nOp[i, SpinDown], {i, 2}];

HermitianOperatorQ[hubbard[4]]
(* True *)

(* ::Section:: *)
(*5. Solving it exactly*)

(* ::Text:: *)
(*FermiWick`Fock` writes the same expressions as sparse matrices: on the Fock space of the sites, or on a sector with a fixed number of particles N and spin Sz. The dimer with two electrons and Sz = 0, symbolically:*)

dimerSector = FockSector[FockSpace[{1, 2}], "N" -> 2, "Sz" -> 0];
Eigenvalues[Normal[FockMatrix[dimer, dimerSector]]]
(* 0, U and (U -+ Sqrt[U^2 + 16 t^2])/2 *)

(* ::Text:: *)
(*A ring of ten sites at half filling, with U = 4t: a sector of 63 504 states, solved in about a second.*)

ringSector = FockSector[FockSpace[Range[10]], "N" -> 10, "Sz" -> 0];
{energies, vectors} = LowestStates[FockMatrix[hubbard[10] /. {t -> 1, U -> 4}, ringSector], 1];
energies
(* {-5.834...}: the ground-state energy, in units of t *)

ground = First[vectors];
{FockExpectation[SdotS[1, 2], ground, ringSector], FockExpectation[DOp[1], ground, ringSector]}
(* about {-0.32, 0.10}: antiferromagnetic neighbours, and few doubly occupied sites *)

(* ::Section:: *)
(*6. Expectation values by Wick's theorem*)

(* ::Text:: *)
(*Approximate states are Gaussian: fixed by their contractions, one amplitude each. GaussianState[] has them all, and ContractionTable shows them between the modes of some sites: on one site, the densities of each spin and <S^+> off the diagonal.*)

MatrixForm[ContractionTable[GaussianState[], {1}, "Normal"]]
(* chi^{up up}_11, chi^{down down}_11 on the diagonal, chi^{up down}_11 = <S^+_1> and its conjugate off it *)

(* ::Text:: *)
(*A symmetry you assume removes contractions, or the indices that told them apart: with the charge conserved there is no pairing; with S^z conserved the hopping keeps the spin and the pairing joins opposite spins; with spin rotations (SU(2)) there is one hopping chi[i, j] and the singlet pairing Delta[i, j]; with time reversal the amplitudes of opposite spins are conjugates, and with SU(2) real.*)

MatrixForm[ContractionTable[GaussianState["SpinRotationZ"], {1, 2}, "Anomalous"]]
(* S^z conserved: the pairing joins opposite spins *)

(* ::Text:: *)
(*The state of most calculations has a short name: SpinSingletState[] is GaussianState[{"SpinRotation", "Charge"}, "Filling" -> 1/2], invariant under spin rotations, the charge conserved, at half filling.*)

MatrixForm[ContractionTable[SpinSingletState[], {1, 2}, "Normal"]]
(* 1/2 on each site, chi_12 between the sites, nothing between opposite spins *)

WickVEV[dimer, SpinSingletState[]]
(* -2 t (chi_12 + Conjugate[chi_12]) + U/2: the energy of the dimer in this state *)

(* ::Text:: *)
(*Products of any length: the correlation of two exchange terms, eight fermions, in a state with singlet pairing too, has sixty terms.*)

big = WickVEV[SdotS[1, 2] ** SdotS[3, 4], GaussianState["SpinRotation"]];
Length[big]
(* 60 *)

(* ::Text:: *)
(*The symmetries decide what survives: the chirality needs a state that breaks time reversal.*)

WickVEV[SdotScrossS[1, 2, 3], GaussianState[{"SpinRotation", "Charge"}]]
(* -3 Im(chi_12 chi_23 chi_31), written with chi_ji = Conjugate[chi_ij] *)

WickVEV[SdotScrossS[1, 2, 3], GaussianState[{"SpinRotation", "Charge", "TimeReversal"}]]
(* 0 *)

(* ::Section:: *)
(*7. States by rules, and modes with any labels*)

(* ::Text:: *)
(*A state can also be declared by rules for its contractions, with any spin structure: here a hopping g that keeps the spin and a hopping h that flips it.*)

flip = GaussianState[{{cDag[i_, s_], c[j_, r_]} :> If[i === j, KroneckerDelta[s, r]/2, If[s === r, g, h]]}];
WickVEV[SzOp[1] ** SzOp[2], flip]
(* (h^2 - g^2)/2 *)

(* ::Text:: *)
(*Modes can carry any labels; DeclareModeLabels says what each label is. With an orbital, the anticommutator gives a delta of the orbitals:*)

DeclareModeLabels[{"Site", "Orbital", "Spin"}];
NormalOrder[c[1, a, SpinUp] ** cDag[1, b, SpinUp]]
(* KroneckerDelta[a, b] - cDag[1, b, up] c[1, a, up] *)

(* ::Section:: *)
(*8. Mean field*)

(* ::Text:: *)
(*MeanField gives the Hartree-Fock-Bogoliubov Hamiltonian of an operator of any order: two operators kept out in every way, the rest contracted, and the constant that keeps the energy right. The Hubbard interaction in a state with magnetization m along z:*)

magnetic = GaussianState[{"SpinRotationZ", "Charge"}, "Filling" -> 1/2, "Magnetization" -> m];
ToSpinChargeBasis[MeanField[U nOp[1, SpinUp] ** nOp[1, SpinDown], magnetic]]
(* m^2 U - 2 m U Sz_1 + U n_1/2 - U/4 *)

(* ::Text:: *)
(*Six fermions are no harder. The fields are the coefficients of the result: the hopping field that the chirality produces on the bond i -> j, with symbolic sites:*)

fields = OperatorCoefficients[MeanField[-Jc SdotScrossS[i, j, k], SpinSingletState[]]];
fields[cDag[i, SpinUp] ** c[j, SpinUp]]
(* (3 I/4) Jc chi_jk Conjugate[chi_ik], that is (3 I/4) Jc chi_jk chi_ki *)

(* ::Section:: *)
(*9. Wick's theorem meets the exact solution*)

(* ::Text:: *)
(*The ground state of free fermions is Gaussian, so Wick's theorem must be exact there. A ring of six sites with three electrons of each spin: its contractions from the exact ground state, then the symbolic result of Wick's theorem with them, against the exact value.*)

freeSector = FockSector[FockSpace[Range[6]], "N" -> 6, "Sz" -> 0];
free = First[Last[LowestStates[FockMatrix[hubbard[6] /. {t -> 1, U -> 0}, freeSector], 1]]];
chi12 = FockExpectation[cDag[1, SpinUp] ** c[2, SpinUp], free, freeSector]
(* 0.333...: one third *)

{WickVEV[SdotS[1, 2], SpinSingletState[]] /. χ[1, 2] -> chi12, FockExpectation[SdotS[1, 2], free, freeSector]}
(* {-0.1667, -0.1667}: Wick's theorem, and the exact value *)

(* ::Section:: *)
(*10. Gutzwiller projectors*)

(* ::Text:: *)
(*FermiWick`Gutzwiller` handles projected operators: the hopping between two Gutzwiller-projected sites is reduced by a factor, here as a function of the variational angle theta.*)

real = GaussianState[{"SpinRotation", "Charge", "TimeReversal"}, "Filling" -> 1/2];
hop = Sum[PlusHC[cDag[1, s] ** c[2, s]], {s, $Spins}];
Simplify[GutzwillerVEV[LocalSandwich[GutzwillerProjector, {1, 2}, hop], real]/WickVEV[hop, real]]
(* Sin[2 theta]^2: the Gutzwiller reduction of the hopping, 8 d (1 - 2 d) with d = Sin[theta]^2/2 *)

(* ::Section:: *)
(*11. LaTeX*)

ToTeX[ToSpinChargeBasis[OperatorCommutator[SdotS[1, 2], SdotS[2, 3]]]]
(* "-i \mathbf S_{1}\cdot(\mathbf S_{2}\times\mathbf S_{3})" *)

WickTeX[cDag[1, SpinUp] ** cDag[2, SpinDown] ** c[3, SpinDown] ** c[4, SpinUp], SpinSingletState["Pairing" -> Δ]]
(* Wick's theorem written out, with the contractions drawn (LaTeX package simpler-wick) *)

(* ::Section:: *)
(*Where next*)

(* ::Text:: *)
(*The wiki (github.com/n4sd4q/FermiWick/wiki) has the recipes: writing operators, declaring a state, mean field, exact diagonalisation, Gutzwiller projectors, LaTeX. The showcase goes through the physics, and the manual describes every function; ?FermiWick`* lists the functions and ?Name explains each one.*)
