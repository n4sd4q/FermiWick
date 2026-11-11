(* ::Example:: load *)
Needs["FermiWick`"]; Needs["FermiWick`Gutzwiller`"];

(* ::Example:: symbolicspin *)
NormalOrder[c[i, \[Sigma]] ** cDag[j, \[Sigma]]]

(* ::Example:: samesite *)
NormalOrder[c[i, \[Sigma]] ** cDag[i, \[Sigma]]]

(* ::Example:: pauli *)
NormalOrder[cDag[1, "\[UpArrow]"] ** cDag[1, "\[UpArrow]"]]

(* ::Example:: labels *)
DeclareModeLabels[{"Site", "Orbital", "Spin"}];
NormalOrder[c[1, x, SpinUp] ** cDag[1, y, SpinUp] ** c[2, x, SpinDown]]

(* ::Example:: spins *)
NormalOrder[Sum[cDag[1, s] ** c[2, s], {s, $Spins}] ** c[3, SpinUp]]

(* ::Example:: composite *)
ToLadderBasis[SzOp[1]]

(* ::Example:: spincomponents *)
NormalOrder[Sum[SpinOperator[1, a] \[CenterDot] SpinOperator[2, a], {a, 3}] - SdotS[1, 2]]

(* ::Example:: projectors *)
ToSpinChargeBasis[EmptyProjector[1] + SingleOccupationProjector[1] + DOp[1]]

(* ::Example:: define *)
DefineOperator[Hop, 2, (cDag[#1, "\[UpArrow]"] ** c[#2, "\[UpArrow]"] + cDag[#1, "\[DownArrow]"] ** c[#2, "\[DownArrow]"]) &,
  Subscript["T", Row[{#1, #2}]] &];
ToSpinChargeBasis[NormalOrder[2 Hop[1, 2] ** nTotOp[3]]]

(* ::Example:: adjoint *)
DefineOperator[Pair, 2,
  (c[#1, "\[UpArrow]"] ** c[#2, "\[DownArrow]"] - c[#1, "\[DownArrow]"] ** c[#2, "\[UpArrow]"])/Sqrt[2] &,
  Subscript["b", Row[{#1, #2}]] &];
{Dagger[Hop[1, 2]], Dagger[Pair[1, 2]],
  ToSpinChargeBasis[Dagger[Pair[1, 2]] ** Pair[1, 2]]}

(* ::Example:: centerdot *)
ToSpinChargeBasis[SzOp[1] \[CenterDot] SzOp[1]]

(* ::Example:: anticommutator *)
OperatorAnticommutator[c[1, "\[UpArrow]"], cDag[1, "\[UpArrow]"]]

(* ::Example:: tospincharge *)
ToSpinChargeBasis[cDag[1, "\[UpArrow]"] ** c[1, "\[DownArrow]"] ** cDag[2, "\[DownArrow]"] ** c[2, "\[UpArrow]"]]

(* ::Example:: norecognize *)
ToSpinChargeBasis[SdotS[1, 2], "Recognize" -> False]

(* ::Example:: group *)
GroupLocalOperators[cDag[1, "\[UpArrow]"] ** cDag[2, "\[UpArrow]"] ** c[2, "\[UpArrow]"] ** c[1, "\[UpArrow]"]]

(* ::Example:: dagger *)
Dagger[t E^(I \[Alpha]) cDag[1, SpinUp] ** c[2, SpinUp] + \[Chi][1, 2] SpOp[1]]

(* ::Example:: plushc *)
PlusHC[g c[1, SpinUp] ** c[2, SpinDown], "Complex" -> {g}]

(* ::Example:: coefficients text *)
OperatorCoefficients[2 x cDag[1, "\[UpArrow]"] ** c[2, "\[UpArrow]"] + y + c[2, "\[UpArrow]"] ** cDag[1, "\[UpArrow]"]]

(* ::Example:: coefficient *)
OperatorCoefficients[2 x cDag[1, SpinUp] ** c[2, SpinUp], c[2, SpinUp] ** cDag[1, SpinUp]]

(* ::Example:: rulesstate *)
dimer = GaussianState[{{cDag[i_, s_], c[j_, t_]} :> KroneckerDelta[s, t] If[i === j, 1/2, g]}];
WickVEV[nOp[1, "\[UpArrow]"] ** nOp[2, "\[UpArrow]"], dimer]

(* ::Example:: channels *)
MatrixForm[ContractionTable[GaussianState[], {1}, "Normal"]]

(* ::Example:: symmetries *)
MatrixForm[ContractionTable[GaussianState[{"SpinRotation", "Charge"}, "Filling" -> \[Nu]], {1, 2}, "Normal"]]

(* ::Example:: contraction *)
{Contraction[dimer, cDag[1, "\[UpArrow]"], c[2, "\[UpArrow]"]], Contraction[dimer, c[1, "\[UpArrow]"], cDag[2, "\[UpArrow]"]]}

(* ::Example:: singletoptions *)
WickVEV[DOp[1], SpinSingletState["Filling" -> \[Nu], "Pairing" -> \[CapitalDelta], "OnSitePairing" -> \[Eta]]]

(* ::Example:: stateselect *)
WickVEV[nTotOp[1] ** nTotOp[2], StateSelect[SpinSingletState[], LocalPairQ]]

(* ::Example:: maxorder *)
WickVEV[SdotS[1, 2] ** SdotS[2, 3], "MaxOrder" -> 2]

(* ::Example:: verify *)
VerifyNormalOrder[c[1, "\[UpArrow]"] ** cDag[2, "\[UpArrow]"] ** c[2, "\[UpArrow]"] ** cDag[1, "\[UpArrow]"]]

(* ::Example:: meanfield *)
MeanField[nTotOp[1] ** nTotOp[2]]

(* ::Example:: meanfieldquadratic *)
MeanField[cDag[1, "\[UpArrow]"] ** c[2, "\[UpArrow]"]]

(* ::Example:: gutzwillerprojector *)
ToSpinChargeBasis[GutzwillerProjector[1]] /. \[Theta] -> 0

(* ::Example:: localsandwich *)
ToSpinChargeBasis[LocalSandwich[GutzwillerProjector, {1}, DOp[1]]]

(* ::Example:: gutzwillervev *)
GutzwillerVEV[LocalSandwich[GutzwillerProjector, {1}, nTotOp[1]]]

(* ::Example:: peierls *)
DefinePeierlsHopping[t, \[CurlyPhi]];
DefineFlux[\[CurlyPhi], {1, 2, 3}, \[CapitalPhi]];
Simplify[Subscript[t, 1, 2] Subscript[t, 2, 3] Subscript[t, 3, 1]]

(* ::Example:: totex text *)
ToTeX[SdotS[1, 2] - 2 DOp[1] + \[Chi][1, 2] cDag[1, SpinUp] ** c[2, SpinUp]]

(* ::Example:: wickpairings text *)
WickPairings[cDag[1, SpinUp] ** c[2, SpinUp] ** cDag[2, SpinUp] ** c[1, SpinUp]]

(* ::Example:: wicktex tex *)
WickTeX[cDag[1, SpinUp] ** c[2, SpinUp] ** cDag[2, SpinUp] ** c[1, SpinUp], "Values" -> False]

(* ::Example:: fockload *)
Needs["FermiWick`Fock`"]

(* ::Example:: fockspace text *)
FockSpace[{1, 2}, "Orbitals" -> {"s", "p"}]["Modes"]

(* ::Example:: modeindex text *)
{ModeIndex[FockSpace[{1, 2}], c[2, "\[UpArrow]"]], FockSpace[Range[3], "Spins" -> None]["Modes"]}

(* ::Example:: focksector text *)
{FockSector[FockSpace[Range[6]], "N" -> 6, "Sz" -> 0]["Dimension"],
 FockSector[FockSpace[{1, 2}], "Counts" -> {{{c[1, "\[UpArrow]"], c[1, "\[DownArrow]"]}, 2}}]["Dimension"]}

(* ::Example:: fockmatrix *)
Normal[FockMatrix[cDag[1, "\[UpArrow]"] ** c[2, "\[UpArrow]"], FockSector[FockSpace[{1, 2}], "N" -> 1]]]

(* ::Example:: fockmatrixsymbolic *)
Normal[FockMatrix[J SdotS[1, 2], FockSector[FockSpace[{1, 2}], "N" -> 2, "Sz" -> 0]]]

(* ::Example:: fockexpectation *)
sec = FockSector[FockSpace[{1, 2}], "N" -> 2, "Sz" -> 0];
singlet = {0, 1, -1, 0}/Sqrt[2];
{FockExpectation[SdotS[1, 2], singlet, sec], FockExpectation[DOp[1], singlet, sec]}

(* ::Example:: statetoket text *)
StateToKet[sec, singlet]

(* ::Example:: basisstate text *)
{OccupationBasis[sec], Normal[BasisState[sec, {1, 0, 0, 1}]]}

(* ::Example:: spinsquared text *)
Eigenvalues[Normal[SpinSquared[sec]]]

(* ::Example:: lowest text *)
dimer = -Sum[cDag[1, s] ** c[2, s] + cDag[2, s] ** c[1, s], {s, $Spins}] + 4 (DOp[1] + DOp[2]);
N[Round[10^4 First[LowestStates[FockMatrix[dimer, sec], 2]]]/10^4]

(* ::Example:: lehmann text *)
space = FockSpace[{1, 2}];
g = LehmannGreen[FockMatrix[dimer, space], space, {c[1, "\[UpArrow]"]}, "N" -> 2];
{N[Round[10^4 g["GroundEnergy"]]/10^4], g["Degeneracy"], Max[Abs[Total[g["Residues"]] - IdentityMatrix[1]]] < 10^-12}

(* ::Example:: greenmatrix text *)
N[Round[10^4 GreenMatrix[g, 0.5 I]]/10^4]

(* ::Example:: particlehole text *)
With[{ph = ParticleHoleOperator[space]}, Max[Abs[ph . FockMatrix[c[1, "\[UpArrow]"], space] . Transpose[ph] - FockMatrix[cDag[1, "\[UpArrow]"], space]]]]

(* ::Example:: sectorblock text *)
SectorBlock[FockMatrix[DOp[1] + DOp[2], space], sec] == FockMatrix[DOp[1] + DOp[2], sec]
