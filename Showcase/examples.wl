(* ::Example:: car *)
NormalOrder[c[1, "\[UpArrow]"] ** cDag[1, "\[UpArrow]"]]

(* ::Example:: distinctsites *)
NormalOrder[c[i, \[Sigma]] ** cDag[j, \[Sigma]]]

(* ::Example:: samesite *)
NormalOrder[c[i, \[Sigma]] ** cDag[i, \[Sigma]]]

(* ::Example:: pauli *)
NormalOrder[cDag[1, "\[UpArrow]"] ** cDag[2, "\[DownArrow]"] ** cDag[1, "\[UpArrow]"]]

(* ::Example:: commutator1 *)
OperatorCommutator[SzOp[1], cDag[1, "\[UpArrow]"]]

(* ::Example:: anticommutator *)
OperatorAnticommutator[c[1, "\[UpArrow]"], cDag[1, "\[UpArrow]"] + cDag[2, "\[UpArrow]"]]

(* ::Example:: ladder *)
ToLadderBasis[SpOp[1] + SzOp[1]]

(* ::Example:: compact *)
ToSpinChargeBasis[DOp[1] ** nTotOp[1] ** SzOp[2]]

(* ::Example:: spinflip *)
ToSpinChargeBasis[SpOp[1] \[CenterDot] SmOp[1]]

(* ::Example:: szsz *)
ToSpinChargeBasis[SzOp[1] \[CenterDot] SzOp[1]]

(* ::Example:: spinlength *)
ToSpinChargeBasis[SdotS[1, 1]]

(* ::Example:: projectors *)
ToSpinChargeBasis[EmptyProjector[1] + SingleOccupationProjector[1] + DOp[1]]

(* ::Example:: projectorproduct *)
ToSpinChargeBasis[SingleOccupationProjector[1] \[CenterDot] SpOp[1] \[CenterDot] SingleOccupationProjector[1] - SpOp[1]]

(* ::Example:: group *)
GroupLocalOperators[cDag[1, "\[UpArrow]"] ** cDag[2, "\[DownArrow]"] ** c[2, "\[DownArrow]"] ** c[1, "\[UpArrow]"]]

(* ::Example:: decompose *)
ToSpinChargeBasis[cDag[1, "\[UpArrow]"] ** c[1, "\[DownArrow]"] ** cDag[2, "\[DownArrow]"] ** c[2, "\[UpArrow]"]]

(* ::Example:: recognizess *)
ToSpinChargeBasis[Sum[SpinOperator[1, a] \[CenterDot] SpinOperator[2, a], {a, 3}]]

(* ::Example:: norecognize *)
ToSpinChargeBasis[SdotS[1, 2], "Recognize" -> False]

(* ::Example:: xxz *)
ToSpinChargeBasis[J (SpinOperator[1, 1] \[CenterDot] SpinOperator[2, 1] + SpinOperator[1, 2] \[CenterDot] SpinOperator[2, 2]) +
  Jz SzOp[1] \[CenterDot] SzOp[2]]

(* ::Example:: commutator *)
ToSpinChargeBasis[OperatorCommutator[SdotS[1, 2], SdotS[2, 3]]]

(* ::Example:: chiralityorder *)
ToSpinChargeBasis[SdotScrossS[2, 1, 3]]

(* ::Example:: invariance *)
ToSpinChargeBasis[OperatorCommutator[SdotS[1, 2] + SdotS[2, 3] + SdotS[3, 1], SpOp[1] + SpOp[2] + SpOp[3]]]

(* ::Example:: define *)
DefineOperator[Hop, 2, (cDag[#1, "\[UpArrow]"] ** c[#2, "\[UpArrow]"] + cDag[#1, "\[DownArrow]"] ** c[#2, "\[DownArrow]"]) &,
  Subscript["T", Row[{#1, #2}]] &];
ToSpinChargeBasis[OperatorCommutator[Hop[1, 2], nTotOp[2]]]

(* ::Example:: definedagger *)
ToSpinChargeBasis[Dagger[Hop[1, 2]]]

(* ::Example:: dagger *)
Dagger[2 I J cDag[1, SpinUp] ** c[2, SpinDown] + \[Chi][1, 2] SpOp[1] ** Bdag[2]]

(* ::Example:: plushc *)
PlusHC[t E^(I \[Alpha]) cDag[1, SpinUp] ** c[2, SpinUp]]

(* ::Example:: hermitianq *)
{HermitianOperatorQ[SdotScrossS[1, 2, 3]], HermitianOperatorQ[I SdotScrossS[1, 2, 3]]}

(* ::Example:: spinless *)
NormalOrder[c[1] ** cDag[2] ** c[2] ** cDag[1]]

(* ::Example:: orbitals *)
DeclareModeLabels[{"Site", "Orbital", "Spin"}];
NormalOrder[c[1, a, SpinUp] ** cDag[1, b, SpinUp]]

(* ::Example:: momenta *)
DeclareModeLabels[{"Momentum"}];
WickVEV[cDag[k] ** c[k] ** cDag[q] ** c[q], GaussianState[{{cDag[p1_], c[p2_]} :> KroneckerDelta[p1, p2] nk[p1]}]]

(* ::Example:: coefficients *)
OperatorCoefficients[2 x cDag[1, "\[UpArrow]"] ** c[2, "\[UpArrow]"] + y + c[2, "\[UpArrow]"] ** cDag[1, "\[UpArrow]"]]

(* ::Example:: totex text *)
ToTeX[ToSpinChargeBasis[OperatorCommutator[SdotS[1, 2], SdotS[2, 3]]]]

(* ::Example:: peierls *)
DefinePeierlsHopping[t, \[CurlyPhi]];
DefineFlux[\[CurlyPhi], {1, 2, 3}, \[CapitalPhi]];
Simplify[Subscript[t, 1, 2] Subscript[t, 2, 3] Subscript[t, 3, 1]]

(* ::Example:: generallocal *)
MatrixForm[ContractionTable[GaussianState[], {1}, "Normal"]]

(* ::Example:: generalpairing *)
MatrixForm[ContractionTable[GaussianState[], {1, 2}, "Anomalous"]]

(* ::Example:: symz *)
MatrixForm[ContractionTable[GaussianState["SpinRotationZ"], {1, 2}, "Anomalous"]]

(* ::Example:: symsu2tr *)
MatrixForm[ContractionTable[GaussianState[{"SpinRotation", "TimeReversal"}], {1, 2}, "Normal"]]

(* ::Example:: szexchange *)
WickVEV[SdotS[1, 2], GaussianState[{"SpinRotationZ", "Charge"}]]

(* ::Example:: rulesstate *)
dimer = GaussianState[{{cDag[i_, s_], c[j_, t_]} :> KroneckerDelta[s, t] If[i === j, 1/2, g]}];
WickVEV[nOp[1, "\[UpArrow]"] ** nOp[2, "\[UpArrow]"], dimer]

(* ::Example:: contraction *)
{Contraction[dimer, cDag[1, "\[UpArrow]"], c[2, "\[UpArrow]"]], Contraction[dimer, c[1, "\[UpArrow]"], cDag[2, "\[UpArrow]"]]}

(* ::Example:: flipstate *)
flip = GaussianState[{{cDag[i_, s_], c[j_, t_]} :> If[i === j, KroneckerDelta[s, t]/2, If[s === t, g, h]]}];
WickVEV[SzOp[1] ** SzOp[2], flip]

(* ::Example:: singletequal *)
Simplify[WickVEV[SdotScrossS[1, 2, 3], GaussianState[{"SpinRotation", "Charge"}, "Filling" -> 1/2]] -
  WickVEV[SdotScrossS[1, 2, 3], SpinSingletState[]]]

(* ::Example:: singletstate text *)
state = SpinSingletState["Pairing" -> \[CapitalDelta]]

(* ::Example:: exchange *)
WickVEV[SdotS[1, 2], state]

(* ::Example:: chirality *)
WickVEV[SdotScrossS[1, 2, 3], state]

(* ::Example:: doubleocc *)
WickVEV[DOp[1], SpinSingletState["Filling" -> \[Nu], "Pairing" -> \[CapitalDelta], "OnSitePairing" -> \[Eta]]]

(* ::Example:: magnetic text *)
sgn["\[UpArrow]"] = 1; sgn["\[DownArrow]"] = -1;
neel = GaussianState[{{cDag[i_, s_], c[i_, t_]} :> KroneckerDelta[s, t] (1/2 + sgn[s] m[i])}]

(* ::Example:: magneticexchange *)
WickVEV[SdotS[1, 2], neel]

(* ::Example:: magneticsym *)
WickVEV[SdotS[1, 2], StateSelect[GaussianState["SpinRotationZ", "Magnetization" -> (m[#] &)], LocalPairQ]]

(* ::Example:: general *)
WickVEV[SpOp[1] + SzOp[1], GaussianState[]]

(* ::Example:: generalanomalous *)
WickVEV[Bdag[1] ** B[2], GaussianState[]]

(* ::Example:: select *)
WickVEV[SdotS[1, 2], StateSelect[state, ! AnomalousPairQ[##] &]]

(* ::Example:: localonly *)
WickVEV[nTotOp[1] ** nTotOp[2], StateSelect[SpinSingletState[], LocalPairQ]]

(* ::Example:: nonlocalonly *)
WickVEV[nTotOp[1] ** nTotOp[2], StateSelect[SpinSingletState[], ! LocalPairQ[##] &]]

(* ::Example:: maxorder *)
WickVEV[SdotS[1, 2] ** SdotS[2, 3], "MaxOrder" -> 2]

(* ::Example:: verify *)
VerifyNormalOrder[SdotS[1, 2] ** SdotS[2, 3], state]

(* ::Example:: pairings text *)
WickPairings[cDag[1, SpinUp] ** c[2, SpinUp] ** cDag[2, SpinUp] ** c[1, SpinUp]]

(* ::Example:: wicktex tex *)
WickTeX[cDag[1, SpinUp] ** cDag[2, SpinDown] ** c[3, SpinDown] ** c[4, SpinUp], state]

(* ::Example:: wicktexloop tex *)
WickTeX[cDag[1, SpinUp] ** c[2, SpinUp] ** cDag[2, SpinUp] ** c[3, SpinUp] ** cDag[3, SpinUp] ** c[1, SpinUp],
  StateSelect[state, ! LocalPairQ[##] &]]

(* ::Example:: hubbard *)
MeanField[U DOp[1], SpinSingletState["Filling" -> \[Nu], "Pairing" -> \[CapitalDelta], "OnSitePairing" -> \[Eta]]]

(* ::Example:: magnetichubbard *)
MeanField[U DOp[1], neel]

(* ::Example:: density *)
MeanField[nTotOp[1] ** nTotOp[2]]

(* ::Example:: exchangefields *)
mf = MeanField[J SdotS[1, 2], state];
fields = OperatorCoefficients[mf];
{fields[cDag[1, "\[UpArrow]"] ** c[2, "\[UpArrow]"]], fields[cDag[1, "\[UpArrow]"] ** cDag[2, "\[DownArrow]"]], fields[1]}

(* ::Example:: chiralfields *)
mfc = OperatorCoefficients[MeanField[-Jc SdotScrossS[1, 2, 3], state]];
{mfc[cDag[1, "\[UpArrow]"] ** c[2, "\[UpArrow]"]], mfc[cDag[1, "\[UpArrow]"] ** cDag[2, "\[DownArrow]"]]}

(* ::Example:: hermitianmf *)
HermitianOperatorQ[MeanField[-Jc SdotScrossS[1, 2, 3], state]]

(* ::Example:: variational *)
Simplify[WickVEV[MeanField[SdotScrossS[1, 2, 3], state], state] - WickVEV[SdotScrossS[1, 2, 3], state]]

(* ::Example:: quadratic *)
MeanField[cDag[1, "\[UpArrow]"] ** c[2, "\[UpArrow]"] + cDag[2, "\[UpArrow]"] ** c[1, "\[UpArrow]"]]

(* ::Example:: gutzwillerprojector *)
ToSpinChargeBasis[GutzwillerProjector[1]]

(* ::Example:: gutzwillernorm *)
GutzwillerVEV[GutzwillerProjector[1] \[CenterDot] GutzwillerProjector[1]]

(* ::Example:: gutzwillerd *)
GutzwillerVEV[LocalSandwich[GutzwillerProjector, {1}, DOp[1]]]

(* ::Example:: hopping *)
real = GaussianState[{"SpinRotation", "Charge", "TimeReversal"}, "Filling" -> 1/2];
hop = Sum[cDag[1, s] ** c[2, s] + cDag[2, s] ** c[1, s], {s, $Spins}];
Simplify[GutzwillerVEV[LocalSandwich[GutzwillerProjector, {1, 2}, hop], real]/WickVEV[hop, real]]

(* ::Example:: brq *)
Simplify[Sin[2 \[Theta]]^2 - 8 d (1 - 2 d) /. d -> Sin[\[Theta]]^2/2]

(* ::Example:: brd *)
d /. First[Solve[D[eps0 8 d (1 - 2 d) + U d, d] == 0, d]]

(* ::Example:: projectedexchange *)
Simplify[WickVEV[LocalSandwich[GutzwillerProjector, {1, 2}, SdotS[1, 2]]]]

(* ::Example:: fockload *)
Needs["FermiWick`Fock`"]

(* ::Example:: fockspace text *)
space = FockSpace[{1, 2}];
space["Modes"]

(* ::Example:: focksector text *)
sector = FockSector[space, "N" -> 2, "Sz" -> 0];
OccupationBasis[sector]

(* ::Example:: dimerspectrum *)
H = -t Sum[PlusHC[cDag[1, s] ** c[2, s]], {s, $Spins}] + U/2 Sum[(nTotOp[i] - 1)^2, {i, 2}];
Eigenvalues[Normal[FockMatrix[H, sector]]]

(* ::Example:: dimerground *)
{vals, vecs} = Eigensystem[Normal[FockMatrix[H, sector]]];
ground = Normalize[vecs[[First[Ordering[N[vals /. {t -> 1, U -> 1}], 1]]]]];
FullSimplify[StateToKet[sector, ground], t > 0 && U > 0]

(* ::Example:: dimerpartition *)
Z = PartitionFunction[H, space, T, \[Mu]]

(* ::Example:: dimerfilling *)
\[CapitalOmega] = -T Log[Z];
FullSimplify[{-D[\[CapitalOmega], \[Mu]], -D[\[CapitalOmega], {\[Mu], 2}]} /. \[Mu] -> 0, T > 0]

(* ::Example:: dimerentropy *)
S = -D[\[CapitalOmega], T];
{Limit[S /. {\[Mu] -> 0, U -> 4 t}, T -> 0, Assumptions -> t > 0], Limit[S /. \[Mu] -> 0, T -> Infinity]}

(* ::Example:: dimergreen *)
g = LehmannGreen[H, space, {c[1, "\[UpArrow]"], c[2, "\[UpArrow]"]}, "ChemicalPotential" -> \[Mu],
  "Assumptions" -> t > 0 && U > 0 && Abs[\[Mu]] < Sqrt[16 t^2 + U^2]/2 - t];
First[g["GroundSectors"]]

(* ::Example:: dimergreenmatrix *)
P = {{1, 1}, {1, -1}}/Sqrt[2];
MatrixForm[Simplify[P . GreenMatrix[g, z] . P, t > 0 && U > 0]]

(* ::Example:: dimergreenT *)
gT = LehmannGreen[H, space, {c[1, "\[UpArrow]"], c[2, "\[UpArrow]"]}, "Temperature" -> T, "ChemicalPotential" -> \[Mu]];
bonding = Select[Transpose[{gT["Poles"], Z (P . # . P)[[1, 1]] & /@ gT["Residues"]}], ! PossibleZeroQ[Last[#]] &];
MatrixForm[FullSimplify[bonding, T > 0 && t > 0 && U > 0]]

(* ::Example:: ring8 text *)
nb[i_] := Mod[i - 1, 8] + 1;
hop8 = -Sum[cDag[i, s] ** c[nb[i + 1], s] + cDag[nb[i + 1], s] ** c[i, s], {i, 8}, {s, $Spins}];
half = FockSector[FockSpace[Range[8]], "N" -> 8, "Sz" -> 0];
{e, v} = LowestStates[FockMatrix[hop8 + 4 Sum[DOp[i], {i, 8}], half], 1];
N[Round[10^4 {First[e]/8, FockExpectation[DOp[1], First[v], half],
    FockExpectation[SdotS[1, 2], First[v], half]}]/10^4]

(* ::Example:: lehmann text *)
ring4 = -Sum[PlusHC[cDag[i, s] ** c[Mod[i, 4] + 1, s]], {i, 4}, {s, $Spins}] + 2 Sum[(nTotOp[i] - 1)^2, {i, 4}];
g4 = LehmannGreen[ring4, FockSpace[Range[4]], {c[1, "\[UpArrow]"]}, "Temperature" -> 1/2, "ChemicalPotential" -> 0];
{N[Round[10^4 g4["GrandPotential"]]/10^4], Length[g4["Poles"]], N[Round[10^4 Im[GreenMatrix[g4, I Pi/2][[1, 1]]]]/10^4]}

(* ::Example:: wickexact text *)
six = FockSector[FockSpace[Range[8]], "N" -> 6, "Sz" -> 0];
v0 = First[Last[LowestStates[FockMatrix[hop8, six], 1]]];
pairs = Flatten[Table[{cDag[i, s], c[j, s]}, {i, 8}, {j, 8}, {s, $Spins}], 2];
exact = GaussianState[Thread[pairs -> (FockExpectation[#1 ** #2, v0, six] & @@@ pairs)]];
op = SdotS[1, 2] + DOp[1];
Chop[WickVEV[op, exact] - FockExpectation[op, v0, six]]
