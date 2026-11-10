(* FermiWick test suite. Run with  math -script Tests/RunTests.wls  (or TestReport["Tests/FermiWick.wlt"]). *)

testDir = DirectoryName[$TestFileName]; root = ParentDirectory[testDir];
PacletDirectoryLoad[FileNameJoin[{root, "FermiWick"}]];
Needs["FermiWick`"]; Needs["FermiWick`Gutzwiller`"]; Needs["FermiWick`Fock`"];
VerificationTest[MemberQ[$Packages, "FermiWick`"] && MemberQ[$Packages, "FermiWick`Gutzwiller`"], True, TestID -> "load"]

up = "\[UpArrow]"; dn = "\[DownArrow]";
zeroQ[x_] := AllTrue[Values[OperatorCoefficients[x]], Simplify[#] === 0 &];
conj[e_] := e /. {Complex[a_, b_] :> Complex[a, -b], Conjugate[x_] :> x, (h : \[Chi] | \[CapitalDelta])[a__] :> Conjugate[h[a]]};
hop[i_, j_] := cDag[i, up] ** c[j, up] + cDag[i, dn] ** c[j, dn];
pairDag[i_, j_] := cDag[i, up] ** cDag[j, dn] - cDag[i, dn] ** cDag[j, up];
pair[i_, j_] := c[j, dn] ** c[i, up] - c[j, up] ** c[i, dn];
bond[i_, j_, t_, d_] := t hop[i, j] + conj[t] hop[j, i] + d pairDag[i, j] + conj[d] pair[i, j];
x[a_, b_] := If[OrderedQ[{a, b}], \[Chi][a, b], Conjugate[\[Chi][b, a]]]; d[a_, b_] := \[CapitalDelta] @@ Sort[{a, b}]; ds[a_, b_] := Conjugate[d[a, b]];
sing = SpinSingletState[]; pairing = SpinSingletState["Pairing" -> \[CapitalDelta]];

(* ---------------------------------------------------------------- algebra *)
VerificationTest[NormalOrder[c[1, up] ** cDag[1, up]], 1 - cDag[1, up] ** c[1, up], TestID -> "CAR"]
VerificationTest[NormalOrder[c[1, up] ** c[1, up]], 0, TestID -> "Pauli"]
VerificationTest[NormalOrder[(2 c[1, up]) ** (a cDag[2, dn] + b c[3, up])], -2 a cDag[2, dn] ** c[1, up] + 2 b c[1, up] ** c[3, up],
  TestID -> "linearity of **"]
VerificationTest[ToSpinChargeBasis[SdotScrossS[j, i, k]], -SdotScrossS[i, j, k], TestID -> "recognition of the chirality"]
VerificationTest[NormalOrder[SpinOperator[i, 1]\[CenterDot]SpinOperator[j, 1] + SpinOperator[i, 2]\[CenterDot]SpinOperator[j, 2] +
   SpinOperator[i, 3]\[CenterDot]SpinOperator[j, 3] - SdotS[i, j]], 0, TestID -> "S.S from the components"]
VerificationTest[NormalOrder[OperatorCommutator[SzOp[1], c[1, up]]], -c[1, up]/2, TestID -> "commutator"]

(* ---------------------------------------------------------------- mean field of the spin operators and of the Hubbard interaction *)
VerificationTest[Simplify[WickVEV[SdotS[1, 2], pairing] + 3/2 (x[1, 2] x[2, 1] + d[1, 2] ds[1, 2])], 0, TestID -> "<S.S> with pairing"]
z = x[1, 2] x[2, 3] x[3, 1] + x[1, 2] ds[2, 3] d[3, 1] + x[2, 3] ds[3, 1] d[1, 2] + x[3, 1] ds[1, 2] d[2, 3];
chir = WickVEV[SdotScrossS[1, 2, 3], pairing];
VerificationTest[Simplify[chir - 3/(2 I) (z - conj[z])], 0, TestID -> "<S.(SxS)> with pairing"]
VerificationTest[zeroQ[MeanField[SdotS[1, 2], pairing] -
   (bond[1, 2, -3/4 x[2, 1], 3/4 d[1, 2]] + 3/2 (x[1, 2] x[2, 1] + d[1, 2] ds[1, 2]))], True, TestID -> "MeanField[S.S]"]
dt[a_, b_, k_] := 3 I/4 (x[b, k] x[k, a] + ds[b, k] d[k, a]);
dD[a_, b_, k_] := 3 I/4 (x[k, b] d[k, a] - x[k, a] d[b, k]);
VerificationTest[zeroQ[MeanField[-SdotScrossS[1, 2, 3], pairing] -
   (Total[bond[#1, #2, dt[#1, #2, #3], dD[#1, #2, #3]] & @@@ {{1, 2, 3}, {2, 3, 1}, {3, 1, 2}}] + 2 chir)], True,
  TestID -> "MeanField[-S.(SxS)]"]
VerificationTest[zeroQ[MeanField[DOp[1], SpinSingletState["Filling" -> \[Nu], "Pairing" -> \[CapitalDelta], "OnSitePairing" -> \[Eta]]] -
   (\[Nu] nTotOp[1] - \[Eta] cDag[1, up] ** cDag[1, dn] - Conjugate[\[Eta]] c[1, dn] ** c[1, up] - \[Nu]^2 - \[Eta] Conjugate[\[Eta]])], True,
  TestID -> "MeanField of Hubbard U: Hartree + BCS"]
VerificationTest[With[{st = SpinSingletState["Filling" -> \[Nu], "Pairing" -> \[CapitalDelta], "OnSitePairing" -> \[Eta]]},
   AllTrue[{SdotS[1, 2], SdotScrossS[1, 2, 3], DOp[1], nTotOp[1] ** nTotOp[2], SdotS[1, 2] ** SzOp[3]},
    Simplify[WickVEV[MeanField[#, st], st] - WickVEV[#, st]] === 0 &]], True, TestID -> "<MeanField[X]> = <X>"]

(* ---------------------------------------------------------------- states *)
rulesState = GaussianState[{{cDag[i_, s_], c[j_, t_]} :> KroneckerDelta[s, t] If[i === j, 1/2, x[i, j]]}];
VerificationTest[Simplify[WickVEV[SdotScrossS[1, 2, 3], rulesState] - WickVEV[SdotScrossS[1, 2, 3], sing]], 0,
  TestID -> "state from rules = SpinSingletState (<c cDag> from the anticommutator)"]
VerificationTest[Simplify[WickVEV[SdotS[1, 2], StateSelect[pairing, ! AnomalousPairQ[##] &]] - WickVEV[SdotS[1, 2], sing]], 0,
  TestID -> "StateSelect without anomalous contractions"]
VerificationTest[WickVEV[SdotS[1, 2], StateSelect[sing, LocalPairQ]], 0, TestID -> "StateSelect local only: S.S has no local part"]
VerificationTest[With[{g = WickVEV[SdotS[1, 2], GeneralState[]] /. \[Rho][a_, s_, b_, t_] :> KroneckerDelta[s, t] If[a === b, 1/2, x[a, b]]},
   Simplify[g + 3/2 x[1, 2] x[2, 1]]], 0, TestID -> "GeneralState reduces to the singlet one"]
VerificationTest[GutzwillerVEV[GutzwillerProjector[i]\[CenterDot]GutzwillerProjector[i]], 1, TestID -> "<P_G^2> = 1 at half filling"]

(* ---------------------------------------------------------------- numbers: exact Fock space of 3 sites (FockMatrix) *)
modes = Tuples[{{1, 2, 3}, {up, dn}}]; nm = Length[modes]; sp3 = FockSpace[{1, 2, 3}];
jw = Normal[FockMatrix[c @@ #, sp3]] & /@ modes;              (* the matrices of c on the modes, in this order (FermiWick`Fock`) *)
SeedRandom[7];
hA = RandomComplex[{-1 - I, 1 + I}, {nm, nm}]; hA = hA + ConjugateTranspose[hA];
hB = RandomComplex[{-1 - I, 1 + I}, {nm, nm}]; hB = hB - Transpose[hB];
hFock = Sum[hA[[a, b]] ConjugateTranspose[jw[[a]]] . jw[[b]] + hB[[a, b]] ConjugateTranspose[jw[[a]]] . ConjugateTranspose[jw[[b]]] +
     Conjugate[hB[[a, b]]] jw[[b]] . jw[[a]], {a, nm}, {b, nm}];
omega = With[{e = MatrixExp[-1.3 hFock]}, e/Tr[e]];
mat[op_] := Normal[FockMatrix[op, sp3]];
toMatrix[expr_] := Normal[FockMatrix[expr, sp3]];
numbers = {\[Rho][a_, s_, b_, t_] :> Tr[omega . mat[cDag[a, s]] . mat[c[b, t]]], \[Kappa][a_, s_, b_, t_] :> Tr[omega . mat[c[a, s]] . mat[c[b, t]]]};
gen = GeneralState["Anomalous" -> \[Kappa]];
tests = {SdotS[1, 2], SdotScrossS[1, 2, 3], DOp[1] ** DOp[2], nTotOp[1] ** SzOp[2] ** SpOp[3], Bdag[1] ** B[3] ** nTotOp[2]};
VerificationTest[Max[Abs[(WickVEV[#, gen] /. numbers) - Tr[omega . toMatrix[#]] & /@ tests]] < 10^-10, True,
  TestID -> "WickVEV = exact Fock-space value (general state with pairing)"]
VerificationTest[Max[Abs[(WickVEV[MeanField[#, gen], gen] /. numbers) - Tr[omega . toMatrix[#]] & /@ tests]] < 10^-10, True,
  TestID -> "<MeanField[X]> = exact <X>"]

(* ---------------------------------------------------------------- conveniences *)
VerificationTest[{c[1, SpinUp], $Spins}, {c[1, "\[UpArrow]"], {"\[UpArrow]", "\[DownArrow]"}}, TestID -> "SpinUp, SpinDown, $Spins"]
VerificationTest[With[{h = Sum[cDag[1, s] ** c[2, s], {s, $Spins}]},
   {OperatorCoefficients[3 h, cDag[1, SpinUp] ** c[2, SpinUp]], OperatorCoefficients[3 h, c[2, SpinUp] ** cDag[1, SpinUp]]}], {3, -3},
  TestID -> "OperatorCoefficients of a monomial in any order"]

(* ---------------------------------------------------------------- LaTeX and Wick written out *)
VerificationTest[With[{x = c[1, SpinUp] ** cDag[2, SpinUp] ** c[2, SpinUp] ** cDag[1, SpinUp] ** cDag[3, SpinDown] ** c[3, SpinDown]},
   Simplify[Total[#["Sign"] #["Value"] & /@ WickPairings[x, pairing]] - WickVEV[x, pairing, "NormalOrder" -> False]]], 0,
  TestID -> "WickPairings sum to WickVEV"]
VerificationTest[StringContainsQ[WickTeX[cDag[1, SpinUp] ** c[2, SpinUp] ** cDag[2, SpinUp] ** c[1, SpinUp]], "\\wick{\\c2 c^\\dagger_{1\\uparrow} \\c1 c_{2\\uparrow}"],
  True,
  TestID -> "WickTeX"]
VerificationTest[ToTeX[SdotS[1, 2] + cDag[1, SpinUp] ** c[2, SpinDown] + SzOp[3]], "c^\\dagger_{1\\uparrow}c_{2\\downarrow}+\\mathbf S_{1}\\cdot\\mathbf S_{2}+S^z_{3}",
  TestID -> "ToTeX"]

(* ---------------------------------------------------------------- states from symmetries *)
VerificationTest[Simplify[WickVEV[SdotScrossS[1, 2, 3], GaussianState[{"SpinRotation", "Charge"}, "Filling" -> 1/2]] -
   WickVEV[SdotScrossS[1, 2, 3], sing]], 0, TestID -> "GaussianState[{SpinRotation, Charge}] = SpinSingletState"]
VerificationTest[ContractionTable[GaussianState[{"SpinRotation", "TimeReversal"}], {1, 2}, "Normal"],
  {{\[Chi][1, 1], 0, \[Chi][1, 2], 0}, {0, \[Chi][1, 1], 0, \[Chi][1, 2]}, {\[Chi][1, 2], 0, \[Chi][2, 2], 0}, {0, \[Chi][1, 2], 0, \[Chi][2, 2]}},
  TestID -> "SU(2) and time reversal: real symmetric hopping"]
VerificationTest[ContractionTable[GaussianState["SpinRotationZ"], {1, 2}, "Anomalous"][[1, 4]], \[CapitalDelta][1, 2] + \[CapitalDelta][1, 2, "z"],
  TestID -> "rotations about z: singlet and m = 0 triplet"]
VerificationTest[GaussianState["Foo"], $Failed, {GaussianState::sym}, TestID -> "unknown symmetry"]

(* exact channels of Gibbs states of quadratic Hamiltonians with a given symmetry *)
pauli = {{{0, 1}, {1, 0}}, {{0, -I}, {I, 0}}, {{1, 0}, {0, -1}}}; eps = {{0, 1}, {-1, 0}};
gibbs[hA_, hB_] := With[{h = Sum[hA[[a, b]] ConjugateTranspose[jw[[a]]] . jw[[b]] + hB[[a, b]] ConjugateTranspose[jw[[a]]] . ConjugateTranspose[jw[[b]]] +
       Conjugate[hB[[a, b]]] jw[[b]] . jw[[a]], {a, nm}, {b, nm}]}, With[{e = MatrixExp[-1.1 h]}, e/Tr[e]]];
blockG[w_, i_, j_] := Table[Tr[w . mat[cDag[i, s]] . mat[c[j, t]]], {s, {up, dn}}, {t, {up, dn}}];
blockF[w_, i_, j_] := Table[Tr[w . mat[c[i, s]] . mat[c[j, t]]], {s, {up, dn}}, {t, {up, dn}}];
channels[w_] := Flatten[Table[{
     \[Chi][i, j] -> Tr[blockG[w, i, j]]/2, Table[\[Chi][i, j, a] -> Tr[blockG[w, i, j] . Transpose[pauli[[k]]]]/2, {a, {"x", "y", "z"}}, {k, {Position[{"x", "y", "z"}, a][[1, 1]]}}],
     \[CapitalDelta][i, j] -> Tr[blockF[w, i, j] . (-eps)]/2, Table[\[CapitalDelta][i, j, a] -> Tr[blockF[w, i, j] . (-eps) . pauli[[k]]]/2, {a, {"x", "y", "z"}}, {k, {Position[{"x", "y", "z"}, a][[1, 1]]}}]},
    {i, 3}, {j, i, 3}]];
spinMat[m_] := ArrayFlatten[Table[m[i, j], {i, 3}, {j, 3}]];  (* 6 x 6 with modes (site, spin) *)
SeedRandom[11];
rs := RandomReal[{-1, 1}, {3, 3}];
symm[m_] := m + Transpose[m]; anti[m_] := m - Transpose[m];
VerificationTest[Module[{w = gibbs[hA, hB], ch, gen = GaussianState[]},
   ch = channels[w];
   Max[Abs[(WickVEV[#, gen] /. ch) - Tr[w . toMatrix[#]] & /@ tests]] < 10^-10], True,
  TestID -> "general state in channels = exact Fock-space values"]
VerificationTest[Module[{a = symm[rs], b = symm[rs], w, ch},
   w = gibbs[spinMat[a[[#1, #2]] IdentityMatrix[2] &], spinMat[b[[#1, #2]] eps &]/2];
   ch = Association[channels[w]];
   Max[Abs[Values[KeySelect[ch, Length[#] == 3 &]]]] < 10^-10], True, TestID -> "SU(2) Hamiltonian: no spin channels"]
VerificationTest[Module[{a = symm[rs], az = symm[rs], b = rs, w, ch},
   w = gibbs[spinMat[a[[#1, #2]] IdentityMatrix[2] + az[[#1, #2]] pauli[[3]] &], spinMat[b[[#1, #2]] {{0, 1}, {0, 0}} &]];
   ch = Association[channels[w]];
   Max[Abs[Values[KeySelect[ch, MatchQ[#, _[_, _, "x" | "y"]] &]]]] < 10^-10], True, TestID -> "Sz-conserving Hamiltonian: only z channels"]
VerificationTest[Module[{a = symm[rs], bs = Table[anti[rs], 3], p = symm[rs], w, ch},
   w = gibbs[spinMat[a[[#1, #2]] IdentityMatrix[2] + I Sum[bs[[k, #1, #2]] pauli[[k]], {k, 3}] &], spinMat[p[[#1, #2]] eps &]/2];
   ch = Association[channels[w]];
   Max[Abs[Join[Im[Values[KeySelect[ch, Length[#] == 2 &]]], Re[Values[KeySelect[ch, Length[#] == 3 &]]]]]] < 10^-10], True,
  TestID -> "time-reversal invariant Hamiltonian: real scalar channels, imaginary spin channels"]

(* ---------------------------------------------------------------- Hermitian conjugate *)
VerificationTest[Dagger[2 I J cDag[1, up] ** c[2, dn] + \[Chi][1, 2] SpOp[1] ** Bdag[2]],
  -2 I J cDag[2, dn] ** c[1, up] + Conjugate[\[Chi][1, 2]] B[2] ** SmOp[1], TestID -> "Dagger: order, operators, coefficients"]
VerificationTest[{Dagger[SdotS[1, 2]], Dagger[SdotScrossS[1, 2, 3]], Dagger[t Exp[I \[CurlyPhi]] c[1, up]], Dagger[g c[1, up], "Complex" -> {g}]},
  {SdotS[1, 2], SdotScrossS[1, 2, 3], t Exp[-I \[CurlyPhi]] cDag[1, up], Conjugate[g] cDag[1, up]}, TestID -> "Dagger: Hermitian operators, real and complex symbols"]
VerificationTest[With[{a = x[1, 2] cDag[1, up] ** c[2, dn] + 3 SzOp[1], b = SpOp[2] ** c[3, up] + I DOp[1]},
   AllTrue[Values[OperatorCoefficients[Dagger[a ** b] - Dagger[b] ** Dagger[a]]], Simplify[#] === 0 &]], True, TestID -> "(AB)^dagger = B^dagger A^dagger"]
VerificationTest[{HermitianOperatorQ[MeanField[SdotS[1, 2], pairing]], HermitianOperatorQ[MeanField[-Jc SdotScrossS[1, 2, 3], pairing]],
   HermitianOperatorQ[cDag[1, up] ** c[2, up]], HermitianOperatorQ[PlusHC[t cDag[1, up] ** c[2, up]]]}, {True, True, False, True},
  TestID -> "mean-field Hamiltonians are Hermitian; PlusHC"]
VerificationTest[DefineOperator[HopT, 2, (cDag[#1, up] ** c[#2, up] + cDag[#1, dn] ** c[#2, dn]) &, Subscript["T", Row[{#1, #2}]] &];
  {ToSpinChargeBasis[cDag[2, up] ** c[1, up] + cDag[2, dn] ** c[1, dn]], ToSpinChargeBasis[Dagger[HopT[1, 2]]]}, {HopT[2, 1], HopT[2, 1]},
  TestID -> "recognition with sites out of order keeps the fermion sign"]

(* ---------------------------------------------------------------- modes with any labels *)
VerificationTest[{NormalOrder[c[1] ** cDag[1]], NormalOrder[c[i] ** cDag[j]]}, {1 - cDag[1] ** c[1], -cDag[j] ** c[i]},
  TestID -> "spinless fermions c[i]"]
VerificationTest[DeclareModeLabels[{"Momentum"}];
  {NormalOrder[c[k] ** cDag[q]],
   Simplify[WickVEV[cDag[k] ** c[k] ** cDag[q] ** c[q], GaussianState[{{cDag[k1_], c[k2_]} :> KroneckerDelta[k1, k2] nk[k1]}]] -
     (nk[k] nk[q] + KroneckerDelta[k, q] nk[k] (1 - nk[k]))]},
  {KroneckerDelta[k, q] - cDag[q] ** c[k], 0}, TestID -> "momentum labels: symbolic deltas are kept"]
VerificationTest[DeclareModeLabels[{"Site"}]; NormalOrder[c[k] ** cDag[q]], -cDag[q] ** c[k], TestID -> "back to site labels"]
VerificationTest[DeclareModeLabels[{"Site", "Orbital", "Spin"}]; NormalOrder[c[1, a, up] ** cDag[1, b, up]],
  KroneckerDelta[a, b] - cDag[1, b, up] ** c[1, a, up], TestID -> "orbital labels: symbolic delta"]
modes8 = Tuples[{{1, 2}, {1, 2}, {up, dn}}]; sp8 = FockSpace[{1, 2}, "Orbitals" -> {1, 2}];
jw8 = Normal[FockMatrix[c @@ #, sp8]] & /@ modes8;
mat8[op_] := Normal[FockMatrix[op, sp8]];
toMatrix8[expr_] := Normal[FockMatrix[expr, sp8]];
VerificationTest[Module[{a8, b8, h8, w8, num, ops},
   SeedRandom[5]; a8 = RandomComplex[{-1 - I, 1 + I}, {8, 8}]; a8 = a8 + ConjugateTranspose[a8]; b8 = RandomComplex[{-1 - I, 1 + I}, {8, 8}]; b8 = b8 - Transpose[b8];
   h8 = Sum[a8[[p, q]] ConjugateTranspose[jw8[[p]]] . jw8[[q]] + b8[[p, q]] ConjugateTranspose[jw8[[p]]] . ConjugateTranspose[jw8[[q]]] +
       Conjugate[b8[[p, q]]] jw8[[q]] . jw8[[p]], {p, 8}, {q, 8}];
   w8 = With[{e = MatrixExp[-0.9 h8]}, e/Tr[e]];
   num = {\[Rho][x__] :> Tr[w8 . mat8[cDag @@ Take[{x}, 3]] . mat8[c @@ Drop[{x}, 3]]], \[Kappa][x__] :> Tr[w8 . mat8[c @@ Take[{x}, 3]] . mat8[c @@ Drop[{x}, 3]]]};
   ops = {cDag[1, 1, up] ** c[2, 2, dn] ** cDag[2, 1, dn] ** c[1, 2, up], cDag[1, 2, up] ** cDag[2, 1, dn] ** c[1, 1, dn] ** c[2, 2, up],
     Sum[cDag[1, o, s] ** c[2, o, s], {o, 2}, {s, {up, dn}}] ** Sum[cDag[2, o, s] ** c[1, o, s], {o, 2}, {s, {up, dn}}]};
   Max[Abs[(WickVEV[#, GeneralState["Anomalous" -> \[Kappa]]] /. num) - Tr[w8 . toMatrix8[#]] & /@ ops]] < 10^-10], True,
  TestID -> "two orbitals with spin: Wick = exact Fock-space values"]

(* ---------------------------------------------------------------- the banner and the update check, without the network *)
(* the answer of GitHub read from a saved one (LatestRelease.json, the release 0.2.0, trimmed; its text is UTF-8); the line
   of a newer release (the next patch of the version loaded); the once-a-day rule in a temporary directory; loading in a
   script, which never checks; the setting typed before the package is loaded *)
VerificationTest[FermiWick`Private`newerVersionQ @@@ {{"0.10.0", "0.9.0"}, {"0.9.0", "0.10.0"}, {"v0.2.0", "0.2.0"}, {"1.0.0", "0.99.99"},
    {"0.2.1", "0.2.0"}, {"v0.3.0-rc1", "0.2.0"}, {"0.3", "0.2.0"}, {"version 0.3.0", "0.2.0"}},
  {True, False, False, True, True, False, False, False}, TestID -> "update check: versions compared number by number"]
latestAnswer = FromCharacterCode[BinaryReadList[FileNameJoin[{testDir, "LatestRelease.json"}]], "UTF-8"];
latestRelease = <|"Version" -> "0.2.0", "Paclet" -> "https://github.com/n4sd4q/FermiWick/releases/download/v0.2.0/FermiWick-0.2.0.paclet",
   "Page" -> "https://github.com/n4sd4q/FermiWick/releases/tag/v0.2.0"|>;
VerificationTest[{Max[ToCharacterCode[latestAnswer]] > 127, FermiWick`Private`latestReleaseOf[latestAnswer]},
  {True, latestRelease}, TestID -> "update check: the latest release read from an answer of GitHub (text beyond ASCII)"]
VerificationTest[FermiWick`Private`latestReleaseOf /@ {
    StringReplace[latestAnswer, "\"draft\": false" -> "\"draft\": true"],
    StringReplace[latestAnswer, "\"prerelease\": false" -> "\"prerelease\": true"],
    StringReplace[latestAnswer, "\"tag_name\": \"v0.2.0\"" -> "\"tag_name\": \"v0.3.0-rc1\""],
    StringReplace[latestAnswer, "github.com/n4sd4q/FermiWick/releases/download" -> "example.com/releases/download"],
    StringReplace[latestAnswer, "/FermiWick-0.2.0.paclet\"" -> "/a\\\"]; Quit[]; (\\\".paclet\""],
    "{\"message\": \"API rate limit exceeded\", \"documentation_url\": \"https://docs.github.com/rest\"}", "<html>offline</html>", $Failed},
  ConstantArray[$Failed, 8], TestID -> "update check: no release for a draft, a pre-release, a malformed tag, a paclet elsewhere or with a quote, an error"]
nextVersion = StringRiffle[ToString /@ (FermiWick`Private`versionNumbers[$FermiWickVersion] + {0, 0, 1}), "."];
newerRelease = FermiWick`Private`latestReleaseOf[StringReplace[latestAnswer, "0.2.0" -> nextVersion]];
VerificationTest[With[{paclet = "https://github.com/n4sd4q/FermiWick/releases/download/v" <> nextVersion <> "/FermiWick-" <> nextVersion <> ".paclet"},
   {FermiWick`Private`newerVersionQ[newerRelease["Version"], $FermiWickVersion],
    FermiWick`Private`newerVersionQ[latestRelease["Version"], $FermiWickVersion],
    Cases[FermiWick`Private`loadBanner[], Column[lines_List, ___] :> Length[lines], Infinity],
    Cases[FermiWick`Private`loadBanner[newerRelease], Column[lines_List, ___] :> Length[lines], Infinity],
    Cases[FermiWick`Private`updateNotice[newerRelease], Hyperlink[_, url_] :> url, Infinity] === {"https://github.com/n4sd4q/FermiWick/releases/tag/v" <> nextVersion},
    Cases[FermiWick`Private`updateNotice[newerRelease], Defer[e_] :> HoldForm[e], Infinity] === {HoldForm[PacletInstall[paclet]]},
    FermiWick`Private`updateNoticeText[newerRelease] === "FermiWick " <> nextVersion <> " is available (" <> $FermiWickVersion <>
      " is loaded): PacletInstall[\"" <> paclet <> "\"]"}],
  {True, False, {2}, {3}, True, True, True}, TestID -> "update check: the line of a newer release under the banner, and as text"]
VerificationTest[Module[{dir = CreateDirectory[], r},
   r = {FermiWick`Private`updateCheckDueQ[FileNameJoin[{dir, "a"}], "2026-11-10"], FermiWick`Private`updateCheckDueQ[FileNameJoin[{dir, "a"}], "2026-11-10"],
     FermiWick`Private`updateCheckDueQ[FileNameJoin[{dir, "a"}], "2026-11-11"], StringTrim[ReadString[FileNameJoin[{dir, "a", "LastUpdateCheck.txt"}]]],
     Quiet[FermiWick`Private`updateCheckDueQ[FileNameJoin[{dir, "a", "LastUpdateCheck.txt", "b"}], "2026-11-10"]]};
   DeleteDirectory[dir, DeleteContents -> True]; r],
  {True, False, True, "2026-11-11", True}, TestID -> "update check: at most once a day; a date that cannot be kept means a check"]
(* in fresh kernels (scripts): loading prints nothing and asks nothing; $FermiWickUpdateCheck = False typed before loading is the
   package's setting and the symbol typed is removed; loading again keeps the setting *)
freshKernel[code_String] := Module[{file = FileNameJoin[{$TemporaryDirectory, "fw_check_" <> ToString[RandomInteger[10^9]] <> ".wl"}], out},
   Export[file, "PacletDirectoryLoad[" <> ToString[FileNameJoin[{root, "FermiWick"}], InputForm] <> "];\n" <> code, "Text"];
   out = RunProcess[{First[$CommandLine], "-script", file}, "StandardOutput"]; DeleteFile[file]; StringTrim[out]];
VerificationTest[{freshKernel["Needs[\"FermiWick`\"];"],
   freshKernel["$FermiWickUpdateCheck = False;\nNeeds[\"FermiWick`\"];\nPrint[ToString[{FermiWick`$FermiWickUpdateCheck, Names[\"Global`$FermiWickUpdateCheck\"]}, InputForm]];"],
   freshKernel["Needs[\"FermiWick`\"];\n$FermiWickUpdateCheck = False;\nGet[\"FermiWick`\"];\nPrint[ToString[$FermiWickUpdateCheck, InputForm]];"]},
  {"", "{False, {}}", "False"}, TestID -> "update check: a script loads the package printing nothing; the setting typed before loading is kept"]
