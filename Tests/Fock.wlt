(* Tests of FermiWick`Fock`: the operators of FermiWick as matrices, sectors, exact diagonalisation.
   Run with  math -script Tests/RunTests.wls  (or TestReport["Tests/Fock.wlt"]). *)

testDir = DirectoryName[$TestFileName]; root = ParentDirectory[testDir];
PacletDirectoryLoad[FileNameJoin[{root, "FermiWick"}]];
Needs["FermiWick`Fock`"]; Needs["FermiWick`Gutzwiller`"];
VerificationTest[MemberQ[$Packages, "FermiWick`Fock`"], True, TestID -> "load Fock"]
(* the engine below FockMatrix: strings of mode indices, +alpha for c^dagger_alpha and -alpha for c_alpha *)
ladder = FermiWick`Fock`Private`ladder; fockOperator = FermiWick`Fock`Private`fockOperator;

up = SpinUp; dn = SpinDown;
hubbard[sites_, bonds_, t_, u_] := -t Sum[cDag[i, s] ** c[j, s] + cDag[j, s] ** c[i, s], {b, bonds}, {i, {b[[1]]}}, {j, {b[[2]]}}, {s, $Spins}] +
   u Sum[DOp[i], {i, sites}];
ring[l_] := Table[{i, Mod[i, l] + 1}, {i, l}];

(* ---------------------------------------------------------------- the matrices of the ladder operators *)
(* the convention, written out once with Kronecker products: mode 1 is the most significant bit, c_m = Z (x) ... (x) Z (x) sigma^- (x) 1 ... *)
VerificationTest[Module[{sp = FockSpace[{1, 2, 3}], jw},
   jw = Table[KroneckerProduct @@ Join[ConstantArray[{{1, 0}, {0, -1}}, m - 1], {{{0, 1}, {0, 0}}}, ConstantArray[IdentityMatrix[2], 6 - m]], {m, 6}];
   And @@ Table[Normal[FockMatrix[sp["Modes"][[m]], sp]] == jw[[m]] && Normal[FockMatrix[cDag @@ sp["Modes"][[m]], sp]] == Transpose[jw[[m]]], {m, 6}]],
  True, TestID -> "c and cDag are the Jordan-Wigner matrices"]
VerificationTest[Module[{sp = FockSpace[{1, 2, 3}], id = IdentityMatrix[64], ops, M},
   ops = sp["Modes"]; M[x_] := Normal[FockMatrix[x, sp]];
   {Max@Table[Max@Abs[M[a] . M[cDag @@ b] + M[cDag @@ b] . M[a] - If[a === b, id, 0 id]], {a, ops}, {b, ops}],
    Max@Table[Max@Abs[M[a] . M[b] + M[b] . M[a]], {a, ops}, {b, ops}], Max@Abs[Normal[FockMatrix[Sum[nTotOp[i], {i, 3}], sp]] - Total[M /@ (nOp @@@ ops)]]}],
  {0, 0, 0}, TestID -> "canonical anticommutators, N"]
VerificationTest[With[{sp = FockSpace[{1, 2}, "Orbitals" -> {"a", "b"}]}, {sp["Count"], sp["Modes"][[{1, 2, 3, 8}]], ModeIndex[sp, c[2, "a", dn]],
    ModeIndex[sp, 2, "b", up], Quiet[ModeIndex[sp, c[2, dn]]]}],
  {8, {c[1, "a", up], c[1, "a", dn], c[1, "b", up], c[2, "b", dn]}, 6, 7, $Failed}, TestID -> "modes: site, orbital, spin fastest"]
VerificationTest[With[{sp = FockSpace[{"x", "y"}, "Spins" -> None]}, {sp["Modes"], Normal[FockMatrix[cDag["x"] ** c["y"], sp]]}],
  {{c["x"], c["y"]}, {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 0}}}, TestID -> "spinless modes"]

(* ---------------------------------------------------------------- FockMatrix reads FermiWick *)
sp2 = FockSpace[{1, 2}]; M2[x_] := Normal[FockMatrix[x, sp2]];
ops2 = Join[sp2["Modes"], cDag @@@ sp2["Modes"]];
VerificationTest[SeedRandom[1]; And @@ Table[With[{x = Total[RandomInteger[{-2, 2}, 3] (NonCommutativeMultiply @@ RandomChoice[ops2, #] & /@ {2, 3, 4})]},
    M2[x] === M2[NormalOrder[x]] && M2[x] === M2[ToSpinChargeBasis[x]] && M2[x] === M2[GroupLocalOperators[x]]], {40}],
  True, TestID -> "the matrix of an expression = the matrix of its normal order, of its spin-charge form"]
VerificationTest[ToSpinChargeBasis[c[2, up] ** c[1, dn] - 2 cDag[1, dn] ** c[2, up] ** c[2, up] ** cDag[1, dn]], -c[1, dn] ** c[2, up],
  TestID -> "a product with a vanishing local factor vanishes (ToSpinChargeBasis)"]
VerificationTest[{M2[SpOp[1] ** cDag[2, up] ** c[1, dn]] == M2[SpOp[1]] . M2[cDag[2, up]] . M2[c[1, dn]],
   M2[Dagger[3 I cDag[1, up] ** c[2, dn] ** Bdag[2]]] == ConjugateTranspose[M2[3 I cDag[1, up] ** c[2, dn] ** Bdag[2]]],
   M2[OperatorCommutator[SpOp[1], SmOp[1]]] == 2 M2[SzOp[1]]},
  {True, True, True}, TestID -> "products, Dagger and commutators as matrices"]
VerificationTest[{M2[SzOp[1]] == (M2[nOp[1, up]] - M2[nOp[1, dn]])/2, M2[nTotOp[2]] == M2[nOp[2]],
   M2[SdotS[1, 2]] == Sum[M2[SpinOperator[1, a]] . M2[SpinOperator[2, a]], {a, 3}], M2[DOp[1, 2]] == M2[nOp[1, up] ** nOp[1, dn]] . M2[DOp[2]],
   M2[EmptyProjector[1]] == M2[(1 - nOp[1, up]) ** (1 - nOp[1, dn])], M2[B[1]] == ConjugateTranspose[M2[Bdag[1]]],
   M2[SingleOccupationProjector[2]] == M2[nOp[2] - 2 DOp[2]]},
  ConstantArray[True, 7], TestID -> "composite and defined operators"]
VerificationTest[Simplify[M2[GutzwillerProjector[1] ** GutzwillerProjector[1]] - M2[2 (Cos[\[Theta]]^2 SingleOccupationProjector[1] +
       Sin[\[Theta]]^2 (EmptyProjector[1] + DOp[1]))]] === ConstantArray[0, {16, 16}],
  True, TestID -> "Gutzwiller projector, symbolic angle"]
VerificationTest[Module[{sp = FockSpace[Range[3]], x},
   x = SdotScrossS[1, 2, 3];
   {Max@Abs[Normal[FockMatrix[x, sp] - ConjugateTranspose[FockMatrix[x, sp]]]], Max@Abs[Normal[FockMatrix[OperatorCommutator[x, Sum[SzOp[i], {i, 3}]], sp]]]}],
  {0, 0}, TestID -> "chirality: Hermitian, commutes with S^z"]
VerificationTest[Module[{sp = FockSpace[{1, 2}], s1, s2},
   s1 = FockSector[sp, "N" -> 1, "Sz" -> 1/2]; s2 = FockSector[sp, "N" -> 2, "Sz" -> 0];
   {Dimensions[FockMatrix[cDag[1, dn], s1]], Dimensions[FockMatrix[cDag[2, dn] ** c[1, up], {s1, FockSector[sp, "N" -> 1, "Sz" -> -1/2]}]],
    With[{s0 = FockSector[sp, "N" -> 0]}, Normal[FockMatrix[Bdag[1], {s0, s2}]] == Normal[FockMatrix[Bdag[1], sp]][[s2["Indices"] + 1, s0["Indices"] + 1]]],
    Normal[FockMatrix[5, s2]] == 5 IdentityMatrix[4], Normal[FockMatrix[0 cDag[1, up], s2]] == ConstantArray[0, {4, 4}]}],
  {{4, 2}, {2, 2}, True, True, True}, TestID -> "FockMatrix between sectors, constants"]
VerificationTest[Module[{sp = FockSpace[{1, 2}]}, Simplify[Normal[FockMatrix[hubbard[{1, 2}, {{1, 2}}, tt, uu] - mu Sum[nOp[i], {i, 2}], sp]] -
     Normal[fockOperator[sp, Join[Flatten[Table[{{-tt, {ModeIndex[sp, i, s], -ModeIndex[sp, 3 - i, s]}}}, {i, 2}, {s, $Spins}], 2],
        Table[{uu, {ModeIndex[sp, i, up], ModeIndex[sp, i, dn], -ModeIndex[sp, i, dn], -ModeIndex[sp, i, up]}}, {i, 2}],
        Flatten[Table[{-mu, {ModeIndex[sp, i, s], -ModeIndex[sp, i, s]}}, {i, 2}, {s, $Spins}], 1]]]]] === ConstantArray[0, {16, 16}]],
  True, TestID -> "symbolic coefficients: FockMatrix = FockOperator of the strings"]
VerificationTest[Quiet[{FockMatrix[c[3, up], sp2], FockMatrix[nOp[1, s], sp2], FockSpace[{}], FockSpace[{1}, "Spins" -> {}],
    FockExpectation[SzOp[1], {1, 0}, sp2]}],
  ConstantArray[$Failed, 5], TestID -> "FockMatrix: modes outside the space, symbolic spins"]
VerificationTest[M2[c[1, up] ** foo[2]] == foo[2] M2[c[1, up]], True, TestID -> "FockMatrix: other heads are coefficients, as in FermiWick"]

(* ---------------------------------------------------------------- Wick's theorem, exactly *)
(* (the Gibbs states of quadratic Hamiltonians with pairing are in FermiWick.wlt) the contractions of a free ground state make a
   Gaussian state, whose WickVEV must give the exact expectation values *)
VerificationTest[Module[{l = 6, sp, sec, h, e, v, rules, gs},
   sp = FockSpace[Range[l]]; sec = FockSector[sp, "N" -> 6, "Sz" -> 0];                     (* half filling: a closed shell, unique ground state *)
   h = FockMatrix[hubbard[Range[l], ring[l], 1, 0], sec];
   {e, v} = LowestStates[h, 1];
   rules = Flatten[Table[{cDag[i, s], c[j, t]} -> FockExpectation[cDag[i, s] ** c[j, t], First[v], sec], {i, l}, {j, l}, {s, $Spins}, {t, $Spins}], 3];
   gs = GaussianState[rules];                                              (* <c c> = 0, <c cDag> from the anticommutator *)
   Max[Abs[(WickVEV[#, gs] - FockExpectation[#, First[v], sec]) & /@ {SdotS[1, 2], DOp[1], SdotS[1, 4], nTotOp[1] ** nTotOp[3]}]] < 10^-10],
  True, TestID -> "free ground state: Wick's theorem with the exact contractions"]

(* ---------------------------------------------------------------- sectors *)
VerificationTest[FermiWick`Fock`Private`combinations[10, 4], Select[Range[0, 2^10 - 1], DigitCount[#, 2, 1] == 4 &], TestID -> "enumeration"]
VerificationTest[
  With[{sp = FockSpace[Range[4]]},
   And @@ Table[With[{s = FockSector[sp, Sequence @@ cs]}, FermiWick`Fock`Private`positionsIn[s, s["Indices"]] === Range[s["Dimension"]]],
     {cs, {{"N" -> 3}, {"N" -> 4, "Sz" -> 0}, {"Sz" -> 1}, {"Counts" -> {{{1, 2, 3}, 2}}}, {"N" -> 5, "Sz" -> 1/2, "Counts" -> {{{1, 2}, 1}}}}}]],
  True, TestID -> "positions"]
VerificationTest[
  With[{sp = FockSpace[Range[3], "Orbitals" -> {"a", "b"}], ops = {{2, -5}, {1, 3, -6, -2}, {7, -7}, {4, -9}, {}, {3, -3, 1, -1}, {2, 5, -2, -5}}},
   And @@ Flatten@Table[With[{s = FockSector[sp, Sequence @@ cs]},
       Table[Max[Abs[SectorBlock[ladder[sp, o], s] - ladder[{s, s}, o]]] == 0, {o, ops}]],
     {cs, {{"N" -> 4}, {"N" -> 5, "Sz" -> 1/2}, {"Counts" -> {{{1, 2, 3, 4}, 2}}}, {"N" -> 6, "Sz" -> 0, "Counts" -> {{{1, 2}, 1}}}}}]],
  True, TestID -> "sector operators are blocks (diagonal strings too)"]
VerificationTest[
  With[{sp = FockSpace[Range[3]]}, With[{from = FockSector[sp, "N" -> 2, "Sz" -> 0], to = FockSector[sp, "N" -> 4, "Sz" -> 1]},
    Max[Abs[ladder[sp, {1, 3}][[to["Indices"] + 1, from["Indices"] + 1]] - ladder[{from, to}, {1, 3}]]] == 0 &&
     Dimensions[ladder[from, {1, 3}]] == {to["Dimension"], from["Dimension"]}]],
  True, TestID -> "between sectors"]
VerificationTest[
  With[{sp = FockSpace[Range[3]]}, With[{terms = {{0.3, {1, -3}}, {0.3, {3, -1}}, {-1.2 I, {2, -4}}, {1.2 I, {4, -2}}, {2., {1, -1, 2, -2}}, {0.5, {}}}},
    Max[Abs[fockOperator[sp, terms] - Total[#[[1]] ladder[sp, #[[2]]] & /@ terms]]] < 10^-14]],
  True, TestID -> "FockOperator"]
VerificationTest[With[{sp = FockSpace[{1, 2}]}, {FockSector[sp, "Counts" -> {{{c[1, up], c[1, dn]}, 2}}]["Dimension"],
    FockSector[sp, "Counts" -> {{{c[1, up], c[1, dn]}, 2}}] === FockSector[sp, "Counts" -> {{{1, 2}, 2}}], Quiet[FockSector[sp, "Counts" -> {{{c[3, up]}, 1}}]]}],
  {4, True, $Failed}, TestID -> "Counts given by ladder operators"]

(* ---------------------------------------------------------------- physics *)
(* Hubbard dimer, t = 1, U = 3, interaction (U/2) sum (n_i - 1)^2 *)
dimerSpace = FockSpace[{1, 2}];
dimer = FockMatrix[hubbard[{1, 2}, {{1, 2}}, 1, 0] + 3/2 Sum[(nTotOp[i] - 1) ** (nTotOp[i] - 1), {i, 2}], dimerSpace];
VerificationTest[Min[Eigenvalues[N@Normal@dimer]], -1., SameTest -> (Abs[#1 - #2] < 10^-12 &), TestID -> "dimer ground state"]
VerificationTest[Sort[Eigenvalues[N@Normal@SpinSquared[FockSector[dimerSpace, "N" -> 2]]]], {0., 0., 0., 2., 2., 2.},
  SameTest -> (Max[Abs[#1 - #2]] < 10^-12 &), TestID -> "dimer S^2"]
VerificationTest[With[{sec = FockSector[dimerSpace, "N" -> 2]},
   Max[Abs[Normal[SpinSquared[sec] - FockMatrix[Sum[SdotS[i, j], {i, 2}, {j, 2}], sec]]]]], 0, TestID -> "SpinSquared = FockMatrix of S.S"]
VerificationTest[With[{g = LehmannGreen[dimer, dimerSpace, dimerSpace["Modes"]]},
   {g["GroundEnergy"], g["GroundSectors"], g["Degeneracy"], Sort@DeleteDuplicates[Round[g["Poles"], 10^-9]], Max[Abs[Total[g["Residues"]] - IdentityMatrix[4]]] < 10^-12}],
  {-1., {{2, 0}}, 1, {-7/2, -3/2, 3/2, 7/2}, True}, SameTest -> (Abs[#1[[1]] - #2[[1]]] < 10^-12 && #1[[2 ;;]] == #2[[2 ;;]] &), TestID -> "dimer Green's function"]
VerificationTest[With[{g = LehmannGreen[dimer, dimerSpace, {c[1, up], c[2, up]}]}, {g["Modes"], Chop[GreenMatrix[g, 0.3 I] - GreenMatrix[LehmannGreen[dimer, dimerSpace, {1, 3}], 0.3 I]]}],
  {{c[1, up], c[2, up]}, {{0, 0}, {0, 0}}}, TestID -> "LehmannGreen with ladder operators"]

(* magnetic plaquette, t = 1, t0 = 0.5, phi = pi/3, U = 5, N = 4 (value checked against an independent code) *)
VerificationTest[Module[{sp = FockSpace[Range[4]], nb, phi = Pi/3, t0 = 0.5, u = 5., h},
   nb[j_] := Mod[j - 1, 4] + 1;
   h = Sum[-Exp[I phi/4] cDag[nb[j], s] ** c[nb[j + 1], s] - Exp[-I phi/4] cDag[nb[j + 1], s] ** c[nb[j], s] - t0 cDag[nb[j], s] ** c[nb[j + 2], s], {j, 4}, {s, $Spins}] +
     u/2 Sum[(nTotOp[j] - 1) ** (nTotOp[j] - 1), {j, 4}];
   Min@Eigenvalues@Normal@FockMatrix[h, FockSector[sp, "N" -> 4, "Sz" -> 0]]],
  -1.974705592851098, SameTest -> (Abs[Re[#1] - #2] < 10^-10 &), TestID -> "magnetic plaquette"]

(* six-site Hubbard ring, U = 4, half filling: dense and Arnoldi agree *)
ring6 = FockMatrix[hubbard[Range[6], ring[6], 1, 4], FockSector[FockSpace[Range[6]], "N" -> 6, "Sz" -> 0]];
VerificationTest[{First@First@LowestStates[ring6, 1], Block[{$FockDenseLimit = 10}, First@First@LowestStates[ring6, 1]]},
  {-3.6687061788729345, -3.6687061788729345}, SameTest -> (Max[Abs[#1 - #2]] < 10^-9 &), TestID -> "Hubbard ring, dense and Arnoldi"]
VerificationTest[Module[{h = FockMatrix[hubbard[Range[6], ring[6], 1, 4], FockSector[FockSpace[Range[6]], "N" -> 6]]},
   Max[Abs[First[LowestStates[h, 6]] - Block[{$FockDenseLimit = 10}, First[LowestStates[h, 6]]]]] < 10^-8],
  True, TestID -> "Arnoldi finds degenerate multiplets"]
VerificationTest[Module[{h = FockMatrix[hubbard[Range[7], ring[7], 1, 4], FockSector[FockSpace[Range[7]], "N" -> 7]], ex},
   ex = Take[Sort@Eigenvalues[N@Normal@h], 8];                                      (* dimension 3432: Arnoldi; four-fold ground state *)
   {Max[Abs[First[LowestStates[h, 4]] - ex[[;; 4]]]] < 10^-8, Max[Abs[First[LowestStates[h, 8]] - ex]] < 10^-8}],
  {True, True}, TestID -> "Arnoldi, four-fold degenerate levels"]

(* ---------------------------------------------------------------- symmetries *)
VerificationTest[Module[{sp = FockSpace[{1, 2}, "Orbitals" -> {"s", "p"}], ph},
   ph = ParticleHoleOperator[sp];
   Max@Table[Max@Abs[ph . FockMatrix[a, sp] . Transpose[ph] - FockMatrix[cDag @@ a, sp]], {a, sp["Modes"]}]],
  0, TestID -> "particle-hole: P c P^-1 = cDag"]
VerificationTest[Module[{sp = FockSpace[{1, 2}], tr, v},
   tr = TimeReversal[sp]; SeedRandom[3]; v = RandomComplex[{-1 - I, 1 + I}, 16];
   Max@Table[Max[Abs[tr["Apply"][FockMatrix[c[i, $Spins[[p]]], sp] . tr["Inverse"][v]] -
        Sum[MatrixExp[I Pi/2 PauliMatrix[2]][[p, q]] FockMatrix[c[i, $Spins[[q]]], sp], {q, 2}] . v]], {i, 2}, {p, 2}]],
  0., SameTest -> (Abs[#1 - #2] < 10^-10 &), TestID -> "time reversal"]
VerificationTest[With[{sp = FockSpace[{1}]}, {Normal[ParticleHoleOperator[sp, 3, 2]] == IdentityMatrix[4],
    Normal[ParticleHoleOperator[sp, 0, 2]] == Normal[ParticleHoleOperator[sp, 1, 2]], Quiet[TimeReversal[FockSpace[{1}, "Spins" -> None]]]}],
  {True, True, $Failed}, TestID -> "symmetries: ranges, spinless spaces"]

(* ---------------------------------------------------------------- input checks *)
sp4 = FockSpace[{1, 2}]; s4 = FockSector[sp4, "N" -> 2];
VerificationTest[TimeConstrained[{Quiet[ladder[s4, {5}]], Quiet[ladder[sp4, {0}]], Quiet[fockOperator[sp4, {{1, {7, -1}}}]], Quiet[ladder[{s4, FockSector[FockSpace[{1, 2, 3}], "N" -> 2]}, {1}]]}, 10],
  ConstantArray[$Failed, 4], TestID -> "invalid mode indices"]
VerificationTest[Module[{s = FockSector[sp4, "N" -> 2, "Counts" -> {{{1, 2}, 1}}], h, terms},
   terms = {{1, {1, -2}}, {1, {2, -1}}, {1, {3, -4}}, {1, {4, -3}}, {2, {1, -1, 3, -3}}};
   h = fockOperator[s, terms];
   {Dimensions[h] == {s["Dimension"], s["Dimension"]}, Dimensions[ladder[s, {1, -3}]] == {FockSector[sp4, "N" -> 2, "Counts" -> {{{1, 2}, 2}}]["Dimension"], s["Dimension"]},
    Max[Abs[h - SectorBlock[fockOperator[sp4, terms], s]]] == 0}],
  {True, True, True}, TestID -> "Counts sectors"]
VerificationTest[{Block[{$FockFullLimit = 2}, Quiet[FockSector[FockSpace[{"p", "q"}]]]], FockSectorQ[FockSector[FockSpace[{"p", "q"}]]]},
  {$Failed, True}, TestID -> "failed sectors are not cached"]
VerificationTest[Module[{hx, g1, g2, g3},
   hx = dimer + 0.7 FockMatrix[SpOp[1] + SmOp[1], dimerSpace];          (* a transverse field on site 1: conserves N, not Sz *)
   g1 = LehmannGreen[hx, dimerSpace, Range[4]]; g2 = LehmannGreen[hx, dimerSpace, Range[4], "UseSz" -> False];
   g3 = LehmannGreen[SectorBlock[hx, #] &, dimerSpace, Range[4], "UseSz" -> False];
   {Max[Abs[Sort[g1["Poles"]] - Sort[g2["Poles"]]]] < 10^-10, Max[Abs[Total[g1["Residues"]] - IdentityMatrix[4]]] < 10^-10,
    Max[Abs[Sort[g3["Poles"]] - Sort[g2["Poles"]]]] < 10^-10}],
  {True, True, True}, TestID -> "LehmannGreen without Sz conservation, H as a function"]
VerificationTest[Module[{a = FockSector[sp4, "N" -> 1], b = FockSector[sp4, "N" -> 2]},
   {Dimensions[fockOperator[{a, b}, {{1, {1}}, {2, {3}}}]], Dimensions[fockOperator[a, {}]], Dimensions[FockMatrix[Sum[nOp @@ m, {m, sp4["Modes"]}], {a, a}]],
    Length[ladder[FockSector[sp4, "N" -> 1, "Sz" -> 1], {}]], Quiet[ladder[{a, FockSector[FockSpace[{1, 2, 3}], "N" -> 2]}, {1}]]}],
  {{6, 4}, {4, 4}, {4, 4}, 0, $Failed}, TestID -> "explicit pairs, empty terms, empty sectors"]
VerificationTest[{Quiet[SpinSquared[FockSpace[{1, 2}, "Spins" -> None]]], Quiet[BasisState[sp4, {1, 0}]], Quiet[StateToKet[s4, {1, 0}]],
   FockSector[sp4, "Counts" -> {{{1, 1}, 1}}]["Dimension"] == FockSector[sp4, "Counts" -> {{{1}, 1}}]["Dimension"],
   Quiet[FockSector[sp4, "Counts" -> {{{9}, 1}}]], Quiet[OccupationBasis[3]], Quiet[FockSpace[{1, 2}, {1}]],
   Block[{$FockFullLimit = 2}, Quiet[{StateToKet[FockSpace[{"u", "v"}], {1}], FockMatrix[nOp["u", up], FockSpace[{"u", "v"}]]}]]},
  {$Failed, $Failed, $Failed, True, $Failed, $Failed, $Failed, {$Failed, $Failed}}, TestID -> "input checks"]
VerificationTest[{Quiet[LehmannGreen[IdentityMatrix[3], dimerSpace, Range[4]]], Quiet[LehmannGreen[dimer, FockSpace[{1, 2}, "Spins" -> None], {1}, "UseSz" -> True]],
   Quiet[LehmannGreen[fockOperator[#, {{1, {1, 2}}}] &, dimerSpace, Range[4]]], Quiet[LehmannGreen[dimer, dimerSpace, {c[3, up]}]]},
  {$Failed, $Failed, $Failed, $Failed}, TestID -> "LehmannGreen checks its Hamiltonian and modes"]
VerificationTest[With[{sec = FockSector[dimerSpace, "N" -> 2, "Sz" -> 0]}, StateToKet[sec, {0, 1, -1, 0}/Sqrt[2]]],
  Ket[{{"\[DownArrow]"}, {"\[UpArrow]"}}]/Sqrt[2] - Ket[{{"\[UpArrow]"}, {"\[DownArrow]"}}]/Sqrt[2], TestID -> "StateToKet"]
