(* ::Package:: *)

(* FermiWick`Fock`: the operators of FermiWick as matrices on a fermionic Fock space, symmetry sectors, exact
   diagonalisation and Lehmann Green's functions.

   Conventions
   -----------
   The modes of a space are the ladder operators c[l1, l2, ...] of FermiWick: FockSpace[sites] has the modes c[i, s]
   (s = SpinUp, SpinDown), FockSpace[sites, "Orbitals" -> orbs] the modes c[i, o, s], and with "Spins" -> None the modes
   c[i] or c[i, o]. They are numbered site by site, then orbital by orbital, spin fastest:
       alpha = (site - 1) no nsp + (orbital - 1) nsp + spin ,
   the order of Tuples[{sites, orbitals, spins}] (and of ContractionTable).
   A basis state is a string of occupations (n_1, ..., n_M); its index in the full space is
       k = sum_alpha n_alpha 2^(M - alpha)          (0-based; mode 1 is the most significant bit)
   which is the order of Tuples[{0, 1}, M] and of KroneckerProduct of 2x2 blocks in the basis {|0>, |1>}.
   Jordan-Wigner sign: c_alpha |n> = (-1)^(sum_{beta < alpha} n_beta) |n - e_alpha>.

   FockMatrix reads an expression of FermiWick term by term: its composite and defined operators are written with c and
   cDag (ToLadderBasis), and every product, in the order written (no normal ordering), is a string of ladder operators
   applied to the basis states. All operators are SparseArray objects. A FockSector (fixed particle number, fixed S^z,
   or any set of occupation constraints) carries the list of its basis states, in increasing order of k, and every
   operator is built directly between two sectors, without ever forming the 2^M x 2^M matrix. Sectors with fixed N are
   enumerated bit by bit and indexed through the combinatorial number system, so that time and memory scale with the
   dimension of the sector, not with 2^M. All bit manipulations are vectorised (shifts, masks and table look-ups on
   packed arrays). *)

BeginPackage["FermiWick`Fock`", {"FermiWick`"}];

FockSpace::usage = "FockSpace[sites] is the Fock space of spin-1/2 fermions on the sites, with the modes c[i, s]; \
FockSpace[sites, \"Orbitals\" -> {o1, ...}] has the modes c[i, o, s], and \"Spins\" -> None gives spinless modes c[i] (or \
c[i, o]); \"Spins\" -> {s1, s2} names the two spins (s1 is up). The modes are numbered site by site, then orbital, spin \
fastest. space[\"Modes\"] is the list of the modes, space[\"Dimension\"] the dimension 2^M.";
FockSpaceQ::usage = "FockSpaceQ[expr] is True for a FockSpace object.";
FockSector::usage = "FockSector[space, \"N\" -> n, \"Sz\" -> sz, \"Counts\" -> {{modes, n}, ...}] is the subspace of the basis states that \
satisfy all the given constraints: the particle number, the z component of the spin, the number of particles in a set of modes \
(given as ladder operators c[...] or as mode indices). FockSector[space] is the whole space. sector[\"Dimension\"] is its dimension.";
FockSectorQ::usage = "FockSectorQ[expr] is True for a FockSector object.";
ModeIndex::usage = "ModeIndex[space, c[l1, l2, ...]] is the index of a mode of the space (cDag[...] gives the same index); \
ModeIndex[space, l1, l2, ...] is the same with the labels. It gives $Failed with a message for a mode that is not in the space.";
OccupationBasis::usage = "OccupationBasis[space or sector] gives the occupation vectors of the basis states, in order.";
BasisState::usage = "BasisState[space or sector, occupations] gives the basis vector with the given occupation vector.";
StateToKet::usage = "StateToKet[space or sector, vector] writes a state as a combination of kets labelled by the local occupations \
(0, \[UpArrow], \[DownArrow], \[UpArrow]\[DownArrow] for spin 1/2).";

FockMatrix::usage = "FockMatrix[expr, target] is the matrix of an operator of FermiWick: ladder operators c, cDag, composite operators \
(nOp, SzOp, SpOp, DOp, Bdag, ...), operators defined with DefineOperator, products with ** in the order written, and commuting \
coefficients (numbers or symbols). target is a space, a sector (the image sector is found from the change of N, Sz and of every \
\"Counts\" constraint) or a pair {fromSector, toSector}. The labels of every operator must be those of modes of the space.";
FockExpectation::usage = "FockExpectation[expr, v, sector] is the expectation value <v| expr |v> of an operator of FermiWick in the \
normalised state v of the sector (or space).";
SpinSquared::usage = "SpinSquared[target] is S^2 (on the full space or on a sector with fixed N and/or Sz; two spins).";
ParticleHoleOperator::usage = "ParticleHoleOperator[space, first, last] is the unitary P of the particle-hole transformation of the modes \
first..last (all by default): P c_alpha P^-1 = c^dagger_alpha for those modes, with P = (sigma^x on each of them) times \
exp(-I Pi Sum_alpha (alpha - 1) n_alpha) on the image.";
TimeReversal::usage = "TimeReversal[space] gives <|\"U\" -> U, \"Apply\" -> f, \"Inverse\" -> finv|>: the antiunitary T = U K for spin-1/2 \
fermions, with f[v] = U.Conjugate[v].";
SectorBlock::usage = "SectorBlock[H, sector] extracts from a full-space matrix H the block of the sector.";
LowestStates::usage = "LowestStates[H, k] gives {energies, vectors} of the k lowest eigenstates of the Hermitian matrix H (dense below \
$FockDenseLimit, Arnoldi above; for k > 1 Arnoldi is run from several starting vectors and H is diagonalised in their span, \
so that degenerate multiplets are complete).";
LehmannGreen::usage = "LehmannGreen[H, space, modes] computes the T = 0 single-particle Green's function G_ab(z) = <c_a (z - H + E0)^-1 c_b^dagger> \
+ <c_b^dagger (z + H - E0)^-1 c_a> of the number-conserving Hamiltonian H (full-space matrix, or a function sector -> matrix) for the \
modes a, b in modes (ladder operators c[...] or mode indices), averaged over the degenerate ground multiplet. The result is an \
association with \"GroundEnergy\", \"GroundSectors\", \"Degeneracy\", \"Poles\", \"Residues\". Options: \"N\" (particle numbers to search \
for the ground state, default All), \"UseSz\" (work in (N, Sz) sectors; Automatic checks that a matrix H conserves Sz, and a function H on \
the sector N = M/2), \"Tolerance\" (energy window of the ground multiplet, 10^-9), \"ResidueCutoff\" (smaller residues are dropped, 10^-12).";
GreenMatrix::usage = "GreenMatrix[g, z] evaluates the Green's function g returned by LehmannGreen at the complex frequency z.";
ClearFockCache::usage = "ClearFockCache[] forgets the sectors built so far (they are cached by space and constraints).";
$FockFullLimit::usage = "$FockFullLimit is the largest number of modes for which sectors without a fixed particle number, which need the \
whole space, are built (default 24).";
$FockDenseLimit::usage = "$FockDenseLimit is the largest dimension of the matrices that LowestStates diagonalises densely (default 3000).";

Begin["`Private`"];
$FockDenseLimit = 3000;

(* ================================================================= spaces ================================================================= *)
Options[FockSpace] = {"Orbitals" -> None, "Spins" -> {"\[UpArrow]", "\[DownArrow]"}};
FockSpace[sites_List, opts : OptionsPattern[]] := Module[{orbs = OptionValue["Orbitals"], spins = OptionValue["Spins"], lists, m},
   If[! (orbs === None || ListQ[orbs] && orbs =!= {}) || ! (spins === None || ListQ[spins] && spins =!= {}) || sites === {},
    Message[FockSpace::args, {sites, opts}]; Return[$Failed]];
   lists = DeleteCases[{sites, orbs, spins}, None];
   m = Times @@ Length /@ lists;
   FockSpace[<|"Sites" -> sites, "Orbitals" -> orbs, "Spins" -> spins,
     "Shape" -> {Length[sites], If[orbs === None, 1, Length[orbs]], If[spins === None, 1, Length[spins]]}, "Count" -> m,
     "Dimension" -> 2^m, "Labels" -> Tuples[lists], "Modes" -> Apply[c, Tuples[lists], {1}], "Weights" -> 2^(m - Range[m]),
     "Index" -> AssociationThread[Tuples[lists], Range[m]]|>]];
FockSpace[a_Association][key_String] := a[key];
FockSpaceQ[FockSpace[_Association]] = True; FockSpaceQ[_] = False;
FockSpace /: MakeBoxes[s : FockSpace[a_Association], fmt_] :=
  BoxForm`ArrangeSummaryBox[FockSpace, s, None, {BoxForm`SummaryItem[{"Modes: ", a["Count"]}], BoxForm`SummaryItem[{"Dimension: ", a["Dimension"]}]},
   {BoxForm`SummaryItem[{"Sites: ", a["Sites"]}], BoxForm`SummaryItem[{"Orbitals: ", a["Orbitals"]}], BoxForm`SummaryItem[{"Spins: ", a["Spins"]}]}, fmt];
FockSpace::args = "FockSpace[sites, \"Orbitals\" -> orbitals, \"Spins\" -> spins] with a non-empty list of sites, and None or a non-empty \
list for the orbitals and the spins; got `1`.";

(* the labels of a mode are those of c in FermiWick, matched literally *)
modeOf[a_Association, labels_List] := Lookup[a["Index"], Key[labels], Missing[]];
ModeIndex[FockSpace[a_Association], (c | cDag)[labels___]] := ModeIndex[FockSpace[a], labels];
ModeIndex[FockSpace[a_Association], labels___] := With[{p = modeOf[a, {labels}]},
   If[MissingQ[p], Message[ModeIndex::label, c[labels]]; $Failed, p]];
ModeIndex::label = "`1` is not a mode of the space.";
(* mode indices or ladder operators -> indices; a failure throws *)
modeIndices[sp_, l_List] := Map[If[IntegerQ[#], #, With[{p = ModeIndex[sp, #]}, If[p === $Failed, Throw[$Failed, fwFail], p]]] &, l];

(* ================================================================= checks ================================================================= *)
modesOK[m_, ops_List] := VectorQ[ops, IntegerQ[#] && # != 0 && Abs[#] <= m &];
checkOps[head_, m_, ops_] := If[modesOK[m, ops], True, Message[head::mode, ops, m]; False];
FockMatrix::mode = "The operator string `1` must contain nonzero integers with absolute value at most the number of modes, `2`.";
LehmannGreen::mode = "Mode indices `1` must be integers between 1 and `2`.";
posModesOK[head_, m_, al_List] := If[VectorQ[al, IntegerQ[#] && 1 <= # <= m &], True, Message[head::mode, al, m]; False];
spinCheck[head_, t_] := If[spaceOf[t]["Shape"][[3]] == 2, True, Message[head::spin]; False];
SpinSquared::spin = "Spin operators need exactly two spins (the first is up).";

(* ================================================================= bits =================================================================== *)
(* BitGet and DigitCount are not vectorised; shifts, masks and table look-ups on packed arrays are. *)
pop16 = Module[{p = {0}}, Do[p = Join[p, p + 1], {16}]; Developer`ToPackedArray[p]];
popCount[k_] := Module[{n = pop16[[BitAnd[k, 65535] + 1]], r = BitShiftRight[k, 16]},
   While[Max[r] > 0, n += pop16[[BitAnd[r, 65535] + 1]]; r = BitShiftRight[r, 16]]; n];
bitsOf[m_, k_, al_] := BitAnd[BitShiftRight[k, m - al], 1];          (* occupation of mode al in the states k *)
maskOf[m_, modes_List] := Total[2^(m - Union[modes])];
countIn[m_, k_, modes_List] := If[modes === {}, 0 k, popCount[BitAnd[k, maskOf[m, modes]]]];

(* all m-bit integers with n bits set, in increasing order, built bit by bit (never enumerates 2^m states) *)
combinations[m_, n_] := Module[{A = <|0 -> {0}|>},
   If[n < 0 || n > m, Return[Developer`ToPackedArray[{}, Integer]]];
   Do[A = Association@Table[q -> Join[Lookup[A, q, {}], Lookup[A, q - 1, {}] + 2^(mp - 1)], {q, Max[0, n - (m - mp)], Min[n, mp]}], {mp, m}];
   Developer`ToPackedArray[A[n]]];
(* rank of k among the m-bit integers with the same number of bits set, in increasing order (combinatorial number system):
   r = Sum_i Binomial[p_i, i] over the set bits p_1 < p_2 < ... It is read 12 bits at a time: the bits of a chunk with j set bits
   below it give Sum_i Binomial[p_i, j + i], tabulated for every value of the chunk and every j (rankTables). *)
rankTables[m_] := rankTables[m] = Table[Module[{w = Min[12, m - 12 c], v, bits, cum},
     v = Range[0, 2^w - 1]; bits = Table[BitAnd[BitShiftRight[v, q], 1], {q, 0, w - 1}]; cum = Accumulate[bits];
     Developer`ToPackedArray[Flatten[Transpose[Table[Total[Table[bits[[q + 1]] Binomial[12 c + q, j + cum[[q + 1]]], {q, 0, w - 1}]], {j, 0, m}]]]]],
    {c, 0, Ceiling[m/12] - 1}];
combRank[m_, k_] := Module[{cnt = 0 k, r = 0 k, v, tabs = rankTables[m]},
   Do[v = BitAnd[BitShiftRight[k, 12 (ch - 1)], 4095]; r += tabs[[ch]][[v (m + 1) + cnt + 1]]; cnt += pop16[[v + 1]], {ch, Length[tabs]}]; r];

(* ================================================================= sectors ================================================================ *)
$FockFullLimit = 24;
Options[FockSector] = {"N" -> All, "Sz" -> All, "Counts" -> {}};
FockSector[sp : FockSpace[a_Association], opts : OptionsPattern[]] := Catch[
   With[{key = {sp, OptionValue["N"], OptionValue["Sz"], countsOf[sp, OptionValue["Counts"]]}},
    If[KeyExistsQ[$sectorCache, key], $sectorCache[key],
     With[{r = makeSector @@ key}, If[FockSectorQ[r], $sectorCache[key] = r]; r]]], fwFail];     (* failures are not cached *)
(* "Counts": the modes as ladder operators or indices *)
countsOf[sp_, cnt_] := If[MatchQ[cnt, {{_List, _Integer} ...}], {modeIndices[sp, #[[1]]], #[[2]]} & /@ cnt, cnt];
$sectorCache = <||>;
ClearFockCache[] := ($sectorCache = <||>;);
makeSector[sp : FockSpace[a_Association], n_, sz_, cnt_] := Module[{m = a["Count"], idx, keep, lookup},
   If[! (n === All || IntegerQ[n]), Message[FockSector::opt, "N", n]; Return[$Failed]];
   If[! (sz === All || IntegerQ[2 sz]), Message[FockSector::opt, "Sz", sz]; Return[$Failed]];
   If[! MatchQ[cnt, {{{___Integer}, _Integer} ...}] || ! AllTrue[Flatten[cnt[[All, 1]]], 1 <= # <= m &],
    Message[FockSector::counts, cnt, m]; Return[$Failed]];
   If[sz =!= All && a["Shape"][[3]] =!= 2, Message[FockSector::spin]; Return[$Failed]];
   If[n === All && m > $FockFullLimit, Message[FockSector::big, m, $FockFullLimit]; Return[$Failed]];
   idx = If[n === All, Developer`ToPackedArray[Range[0, 2^m - 1]], combinations[m, n]];
   keep = ConstantArray[1, Length[idx]];
   If[sz =!= All, keep *= 1 - Unitize[countIn[m, idx, Range[1, m, 2]] - countIn[m, idx, Range[2, m, 2]] - 2 sz]];
   Do[keep *= 1 - Unitize[countIn[m, idx, cc[[1]]] - cc[[2]]], {cc, cnt}];
   If[Min[keep] == 0, idx = Pick[idx, keep, 1]];
   (* position look-up: through the rank in the fixed-N sector (array of size Binomial[m, n]) or, without N, an array of size 2^m *)
   lookup = Which[
     n =!= All && Length[idx] == Binomial[m, n], None,
     n =!= All, Module[{map = ConstantArray[0, Binomial[m, n]]}, map[[combRank[m, idx] + 1]] = Range[Length[idx]]; Developer`ToPackedArray[map]],
     True, Module[{map = ConstantArray[0, 2^m]}, map[[idx + 1]] = Range[Length[idx]]; Developer`ToPackedArray[map]]];
   FockSector[<|"Space" -> sp, "Indices" -> Developer`ToPackedArray[idx], "Lookup" -> lookup, "Dimension" -> Length[idx],
     "Constraints" -> <|"N" -> n, "Sz" -> sz, "Counts" -> cnt|>|>]];
FockSector::spin = "\"Sz\" needs exactly two spins.";
FockSector::opt = "Invalid value of \"`1`\": `2`.";
FockSector::counts = "\"Counts\" must be a list {{modes, n}, ...} with modes of the space (`2` modes); got `1`.";
FockSector::big = "A sector without fixed \"N\" on `1` modes needs the whole space; this is allowed up to $FockFullLimit = `2` modes.";
FockSector[a_Association][key_String] := a[key];
FockSectorQ[FockSector[_Association]] = True; FockSectorQ[_] = False;
FockSector /: MakeBoxes[s : FockSector[a_Association], fmt_] :=
  BoxForm`ArrangeSummaryBox[FockSector, s, None, {BoxForm`SummaryItem[{"Dimension: ", a["Dimension"]}],
    BoxForm`SummaryItem[{"Constraints: ", Select[a["Constraints"], # =!= All && # =!= {} &]}]}, {}, fmt];

(* 1-based positions of the states k in the sector, 0 for the states that are not in it *)
positionsIn[s_FockSector, k_] := Module[{m = s["Space"]["Count"], n = s["Constraints"]["N"], lk = s["Lookup"], ok, r},
   If[n === All, lk[[k + 1]],
    ok = 1 - Unitize[popCount[k] - n];
    r = Clip[combRank[m, k], {0, Binomial[m, n] - 1}];
    ok If[lk === None, r + 1, lk[[r + 1]]]]];

toSector[s_FockSector] := s;
toSector[sp_FockSpace] := With[{r = FockSector[sp]}, If[FockSectorQ[r], r, Throw[$Failed, fwFail]]];   (* caught by the public functions *)
spaceOf[s_FockSector] := s["Space"];
spaceOf[sp_FockSpace] := sp;
spaceOf[{s_, _}] := spaceOf[s];

(* sector reached from s by an operator that changes N by dn and Sz by dsz; "Counts" constraints are dropped *)
shiftSector[s_FockSector, dn_, dsz_] := Module[{cs = s["Constraints"]},
   If[(cs["N"] === All || dn == 0) && (cs["Sz"] === All || dsz == 0) && cs["Counts"] === {}, s,
    FockSector[s["Space"], "N" -> If[cs["N"] === All, All, cs["N"] + dn], "Sz" -> If[cs["Sz"] === All, All, cs["Sz"] + dsz]]]];
(* sector reached from s by the string ops: N, Sz and every count are shifted by the exact change the string produces *)
targetSector[s_FockSector, ops_List] := Module[{cs = s["Constraints"], dn = opDeltaN[ops], dsz = opDeltaSz[s["Space"], ops], dc},
   dc = Table[Total[Sign[Select[ops, MemberQ[cc[[1]], Abs[#]] &]]], {cc, cs["Counts"]}];
   If[(cs["N"] === All || dn == 0) && (cs["Sz"] === All || dsz == 0) && AllTrue[dc, # == 0 &], s,
    FockSector[s["Space"], "N" -> If[cs["N"] === All, All, cs["N"] + dn], "Sz" -> If[cs["Sz"] === All, All, cs["Sz"] + dsz],
     "Counts" -> MapThread[{#1[[1]], #1[[2]] + #2} &, {cs["Counts"], dc}]]]];

OccupationBasis[t : (_FockSpace | _FockSector)] := With[{s = toSector[t]}, IntegerDigits[s["Indices"], 2, s["Space"]["Count"]]];
BasisState[t : (_FockSpace | _FockSector), occ_List] := With[{s = toSector[t]},
   If[! (VectorQ[occ, MatchQ[#, 0 | 1] &] && Length[occ] == s["Space"]["Count"]), Message[BasisState::occ, occ, s["Space"]["Count"]]; Return[$Failed]];
   With[{p = positionsIn[s, {FromDigits[occ, 2]}][[1]]}, If[p == 0, Message[BasisState::out]; $Failed, SparseArray[{p} -> 1, {s["Dimension"]}]]]];
BasisState::occ = "`1` is not a list of `2` occupations 0 or 1.";
BasisState::out = "The state is not in the sector.";

(* ================================================================= ladder strings ========================================================= *)
(* apply the operator string ops (right to left) to the states cur0; pos0 is carried along (the source positions).
   Returns {pos kept, image states, sign}. *)
applyString[m_, ops_List, cur0_, pos0_] := Module[{pos = pos0, cur = cur0, sg = ConstantArray[1, Length[cur0]], b, sel, a, par},
   Do[
    a = Abs[o];
    b = BitAnd[BitShiftRight[cur, m - a], 1];
    sel = If[o > 0, 1 - b, b];                                       (* creation needs an empty mode, annihilation an occupied one *)
    If[Min[sel] == 0, {pos, cur, sg} = Pick[#, sel, 1] & /@ {pos, cur, sg}];
    par = BitAnd[popCount[BitShiftRight[cur, m - a + 1]], 1];       (* Jordan-Wigner: modes beta < alpha are the higher bits *)
    sg = sg (1 - 2 par);
    cur = If[o > 0, cur + 2^(m - a), cur - 2^(m - a)],
    {o, Reverse[ops]}];
   {pos, cur, sg}];

opDeltaN[ops_] := Total[Sign[ops]];
opDeltaSz[sp_, ops_] := If[sp["Shape"][[3]] == 2, Total[Sign[#] If[OddQ[Abs[#]], 1/2, -1/2] & /@ ops], 0];
(* a string that creates and annihilates every mode equally often maps each basis state on itself *)
diagonalStringQ[ops_List] := AllTrue[Values[GroupBy[ops, Abs, Total[Sign[#]] &]], # == 0 &];

(* matrix elements of the string ops from the sector "from" to the sector "to": {{row, column} pairs, values} *)
ladderRules[{from_FockSector, to_FockSector}, ops_List] := Module[{m = from["Space"]["Count"], src, img, sg, rows, keep},
   If[from["Dimension"] == 0 || to["Dimension"] == 0, Return[{{}, {}}]];
   {src, img, sg} = applyString[m, ops, from["Indices"], Range[from["Dimension"]]];
   If[src === {}, Return[{{}, {}}]];
   If[from === to && diagonalStringQ[ops], Return[{Transpose[{src, src}], sg}]];          (* no look-up of the images *)
   rows = positionsIn[to, img];
   keep = Unitize[rows];
   If[Total[keep] == 0, {{}, {}}, {Transpose[{Pick[rows, keep, 1], Pick[src, keep, 1]}], Pick[sg, keep, 1]}]];
(* SparseArray that adds repeated entries instead of keeping the first one *)
summedSparse[pos_, vals_, dims_] := If[pos === {}, SparseArray[{}, dims],
   Module[{old = "TreatRepeatedEntries" /. ("SparseArrayOptions" /. SystemOptions["SparseArrayOptions"])},
    Internal`WithLocalSettings[SetSystemOptions["SparseArrayOptions" -> {"TreatRepeatedEntries" -> 1}],
     SparseArray[pos -> vals, dims], SetSystemOptions["SparseArrayOptions" -> {"TreatRepeatedEntries" -> old}]]]];
ladderMatrix[{from_FockSector, to_FockSector}, ops_List] := With[{r = ladderRules[{from, to}, ops]},
   If[r[[1]] === {}, SparseArray[{}, {to["Dimension"], from["Dimension"]}], SparseArray[r[[1]] -> r[[2]], {to["Dimension"], from["Dimension"]}]]];
pairQ[{from : (_FockSpace | _FockSector), to : (_FockSpace | _FockSector)}] := If[spaceOf[from] === spaceOf[to], True,
   Message[FockMatrix::space]; False];
pairQ[_] = False;
FockMatrix::space = "The two sectors belong to different spaces.";
ladder[{from : (_FockSpace | _FockSector), to : (_FockSpace | _FockSector)}, ops_List] :=
  If[pairQ[{from, to}] && checkOps[FockMatrix, spaceOf[from]["Count"], ops], ladderMatrix[{toSector[from], toSector[to]}, ops], $Failed];
ladder[t : (_FockSpace | _FockSector), ops_List] := If[! checkOps[FockMatrix, spaceOf[t]["Count"], ops], $Failed,
   With[{s = toSector[t]}, ladderMatrix[{s, targetSector[s, ops]}, ops]]];

(* fockOperator[target, {{coefficient, string}, ...}]: a sum of strings in one pass, all matrix elements collected and summed in a
   single SparseArray (ladder[target, string] is one string).
   The target is the sector reached by the strings; when they reach different sectors only through "Counts" constraints, the
   target is the (N, Sz) sector without them. *)
termsOK[m_, terms_] := If[MatchQ[terms, {{_, _List} ...}] && AllTrue[terms[[All, 2]], modesOK[m, #] &], True,
   Message[FockMatrix::terms]; False];
FockMatrix::terms = "Terms must be {coefficient, operator string} with nonzero mode indices up to the number of modes.";
fockOperator[{from : (_FockSpace | _FockSector), to : (_FockSpace | _FockSector)}, terms_List] :=
  If[pairQ[{from, to}] && termsOK[spaceOf[from]["Count"], terms], sumOfStrings[toSector[from], toSector[to], terms], $Failed];
fockOperator[t : (_FockSpace | _FockSector), terms_List] := Module[{s, targets, dn, dsz},
   If[! termsOK[spaceOf[t]["Count"], terms], Return[$Failed]];
   s = toSector[t];
   If[terms === {}, Return[SparseArray[{}, {s["Dimension"], s["Dimension"]}]]];
   targets = DeleteDuplicates[targetSector[s, #[[2]]] & /@ terms];
   If[Length[targets] > 1,
    {dn, dsz} = {DeleteDuplicates[opDeltaN /@ terms[[All, 2]]], DeleteDuplicates[opDeltaSz[s["Space"], #] & /@ terms[[All, 2]]]};
    If[Length[dn] > 1 || (s["Constraints"]["Sz"] =!= All && Length[dsz] > 1), Message[FockMatrix::mixed]; Return[$Failed]];
    targets = {shiftSector[s, First[dn], First[dsz]]}];
   sumOfStrings[s, First[targets], terms]];
sumOfStrings[s_FockSector, to_FockSector, terms_] := Module[{parts},
   parts = Select[Table[With[{r = ladderRules[{s, to}, tm[[2]]]}, {r[[1]], tm[[1]] r[[2]]}], {tm, terms}], #[[1]] =!= {} &];
   If[parts === {}, SparseArray[{}, {to["Dimension"], s["Dimension"]}],
    summedSparse[Join @@ parts[[All, 1]], Join @@ parts[[All, 2]], {to["Dimension"], s["Dimension"]}]]];
FockMatrix::mixed = "All the terms must change N (and Sz, on a sector with fixed Sz) by the same amount: give the target {fromSector, toSector}.";

diag[s_FockSector, v_] := With[{d = s["Dimension"], nz = Unitize[v]},       (* Band[] with a long list of values is slow *)
   If[Total[nz] == 0, SparseArray[{}, {d, d}], With[{r = Pick[Range[d], nz, 1]}, SparseArray[Transpose[{r, r}] -> Pick[v, nz, 1], {d, d}]]]];
spinZ[t : (_FockSpace | _FockSector)] := If[! spinCheck[SpinSquared, t], $Failed,
   With[{s = toSector[t], m = spaceOf[t]["Count"]}, diag[s, (countIn[m, s["Indices"], Range[1, m, 2]] - countIn[m, s["Indices"], Range[2, m, 2]])/2]]];
spinPlus[s_FockSector] := fockOperator[{s, shiftSector[s, 0, 1]}, Table[{1, {a, -(a + 1)}}, {a, 1, s["Space"]["Count"], 2}]];
SpinSquared[t : (_FockSpace | _FockSector)] := If[! spinCheck[SpinSquared, t], $Failed,
   With[{s = toSector[t]}, With[{sp = spinPlus[s], sz = spinZ[s]}, ConjugateTranspose[sp] . sp + sz . sz + sz]]];

(* ================================================================= operators of FermiWick ================================================= *)
(* expr -> list of {coefficient, string}: ToLadderBasis writes the composite and defined operators with c and cDag, and iTermList
   gives every product in the order written *)
stringOf[sp_, op_] := Which[
   MatchQ[op, _cDag], {modeIndices[sp, {op}][[1]]},
   MatchQ[op, _c], {-modeIndices[sp, {op}][[1]]},
   True, Message[FockMatrix::op, op]; Throw[$Failed, fwFail]];
fockTerms[sp_, expr_] := Select[{#[[1]], Join @@ (stringOf[sp, #] & /@ #[[2]])} & /@ FermiWick`Private`iTermList[ToLadderBasis[expr]],
   #[[1]] =!= 0 &];
FockMatrix::op = "`1` is not a ladder operator of FermiWick, and ToLadderBasis does not write it with c and cDag.";
FockMatrix[expr_, t : (_FockSpace | _FockSector | {_FockSpace | _FockSector, _FockSpace | _FockSector})] :=
  With[{terms = fockTerms[spaceOf[t], expr]},
   If[terms === {}, With[{s = toSector[If[ListQ[t], First[t], t]]}, SparseArray[{}, {s["Dimension"], s["Dimension"]}]], fockOperator[t, terms]]];
FockExpectation[expr_, v_?VectorQ, t : (_FockSpace | _FockSector)] := With[{s = toSector[t]},
   If[Length[v] != s["Dimension"], Message[FockExpectation::dim, Length[v], s["Dimension"]]; $Failed,
    Conjugate[v] . FockMatrix[expr, {s, s}] . v]];
FockExpectation::dim = "The vector has length `1`, the space or sector dimension `2`.";

(* ================================================================= symmetries ============================================================= *)
ParticleHoleOperator[sp : FockSpace[a_Association], first0_Integer : 1, last_ : Automatic] := Module[{m = a["Count"], first = Max[first0, 1], l, idx, flip, sgn, img},
   l = If[last === Automatic, m, Min[last, m]];
   If[first > l, Return[IdentityMatrix[2^m, SparseArray]]];               (* empty range *)
   idx = Developer`ToPackedArray[Range[0, 2^m - 1]];
   flip = Total[a["Weights"][[first ;; l]]];
   (* exp(-i pi sum_alpha (alpha - 1) n_alpha) on the image, times sigma_x on each mode *)
   img = BitXor[idx, flip];
   sgn = 1 - 2 BitAnd[Total[Table[(al - 1) bitsOf[m, img, al], {al, first, l}]], 1];
   SparseArray[Transpose[{img + 1, idx + 1}] -> sgn, {2^m, 2^m}]];

TimeReversal[sp : FockSpace[a_Association]] := Module[{npair, cup, cdn, g, uloc, U},
   If[a["Shape"][[3]] =!= 2, Message[TimeReversal::spin]; Return[$Failed]];
   npair = a["Shape"][[1]] a["Shape"][[2]];
   cup = KroneckerProduct[{{0, 1}, {0, 0}}, IdentityMatrix[2]];
   cdn = KroneckerProduct[{{1, 0}, {0, -1}}, {{0, 1}, {0, 0}}];
   g = Sum[ConjugateTranspose[{cup, cdn}[[s1]]] . {cup, cdn}[[s2]] PauliMatrix[2][[s1, s2]], {s1, 2}, {s2, 2}];
   uloc = MatrixExp[-I Pi/2 g];
   U = If[npair == 1, SparseArray[uloc], SparseArray[KroneckerProduct @@ ConstantArray[SparseArray[uloc], npair]]];
   <|"U" -> U, "Apply" -> (U . Conjugate[#] &), "Inverse" -> (Conjugate[ConjugateTranspose[U] . #] &)|>];
TimeReversal::spin = "Time reversal needs exactly two spins.";

(* ================================================================= states and kets ======================================================= *)
localSymbols = {"0", "\[DownArrow]", "\[UpArrow]", "\[UpArrow]\[DownArrow]"};
StateToKet[t : (_FockSpace | _FockSector), v_?VectorQ] := Module[{s = toSector[t], sp = spaceOf[t], sh, nz, label},
   If[Length[v] != s["Dimension"], Message[StateToKet::dim, Length[v], s["Dimension"]]; Return[$Failed]];
   sh = sp["Shape"]; nz = Most[ArrayRules[SparseArray[v]]];
   label[k_] := Module[{o = IntegerDigits[s["Indices"][[k]], 2, sp["Count"]]},
     If[sh[[3]] == 2, Partition[localSymbols[[FromDigits[#, 2] + 1]] & /@ Partition[o, 2], sh[[2]]], Partition[o, sh[[2]] sh[[3]]]]];
   Total[(#[[2]] Ket[label[#[[1, 1]]]]) & /@ nz]];
StateToKet::dim = FockExpectation::dim;

(* ================================================================= diagonalisation and Green's functions =================================== *)
SectorBlock[H_, s_FockSector] := H[[s["Indices"] + 1, s["Indices"] + 1]];

LowestStates[H_, k_Integer] := Module[{n = Length[H], vals, vecs, ord},
   If[n == 0, Return[{{}, {}}]];
   If[n <= $FockDenseLimit || k >= n - 1,
    {vals, vecs} = Eigensystem[N@Normal[(H + ConjugateTranspose[H])/2]]; ord = Ordering[Re@vals, Min[k, n]];
    {Re@vals[[ord]], vecs[[ord]]},
    If[k == 1, arnoldi[H, 1, Automatic], lowestComplete[H, k]]]];
(* Arnoldi from a given (or the default) starting vector: the kk lowest Ritz pairs *)
arnoldi[H_, kk_, v0_] := Module[{vals, vecs, ord},
   {vals, vecs} = Eigensystem[-N[H], kk, Method -> Join[{"Arnoldi", "Criteria" -> "RealPart"}, If[v0 === Automatic, {}, {"StartingVector" -> v0}]]];
   ord = Ordering[-Re@vals]; {-Re@vals[[ord]], vecs[[ord]]}];
(* A Krylov space built from one vector contains a single copy of each exactly degenerate eigenvalue, so Arnoldi finds the other
   copies only by accident. For k > 1: run Arnoldi from several random starting vectors, collect the vectors found, and diagonalise H
   in their span (Rayleigh-Ritz), until the k lowest values stop changing. *)
lowestComplete[H_, k_] := Module[{n = Length[H], kk, basis = {}, new, hs, e, v, prev = None, runs = 0, stable = 0, Hn = N[H]},
   kk = Min[k + 5, n - 2];
   While[stable < 2 && runs < 12,
    runs++;
    new = arnoldi[Hn, kk, If[runs == 1, Automatic, RandomReal[{-1, 1}, n]]][[2]];
    basis = Select[Orthogonalize[Join[basis, new]], Norm[#] > 0.5 &];
    hs = Conjugate[basis] . (Hn . Transpose[basis]);
    {e, v} = Eigensystem[(hs + ConjugateTranspose[hs])/2];
    e = Re[e]; With[{o = Ordering[e]}, e = e[[o]]; v = v[[o]]];
    If[prev =!= None && Max[Abs[e[[;; k]] - prev]] < 10^-9 Max[1, Abs[e[[1]]]], stable++, stable = 0];
    prev = e[[;; k]]];
   {e[[;; k]], v[[;; k]] . basis}];

Options[LehmannGreen] = {"Tolerance" -> 10^-9, "ResidueCutoff" -> 10^-12, "N" -> All, "UseSz" -> Automatic};
LehmannGreen[H_, sp_FockSpace, modes0_List, opts : OptionsPattern[]] := Module[{m = sp["Count"], tol = OptionValue["Tolerance"], useSz, keys,
    sec, block, eig, E0, gs, deg, spinOf, poles = {}, res = {}, modes = modeIndices[sp, modes0], isMatrix = MatrixQ[H]},
   If[! posModesOK[LehmannGreen, m, modes], Return[$Failed]];
   (* Sz sectors are used only if H conserves Sz (checked when H is a full-space matrix; a function H must conserve it,
      or be called with "UseSz" -> False) *)
   If[isMatrix && Dimensions[H] =!= {2^m, 2^m}, Message[LehmannGreen::dims, Dimensions[H], 2^m]; Return[$Failed]];
   If[TrueQ[OptionValue["UseSz"]] && sp["Shape"][[3]] =!= 2, Message[LehmannGreen::spin]; Return[$Failed]];
   useSz = Which[OptionValue["UseSz"] =!= Automatic, TrueQ[OptionValue["UseSz"]], sp["Shape"][[3]] =!= 2, False,
     isMatrix, With[{sz = spinZ[sp]}, Max[Abs[H . sz - sz . H]] < 10^-12],
     (* a function H: test the sector with N = M/2 and any Sz *)
     True, With[{s = FockSector[sp, "N" -> Floor[m/2]]}, With[{b = H[s]}, MatrixQ[b] && Dimensions[b] == {s["Dimension"], s["Dimension"]} &&
        With[{sz = spinZ[s]}, Max[Abs[b . sz - sz . b]] < 10^-12]]]];
   If[useSz && isMatrix && OptionValue["UseSz"] === True && With[{sz = spinZ[sp]}, Max[Abs[H . sz - sz . H]] >= 10^-12],
    Message[LehmannGreen::sz]];
   spinOf[a_] := If[useSz, If[OddQ[a], 1/2, -1/2], 0];
   keys = Flatten[Table[If[useSz, Table[{n, sz}, {sz, -Min[n, m - n]/2, Min[n, m - n]/2}], {{n, All}}],
      {n, If[OptionValue["N"] === All, Range[0, m], Flatten[{OptionValue["N"]}]]}], 1];
   sec[k_] := sec[k] = FockSector[sp, "N" -> k[[1]], "Sz" -> k[[2]]];
   block[k_] := With[{s = sec[k]}, With[{b = If[isMatrix, SectorBlock[H, s], H[s]]},
      If[MatrixQ[b, NumericQ] && Dimensions[b] == {s["Dimension"], s["Dimension"]}, b,
       Message[LehmannGreen::block, s["Constraints"]]; Throw[$Failed, fwFail]]]];
   eig[k_] := eig[k] = If[sec[k]["Dimension"] == 0, {{}, {}},
      Module[{b = block[k], vals, vecs}, {vals, vecs} = Eigensystem[N@Normal[(b + ConjugateTranspose[b])/2]]; {Re@vals, vecs}]];
   E0 = Min[Flatten[eig[#][[1]] & /@ keys]];
   gs = Flatten[Table[{k, i}, {k, keys}, {i, Flatten@Position[eig[k][[1]], e_ /; e < E0 + tol]}], 1];
   deg = Length[gs];
   Do[Module[{k = g[[1]], v = eig[g[[1]]][[2, g[[2]]]], targets},
      (* particle (+1) and hole (-1) excitations, grouped by target sector *)
      Do[With[{n = k[[1]] + pm},
         If[0 <= n <= m,
          targets = GatherBy[modes, If[useSz, k[[2]] + pm spinOf[#], All] &];
          Do[Module[{tk = {n, If[useSz, k[[2]] + pm spinOf[First[grp]], All]}, e1, w1, amp},
             If[sec[tk]["Dimension"] > 0,
              {e1, w1} = eig[tk];
              amp = Conjugate[w1] . Transpose[Table[If[MemberQ[grp, a], ladder[{sec[k], sec[tk]}, {pm a}] . v,
                   ConstantArray[0, sec[tk]["Dimension"]]], {a, modes}]];
              poles = Join[poles, pm (e1 - E0)];
              res = Join[res, Table[If[pm == 1, Outer[Times, Conjugate[amp[[q]]], amp[[q]]], Outer[Times, amp[[q]], Conjugate[amp[[q]]]]]/deg,
                 {q, Length[e1]}]]]], {grp, targets}]]], {pm, {1, -1}}]], {g, gs}];
   With[{keep = Flatten@Position[Max[Abs[#]] & /@ res, x_ /; x > OptionValue["ResidueCutoff"]]},
    <|"GroundEnergy" -> E0, "GroundSectors" -> Union[gs[[All, 1]]], "Degeneracy" -> deg, "Modes" -> modes0,
      "Poles" -> poles[[keep]], "Residues" -> res[[keep]]|>]];
LehmannGreen::dims = "H has dimensions `1`; the space has dimension `2`.";
LehmannGreen::spin = "\"UseSz\" -> True needs exactly two spins.";
LehmannGreen::block = "H does not give a numerical square matrix on the sector `1`.";
LehmannGreen::sz = "\"UseSz\" -> True, but H does not conserve Sz: the result is wrong.";
GreenMatrix[g_Association, z_] := Total[MapThread[#2/(z - #1) &, {g["Poles"], g["Residues"]}]];

(* ================================================================= argument guards ======================================================== *)
(* Every public function: a failure deep inside (e.g. a sector that cannot be built, a mode that is not in the space) returns
   $Failed from the outermost call, and a call that matches no definition gives a message and $Failed instead of staying
   unevaluated. *)
General::args = "Invalid arguments in `1`.";
FockSpace[args___] /; ! MatchQ[{args}, {_Association} | {_List, OptionsPattern[]}] := (Message[FockSpace::args, {args}]; $Failed);
$fockPublic = {ModeIndex, OccupationBasis, BasisState, StateToKet, FockMatrix, FockExpectation, SpinSquared, ParticleHoleOperator, TimeReversal, SectorBlock, LowestStates, LehmannGreen, GreenMatrix};
Scan[Function[f,
   f[args___] := (Message[f::args, HoldForm[f[args]]]; $Failed);
   DownValues[f] = Prepend[DownValues[f], HoldPattern[f[args___]] :> Block[{$fwInside = True}, Catch[f[args], fwFail]] /; ! TrueQ[$fwInside]]],
  $fockPublic];

End[];
EndPackage[];
