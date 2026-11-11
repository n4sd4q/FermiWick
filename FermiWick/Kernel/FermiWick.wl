(* ::Package:: *)

(* FermiWick: symbolic algebra of fermion ladder operators, Wick's theorem in declarable Gaussian states,
   mean-field decoupling of any order.

   Conventions
   -----------
   c[i, s] and cDag[i, s] annihilate and create a fermion on site i with spin s ("\[UpArrow]", "\[DownArrow]" or a symbol).
   Sites are any expressions; two different symbolic sites are different sites (KroneckerDelta[i, j] -> 0).
   Products of operators are written with ** (in the given order) or with \[CenterDot] (the normal-ordered product).
   Normal order: all cDag to the left of all c. A Gaussian state is given by its contractions <O1 O2>;
   WickVEV sums over the pairings with the sign of the permutation, MeanField keeps two operators out in every way.

   The normal ordering separates a product site by site: each factor acts on one site and is even (n, Sz, S+-, D, B,
   projectors) or odd (c, cDag), so a product is reordered with the fermion sign only and each site is computed once
   (on a placeholder site, cached). Composite operators stay compact until they have to be expanded. *)

BeginPackage["FermiWick`"];

(* $FermiWickUpdateCheck keeps its setting: the package's own when the package is loaded again, or that of the symbol made
   by typing the name before the package was first loaded (a symbol of Global` or of another context), which is removed
   before the package's symbol is made, so that neither shadows the other. The symbols are named by strings, since a name
   in this code would make the package's symbol at once *)
With[{FermiWick`Private`others = Select[Names["*`$FermiWickUpdateCheck"], StringContainsQ[#, "`"] && ! StringStartsQ[#, "FermiWick`"] &]},
  With[{FermiWick`Private`setting = Join @@ (If[ToExpression[#, InputForm, ValueQ], {ToExpression[#]}, {}] &) /@ FermiWick`Private`others},
   Scan[Quiet[Remove[#]] &, FermiWick`Private`others];
   If[FermiWick`Private`setting =!= {} && ! ToExpression["FermiWick`$FermiWickUpdateCheck", InputForm, ValueQ],
    ToExpression["FermiWick`$FermiWickUpdateCheck = " <> ToString[Last[FermiWick`Private`setting], InputForm]]]]];

$FermiWickVersion::usage = "$FermiWickVersion is the version of FermiWick.";
FermiWickUpdateCheck::usage = "FermiWickUpdateCheck[] asks GitHub (api.github.com) for the latest release of FermiWick and gives \
<|\"InstalledVersion\" -> ..., \"LatestVersion\" -> ..., \"InstallLine\" -> ...|>, the version loaded, the latest release and the \
PacletInstall line of its paclet; when the latest release is newer, it also prints that line. It gives $Failed, with a message, when \
the release cannot be read. Loading FermiWick in a notebook does the same check by itself, at most once a day (see $FermiWickUpdateCheck).";
$FermiWickUpdateCheck::usage = "$FermiWickUpdateCheck is True (the default) when loading FermiWick in a notebook checks for a newer \
release: at most once a day, in the background, it asks GitHub (api.github.com) for the latest release and, if it is newer than \
$FermiWickVersion, adds to the banner the PacletInstall line of that release; any failure is silent, and scripts never check. Set \
$FermiWickUpdateCheck = False before Needs[\"FermiWick`\"] (in init.m to turn it off for good) to load without it. \
FermiWickUpdateCheck[] checks by hand.";

c::usage = "c[l1, l2, ...] annihilates a fermion in the mode with labels l1, l2, ...; the meaning of each label is set by \
DeclareModeLabels. By default c[i] is a spinless fermion on site i and c[i, s] a fermion on site i with spin s (SpinUp = \"\[UpArrow]\", \
SpinDown = \"\[DownArrow]\" or a symbol).";
DeclareModeLabels::usage = "DeclareModeLabels[{type1, type2, ...}] declares the meaning of the labels of the operators c[l1, l2, ...] with that \
number of labels. \"Site\": the site, used for locality (different symbolic sites are different sites); \"Spin\": a spin 1/2 (SpinUp, \
SpinDown or a symbol, split when needed); any other name (\"Orbital\", \"Momentum\", ...): a generic label, compared with KroneckerDelta. \
Defaults: {\"Site\"}, {\"Site\", \"Spin\"}, and {\"Site\", \"Label\", ...} for more labels.";
$ModeLabels::usage = "$ModeLabels is the association number of labels -> types of the labels (set with DeclareModeLabels).";
SpinUp::usage = "SpinUp is the spin \"\[UpArrow]\": c[1, SpinUp] is c[1, \"\[UpArrow]\"].";
SpinDown::usage = "SpinDown is the spin \"\[DownArrow]\".";
$Spins::usage = "$Spins is the list of the spin values, {\"\[UpArrow]\", \"\[DownArrow]\"}: Sum[cDag[i, s] ** c[j, s], {s, $Spins}] is the hopping summed over the spin.";
cDag::usage = "cDag[l1, l2, ...] creates a fermion in the mode with labels l1, l2, ... (see c).";
nOp::usage = "nOp[i, s] is the number of fermions on site i with spin s; nOp[i] is nOp[i, \"\[UpArrow]\"] + nOp[i, \"\[DownArrow]\"].";
nTotOp::usage = "nTotOp[i] is the number of fermions on site i.";
SzOp::usage = "SzOp[i] is the z component of the spin on site i.";
SpOp::usage = "SpOp[i] is the spin raising operator on site i, cDag[i, \[UpArrow]] c[i, \[DownArrow]].";
SmOp::usage = "SmOp[i] is the spin lowering operator on site i, cDag[i, \[DownArrow]] c[i, \[UpArrow]].";
DOp::usage = "DOp[i] is the double occupancy of site i, n[i, \[UpArrow]] n[i, \[DownArrow]]; DOp[i, j, ...] is the product over the sites.";
Bdag::usage = "Bdag[i] creates a doubly occupied site, cDag[i, \[UpArrow]] cDag[i, \[DownArrow]].";
B::usage = "B[i] destroys a doubly occupied site, c[i, \[DownArrow]] c[i, \[UpArrow]].";
SpinOperator::usage = "SpinOperator[i, a] is the component a = 1, 2, 3 of the spin on site i.";
SdotS::usage = "SdotS[i, j] is the operator S_i\[CenterDot]S_j.";
SdotScrossS::usage = "SdotScrossS[i, j, k] is the scalar chirality S_i\[CenterDot](S_j\[Times]S_k).";
EmptyProjector::usage = "EmptyProjector[i] projects site i on the empty state; EmptyProjector[i, j, ...] is the product over the sites.";
SingleOccupationProjector::usage = "SingleOccupationProjector[i] projects site i on the singly occupied states; with several sites it is the product.";

NormalOrder::usage = "NormalOrder[expr] writes expr in normal order (cDag to the left of c). Composite operators are written with c, cDag \
first; a symbolic spin is split into \[UpArrow] and \[DownArrow] with KroneckerDelta.";
Dagger::usage = "Dagger[expr] is the Hermitian conjugate: products reversed, c <-> cDag, SpOp <-> SmOp, B <-> Bdag, coefficients conjugated. \
The amplitudes \[Chi], \[CapitalDelta], \[Rho], \[Kappa] are complex; the other symbols (couplings, angles, phases) are taken as real, unless listed in the option \
\"Complex\" -> {...}. Operators defined with DefineOperator are recognized as Hermitian when they are.";
PlusHC::usage = "PlusHC[expr] is expr + Dagger[expr], expr plus its Hermitian conjugate (\"+ h.c.\"). Options as Dagger.";
HermitianOperatorQ::usage = "HermitianOperatorQ[expr] is True when expr equals its Hermitian conjugate. Options as Dagger.";
OperatorCommutator::usage = "OperatorCommutator[a, b] is the commutator a b - b a, in normal order.";
OperatorAnticommutator::usage = "OperatorAnticommutator[a, b] is the anticommutator a b + b a, in normal order.";
ToLadderBasis::usage = "ToLadderBasis[expr] writes the composite and the defined operators of expr with c and cDag.";
ToSpinChargeBasis::usage = "ToSpinChargeBasis[expr] writes expr with the local operators n, Sz, S+, S-, D, B, Bdag and recognizes the defined \
operators with \"Recognize\" -> True (S_i\[CenterDot]S_j, S_i\[CenterDot](S_j\[Times]S_k), ...). ToSpinChargeBasis[expr, \"Recognize\" -> False] does not \
recognize them.";
GroupLocalOperators::usage = "GroupLocalOperators[expr] normal-orders expr and groups the ladder operators site by site, with the fermion sign.";
OperatorCoefficients::usage = "OperatorCoefficients[expr] is the association monomial -> coefficient of expr in normal order; the constant \
has the key 1. OperatorCoefficients[expr, m] is the coefficient of the monomial m, written in any order (m is normal-ordered first).";
OperatorCoefficients::mono = "`1` is not a single monomial in normal order.";
CleanDeltas::usage ="CleanDeltas[expr] simplifies the KroneckerDelta of sites and spins: two different symbolic sites are different, and \
sums over the spin are done.";
ResetNormalOrder::usage = "ResetNormalOrder[] clears the caches of the normal ordering.";
DefineOperator::usage = "DefineOperator[head, n, definition, display] defines an operator on n sites. definition is a function of the sites \
that writes it with known operators; display is a function of the sites that gives its form on screen. The operator stays compact: \
NormalOrder, WickVEV and MeanField expand it, ToSpinChargeBasis recognizes it again (option \"Recognize\", default True for n >= 2; \
\"Required\" -> All asks the whole operator to be present). Its adjoint stays compact too: Dagger[head[...]] is head itself when Hermitian, head \
with the sites permuted when that is the adjoint (a hopping), otherwise AdjointOperator[head][...]. A one-site operator on several sites is the product: P[i, j] = P[i]**P[j].";
$DefinedOperators::usage = "$DefinedOperators is the association of the operators defined with DefineOperator.";
AdjointOperator::usage = "AdjointOperator[head][sites] is the adjoint of the operator head[sites] defined with DefineOperator, when it is \
neither head itself (Hermitian) nor head with the sites permuted (as a hopping): Dagger[head[i, j]] gives it, it is displayed with a \
dagger, expands like head, and ToSpinChargeBasis recognizes it.";
$SpinSymbols::usage = "$SpinSymbols is the list of the symbols used as spins (updated automatically).";

GaussianState::usage = "GaussianState[symmetries] is the Gaussian state with the given symmetries, written with spin channels: \
<cDag[i, s] c[j, t]> = \[Chi][i, j] \[Delta]_st + Sum_a \[Chi][i, j, a] (\[Sigma]^a)_ts and <c[i, s] c[j, t]> = ((\[CapitalDelta][i, j] + Sum_a \[CapitalDelta][i, j, a] \[Sigma]^a) \[Epsilon])_st, \
a = \"x\", \"y\", \"z\", \[Epsilon]_\[UpArrow]\[DownArrow] = 1. The symmetries set channels to zero or make them real: \"SpinRotation\" (SU(2): only \[Chi], \[CapitalDelta]), \
\"SpinRotationZ\" (rotations about z), \"Charge\" (particle number: no pairing), \"TimeReversal\" (\[Chi], \[CapitalDelta] real; \[Chi]^a, \[CapitalDelta]^a imaginary). \
GaussianState[] has no symmetry. Hermiticity (\[Chi][j, i] = \[Chi][i, j]\[Conjugate]) is always imposed (option \"Hermitian\"). On site: \[Chi][i, i] = <n_i\[Sigma]> averaged \
over \[Sigma] (option \"Filling\"), \[Chi][i, i, a] = <S^a_i> (option \"Magnetization\"), \[CapitalDelta][i, i] = <c[i, \[UpArrow]] c[i, \[DownArrow]]> (option \"OnSitePairing\"). \
GaussianState[{{op1, op2} :> value, ...}] is the state with the contractions given by rules instead. state[\"Symmetries\"], \
state[\"Description\"], state[\"Contraction\"][op1, op2].";
GaussianState::sym = "Unknown symmetry `1`; the symmetries are `2`.";
SpinSingletState::usage = "SpinSingletState[] is GaussianState[{\"SpinRotation\", \"Charge\"}, \"Filling\" -> 1/2]: <cDag[i, s] c[j, t]> = \[Delta]_st \[Chi][i, j]. \
With \"Pairing\" -> \[CapitalDelta] the charge symmetry is dropped: <c[i, s] c[j, t]> = \[Epsilon]_st \[CapitalDelta][i, j] (on site \"OnSitePairing\", default 0). \
Options: \"Hopping\" -> \[Chi], \"Filling\" -> 1/2, \"Pairing\" -> None, \"OnSitePairing\" -> 0, \"Hermitian\" -> True.";
GeneralState::usage = "GeneralState[] is the general Gaussian state written with spin components instead of channels: \
<cDag[i, s] c[j, t]> = \[Rho][i, s, j, t] (Hermitian: \[Rho][j, t, i, s] = \[Rho][i, s, j, t]\[Conjugate]) and, with \"Anomalous\" -> \[Kappa], <c[i, s] c[j, t]> = \
\[Kappa][i, s, j, t] (antisymmetric) and <cDag[i, s] cDag[j, t]> = \[Kappa][j, t, i, s]\[Conjugate]. Options: \"Normal\" -> \[Rho], \"Anomalous\" -> None, \"Hermitian\" -> True.";
ContractionTable::usage = "ContractionTable[state, {i, j, ...}] shows all the contractions between the modes of the sites: \
<cDag_a c_b> and <c_a c_b>. ContractionTable[state, sites, \"Normal\"] and ContractionTable[state, sites, \"Anomalous\"] give the matrices, \
with the modes in the order of Tuples[{sites, $Spins}].";
StateSelect::usage = "StateSelect[state, crit] is the state whose contractions are those of state for the pairs {op1, op2} with crit[op1, op2] \
True, and 0 for the others; e.g. StateSelect[state, ! AnomalousPairQ[##] &]. It is a filter, not a symmetry: it need not give a physical state.";
LocalPairQ::usage = "LocalPairQ[op1, op2] is True when the two ladder operators act on the same site.";
AnomalousPairQ::usage = "AnomalousPairQ[op1, op2] is True when the two ladder operators are both c or both cDag.";
Contraction::usage = "Contraction[state, op1, op2] is the contraction <op1 op2> in state; Contraction[op1, op2] uses $State.";
$State::usage = "$State is the Gaussian state used when none is given (default SpinSingletState[]).";
\[Chi]::usage = "\[Chi][i, j] = (1/2) Sum_s <cDag[i, s] c[j, s]> is the hopping per spin (the charge channel); \[Chi][i, j, a] is the spin channel a = \"x\", \"y\", \"z\", \
and \[Chi][i, i, a] = <S^a_i>.";
\[CapitalDelta]::usage = "\[CapitalDelta][i, j] is the singlet pairing, = <c[i, \[UpArrow]] c[j, \[DownArrow]]> when there is no triplet part, symmetric in the sites; \
\[CapitalDelta][i, j, a] is the triplet pairing (the d vector), antisymmetric in the sites.";
\[Rho]::usage = "\[Rho][i, s, j, t] = <cDag[i, s] c[j, t]> in GeneralState.";
\[Kappa]::usage = "\[Kappa][i, s, j, t] = <c[i, s] c[j, t]> in GeneralState, written with (i, s) before (j, t) in canonical order.";

WickVEV::usage = "WickVEV[expr, state] is the expectation value of expr in the Gaussian state, by Wick's theorem; WickVEV[expr] uses $State. \
Options: \"NormalOrder\" -> False applies Wick without normal ordering first; \"MaxOrder\" -> m keeps only the terms with at most m contractions \
between different sites.";
VerifyNormalOrder::usage = "VerifyNormalOrder[expr, state] gives 0 when WickVEV with and without normal ordering agree, and the difference otherwise.";
MeanField::usage = "MeanField[expr, state] is the mean-field (Hartree\[Dash]Fock\[Dash]Bogoliubov) Hamiltonian of expr in the Gaussian state: every \
normal-ordered product O1...O2n becomes Sum_{p<q} (-1)^(p+q+1) Op Oq <rest> - (n-1) <O1...O2n>, two operators kept out in every way and the \
rest replaced by its expectation value. The result is quadratic, with the same expectation value and the same first variation as expr \
(the variational principle of HFB). MeanField[expr] uses $State.";

ToTeX::usage = "ToTeX[expr] is the LaTeX of an expression of FermiWick: ladder and local operators (c^\\dagger_{1\\uparrow}, S^z_1, D_1, ...), \
spin operators, amplitudes (\[Chi]_{ij}, \[CapitalDelta]_{ij}, \[Rho], \[Kappa]) and spins.";
WickPairings::usage = "WickPairings[product, state] lists the pairings of a product of ladder operators that do not vanish in the state, \
in Wick's theorem applied to the product as written: each is <|\"Sign\" -> \[PlusMinus]1, \"Pairs\" -> {{p, q}, ...}, \"Value\" -> the product of \
the contractions|>. WickPairings[product] uses $State.";
WickTeX::usage = "WickTeX[product, state] is the LaTeX of Wick's theorem for a product of ladder operators: <O1 ... O2n> = the sum of \
its pairings drawn as contractions = their values = the result. The contractions use the LaTeX package simpler-wick \
(\\usepackage{simpler-wick}; compile twice). Options: \"Values\" -> True, \"Result\" -> True. WickTeX[product] uses $State.";
WickTeX::mono = "`1` is not a single product of ladder operators c, cDag (write composite operators with ToLadderBasis first).";
WickTeX::levels = "More than 9 nested contractions: simpler-wick cannot draw them.";

DefinePeierlsHopping::usage ="DefinePeierlsHopping[t, \[CurlyPhi]] makes Subscript[t, a, b] = t E^(I Subscript[\[CurlyPhi], a, b]), with antisymmetric phases \
Subscript[\[CurlyPhi], b, a] = -Subscript[\[CurlyPhi], a, b].";
DefineFlux::usage = "DefineFlux[\[CurlyPhi], {i, j, k, ...}, \[CapitalPhi]] makes the sum of the phases around the oriented loop i\[Rule]j\[Rule]k\[Rule]...\[Rule]i equal to \[CapitalPhi], \
eliminating one free phase (gauge fixing). It can be called for several loops, with common sides.";

DefineOperator::nokey = "No characteristic term of `1` involves all its sites: the operator is expanded but not recognized.";
DefineFlux::dep = "The loop `1` has no free phase: its flux is fixed by the previous loops (sum of the phases = `2`).";
WickVEV::nostate = "`1` is not a GaussianState.";

Begin["`Private`"];

$FermiWickVersion = "0.3.0";
SpinString = {"\[UpArrow]", "\[DownArrow]"};
SpinUp = "\[UpArrow]"; SpinDown = "\[DownArrow]"; $Spins = SpinString;
$ModeLabels = <|1 -> {"Site"}, 2 -> {"Site", "Spin"}|>;
$LabelSymbols = {};
(*True once modes with other labels than {site, spin} appear (or labels are declared): only then the generic paths run*)
$genericModes = False;
iLabelTypes[n_Integer] := Lookup[$ModeLabels, n, Prepend[ConstantArray["Label", n - 1], "Site"]];
iSitePos[n_Integer] := iSitePos[n] = Replace[FirstPosition[iLabelTypes[n], "Site"], {{p_} :> p, _ :> None}];
iGenericQ[n_Integer] := iGenericQ[n] = ! SubsetQ[{"Site", "Spin"}, iLabelTypes[n]];
DeclareModeLabels[types : {__String}] := (
   $ModeLabels[Length[types]] = types; Clear[iSitePos, iGenericQ]; $genericModes = True;
   iSitePos[n_Integer] := iSitePos[n] = Replace[FirstPosition[iLabelTypes[n], "Site"], {{p_} :> p, _ :> None}];
   iGenericQ[n_Integer] := iGenericQ[n] = ! SubsetQ[{"Site", "Spin"}, iLabelTypes[n]];
   ResetNormalOrder[]; types);
$SpinSymbols = {};

(*Patterns of the operators, extended by DefineOperator*)
$ladderPattern = _c | _cDag;
$compositePattern = _nOp | _DOp | _Bdag | _B | _SpOp | _SmOp | _SzOp | _nTotOp;
$opPattern = _c | _cDag | _nOp | _DOp | _Bdag | _B | _SpOp | _SmOp | _SzOp | _nTotOp;
$diagonalPattern = _SzOp | _nTotOp | _nOp | _DOp; (*operators that do not change the quantum numbers*)
$CenterDotChunk = 5000; (*above this estimate of terms, \[CenterDot] normal-orders in pieces*)
$DefinedOperators = <||>;

(* ========================================================================================== *)
(*  1. Products                                                                               *)
(* ========================================================================================== *)

iScalarQ[x_] := FreeQ[x, $opPattern];

(*Expansion of **, done structurally: iTerms[expr] is the list of {coefficient, list of operators}, a product being the
  Cartesian product of the terms of its factors (no pattern matching on the Flat head)*)
iTerms[x_Plus] := Join @@ (iTerms /@ List @@ x);
iTerms[x_NonCommutativeMultiply] :=
  If[AllTrue[List @@ x, MatchQ[#, $opPattern] &], {{1, List @@ x}}, iProductTerms[iTerms /@ List @@ x]];
iTerms[x_Times] :=
  With[{sc = Select[x, iScalarQ], op = Select[x, ! iScalarQ[#] &]},
   Which[
    op === 1, {{x, {}}},
    Head[op] === Times && sc === 1, (Message[iExtractOperators::times, x]; iProductTerms[iTerms /@ List @@ op]),
    True, Map[{sc #[[1]], #[[2]]} &, iTerms[op]]]];
iTerms[Power[b_, k_Integer?Positive]] /; ! iScalarQ[b] := iProductTerms[ConstantArray[iTerms[b], k]];
iTerms[x_?iScalarQ] := {{x, {}}};
iTerms[x_] := {{1, {x}}};
iProductTerms[lists_List] :=
  Fold[Function[{acc, l}, Flatten[Outer[{#1[[1]] #2[[1]], Join[#1[[2]], #2[[2]]]} &, acc, l, 1], 1]], {{1, {}}}, lists];
iNCExpand[expr_] := Expand[Total[Map[#[[1]] iNCM[#[[2]]] &, iTerms[expr]]]];

iNCM[{}] = 1;
iNCM[{x_}] := x;
iNCM[l_List] /; MemberQ[l, 0] := 0;          (*a local factor that vanishes (c c on one mode) makes the product vanish*)
iNCM[l_List] := NonCommutativeMultiply @@ l;

Format[cDag[a__], StandardForm] := Interpretation[Subsuperscript["c", Row[{a}, If[Length[{a}] > 2, ",", ""]], "\[Dagger]"], cDag[a]];
Format[c[a__], StandardForm] := Interpretation[Subscript["c", Row[{a}, If[Length[{a}] > 2, ",", ""]]], c[a]];
Format[nOp[site_, spin_], StandardForm] := Interpretation[Subscript["n", Row[{site, spin}]], nOp[site, spin]];
Format[nOp[site_], StandardForm] := Interpretation[Subscript["n", site], nOp[site]];
Format[DOp[site_], StandardForm] := Interpretation[Subscript["D", site], DOp[site]];
Format[Bdag[site_], StandardForm] := Interpretation[Subsuperscript["B", site, "\[Dagger]"], Bdag[site]];
Format[B[site_], StandardForm] := Interpretation[Subscript["B", site], B[site]];
Format[SpOp[site_], StandardForm] := Interpretation[Subsuperscript["S", site, "+"], SpOp[site]];
Format[SmOp[site_], StandardForm] := Interpretation[Subsuperscript["S", site, "-"], SmOp[site]];
Format[SzOp[site_], StandardForm] := Interpretation[Subsuperscript["S", site, "z"], SzOp[site]];
Format[nTotOp[site_], StandardForm] := Interpretation[Subscript["n", site], nTotOp[site]];

(*In products: brackets around sums and multi-site operators*)
fmtArg[x_Plus] := Row[{"(", x, ")"}];
fmtArg[x_ /; KeyExistsQ[$DefinedOperators, Head[x]] && $DefinedOperators[Head[x]]["Arity"] >= 2] := Row[{"(", x, ")"}];
fmtArg[x_] := x;
Unprotect[NonCommutativeMultiply];
Format[NonCommutativeMultiply[args__]] :=
  With[{disp = Row[fmtArg /@ {args}, ""]}, Interpretation[disp, NonCommutativeMultiply[args]]];
Protect[NonCommutativeMultiply];

(* ========================================================================================== *)
(*  2. Canonical anticommutation relations and normal order                                  *)
(* ========================================================================================== *)

c /: Power[c[___], k_Integer /; k >= 2] := 0;
cDag /: Power[cDag[___], k_Integer /; k >= 2] := 0;

(*Different sites (also symbolic ones) anticommute exactly*)
iAnticommutator[cDag[___], cDag[___]] := 0;
iAnticommutator[c[___], c[___]] := 0;
iAnticommutator[cDag[siteA_, spinA_], c[siteB_, spinB_]] /; iSitePos[2] === 1 && ! iGenericQ[2] := If[siteA === siteB, KroneckerDelta[spinA, spinB], 0];
iAnticommutator[c[siteA_, spinA_], cDag[siteB_, spinB_]] /; iSitePos[2] === 1 && ! iGenericQ[2] := If[siteA === siteB, KroneckerDelta[spinA, spinB], 0];
iAnticommutator[cDag[a___], c[b___]] := iModeDelta[{a}, {b}];
iAnticommutator[c[a___], cDag[b___]] := iModeDelta[{a}, {b}];
(*delta of two modes: sites equal or different (no delta), spins and generic labels by KroneckerDelta*)
iModeDelta[a_List, b_List] := If[Length[a] =!= Length[b], 0, Times @@ MapThread[iLabelDelta, {iLabelTypes[Length[a]], a, b}]];
iLabelDelta["Site", x_, y_] := If[x === y, 1, 0];
iLabelDelta["Spin", x_, y_] := KroneckerDelta[x, y];
iLabelDelta[_, x_, y_] := If[x === y, 1, iRegisterLabels[{x, y}]; KroneckerDelta[x, y]];
iRegisterLabels[l_List] := With[{v = Select[l, ! NumericQ[#] && ! StringQ[#] &]}, If[v =!= {}, $LabelSymbols = Union[$LabelSymbols, v]]];
(*the generic labels of the operators of an expression (only operators with generic labels are looked at)*)
iScanLabels[expr_] :=
  If[$genericModes,
   iRegisterLabels[Flatten[Cases[expr, (c | cDag)[a__] /; iGenericQ[Length[{a}]] :>
       Pick[{a}, ! MemberQ[{"Site", "Spin"}, #] & /@ iLabelTypes[Length[{a}]]], Infinity]]]];

iOpTypeRank[cDag[__]] = 0;
iOpTypeRank[c[__]] = 1;
iOpTypeRank[_] = 2;

iOpLessEqual[op1_, op2_] :=
  With[{t1 = iOpTypeRank[op1], t2 = iOpTypeRank[op2]},
   Which[t1 < t2, True, t1 > t2, False, True, OrderedQ[{op1, op2}]]];

(*\[CenterDot] is the normal-ordered product. The product is normal-ordered as a whole, so that the local operators separate
  site by site without being expanded; only when sums would make too many terms it is normal-ordered in pieces*)
Unprotect[CenterDot];
Clear[CenterDot];
CenterDot[a_] := NormalOrder[a];
CenterDot[args__] :=
  Module[{acc = 1, accSize = 1, sz},
   Scan[
    Function[arg,
     sz = iTermCount[arg];
     If[sz > 1 && accSize > 1 && accSize sz > $CenterDotChunk,
      acc = NormalOrder[acc]; accSize = iTermCount[acc]];
     acc = acc ** arg; accSize = accSize sz],
    {args}];
   NormalOrder[acc]];
Protect[CenterDot];

OperatorCommutator[opA_, opB_] := NormalOrder[opA\[CenterDot]opB - opB\[CenterDot]opA];
OperatorAnticommutator[opA_, opB_] := NormalOrder[opA\[CenterDot]opB + opB\[CenterDot]opA];

iTermCount[x_Plus] := Length[x];
iTermCount[x_Times | x_NonCommutativeMultiply] := Times @@ (iTermCount /@ List @@ x);
iTermCount[_] := 1;

iExtractOperators::times = "Commutative product (Times) of operators in `1`: the canonical order of Mathematica is used. Use ** or \[CenterDot] for a definite order.";
iExtractOperators[term_] :=
  Module[{scalarCoeff, cOps, opList, unpackRule},
   unpackRule = {Power[b : $opPattern, k_Integer /; k >= 2] :> Sequence @@ ConstantArray[b, k]};
   Which[
    Head[term] === NonCommutativeMultiply,
    {1, List @@ term /. unpackRule},
    Head[term] === Times,
    scalarCoeff = Select[term, FreeQ[#, $opPattern] &];
    cOps = Select[term, ! FreeQ[#, $opPattern] &];
    opList = Switch[Head[cOps],
      NonCommutativeMultiply, List @@ cOps,
      Times, (Message[iExtractOperators::times, term];
       Flatten[Replace[List @@ cOps, x_NonCommutativeMultiply :> List @@ x, {1}]]),
      _, {cOps}];
    {scalarCoeff, opList /. unpackRule},
    True, {1, {term} /. unpackRule}]];

(*Expression -> list of {coefficient, ordered list of operators}, one entry per monomial (coefficients expanded, zeros dropped)*)
iTermList[expr_] :=
  DeleteCases[KeyValueMap[{If[Length[#2] == 1, First[#2], Expand[Total[#2]]], #1} &, GroupBy[iTerms[expr], Last -> First]], {0, _}];

(*Applies f to the list of operators of every term, keeping the coefficients*)
iMapOpsTerms[f_, expr_] := Total[Map[If[#[[2]] === {}, #[[1]], #[[1]] f[#[[2]]]] &, iTermList[expr]]];

OperatorCoefficients[expr_] := GroupBy[iTermList[NormalOrder[expr]], (iNCM[Last[#]] &) -> First, Total];
OperatorCoefficients[expr_, mono_] :=
  With[{t = iTermList[NormalOrder[mono]]},
   If[Length[t] == 1, Lookup[OperatorCoefficients[expr], Key[iNCM[t[[1, 2]]]], 0]/t[[1, 1]],
    Message[OperatorCoefficients::mono, mono]; $Failed]];

(*--- the local engine ---*)
iResetLocalCaches[] := (
   Clear[iLocalNO, iLocalTerms, iLocalSC, iDefinedParity];
   iLocalNO[key_List] := iLocalNO[key] = iNormalOrderLadder[iNCM[key]];
   iLocalTerms[key_List] := iLocalTerms[key] = iTermList[iLocalNO[key]];
   iLocalSC[key_List] := iLocalSC[key] = iSCFromLadder[iLocalNO[key]];
   iDefinedParity[h_] := iDefinedParity[h] =
     With[{par = DeleteDuplicates[Mod[Length[#[[2]]], 2] & /@ iTermList[iNonlocalLadder[h, {iSitePH}]]]},
      Which[par === {}, 0, Length[par] == 1, First[par], True, $Failed]]);

(*Local product on site s, computed on the placeholder and moved to s*)
iLocalNOAt[s_, fs_List] := iLocalNO[fs /. s -> iSitePH] /. iSitePH -> s;
iLocalTermsAt[s_, fs_List] := iLocalTerms[fs /. s -> iSitePH] /. iSitePH -> s;
iLocalSCAt[s_, fs_List] := iLocalSC[fs /. s -> iSitePH] /. iSitePH -> s;

(*{site, parity} of a factor, or $Failed if it is not a local operator of definite parity*)
$parityOf = <|c -> 1, cDag -> 1, nOp -> 0, nTotOp -> 0, SzOp -> 0, SpOp -> 0, SmOp -> 0, DOp -> 0, Bdag -> 0, B -> 0|>;
iFactorSiteParity[f_] :=
  With[{p = Lookup[$parityOf, Head[f], None]},
   Which[
    p === 1 && Length[f] === 2 && ! $genericModes, {First[f], 1},
    p === 1, ($genericModes = True; With[{q = iSitePos[Length[f]]}, If[q === None, $Failed, {f[[q]], 1}]]),
    p =!= None, {First[f], p},
    KeyExistsQ[$DefinedOperators, Head[f]] && $DefinedOperators[Head[f]]["Arity"] === 1 && Length[f] === 1,
    With[{par = iDefinedParity[Head[f]]}, If[par === $Failed, $Failed, {First[f], par}]],
    True, $Failed]];

(*Reorders the factors by site: {sign, {{site, factors on that site}, ...}} or $Failed*)
iSiteFactorize[factors_List] :=
  Module[{info = iFactorSiteParity /@ factors, sites, perm, newPos, sign},
   If[MemberQ[info, $Failed], $Failed,
    sites = info[[All, 1]];
    perm = Ordering[Transpose[{sites, Range[Length[factors]]}]];
    newPos = Ordering[perm];
    sign = Signature[newPos[[Flatten[Position[info[[All, 2]], 1]]]]];
    {sign, Map[{sites[[First[#]]], factors[[#]]} &, SplitBy[perm, sites[[#]] &]]}]];

(*Combines the terms of local normal-ordered polynomials on different sites: no contractions, only reordering and sign.
  Result: list of {coefficient, operators}*)
iCombineGlobal[termLists_List] :=
  If[MemberQ[termLists, {}], {},
   Map[
    Function[combo,
     With[{ops = Join @@ combo[[All, 2]]},
      With[{perm = Ordering[Transpose[{iOpTypeRank /@ ops, ops}]]},
       {(Times @@ combo[[All, 1]]) Signature[perm], ops[[perm]]}]]],
    Tuples[termLists]]];

(*A list of terms -> the expanded sum*)
iFromTerms[terms_List] :=
  Total[Flatten[Map[With[{op = iNCM[#[[2]]]}, If[Head[#[[1]]] === Plus, op List @@ #[[1]], {#[[1]] op}]] &, terms]]];

(*Merges the terms with the same operators; coefficients expanded, deltas cleaned, zeros dropped*)
iCleanCoef[x_] := With[{y = Expand[x]}, If[FreeQ[y, KroneckerDelta], y, Expand[CleanDeltas[y]]]];
iMergeTerms[terms_List] :=
  DeleteCases[KeyValueMap[{iCleanCoef[If[Length[#2] == 1, First[#2], Total[#2]]], #1} &, GroupBy[terms, Last -> First]], {0, _}];

(*Plain normal ordering (everything in c, cDag): used for the local products and as a fallback*)
iNormalOrderLadder[expr_] :=
  CleanDeltas[Expand[iMapOpsTerms[iNormalOrderList,
     If[FreeQ[expr, $compositePattern], expr, ToLadderBasis[expr]]]]];

iSortedLadderQ[fs_List] :=
  AllTrue[fs, MatchQ[#, _c | _cDag] &] && DuplicateFreeQ[fs] && OrderedQ[Transpose[{iOpTypeRank /@ fs, fs}]];

(*One term {coef, factors} -> its normal-ordered terms*)
iNOTerm[{coef_, factors_List}] :=
  Which[
   factors === {} || iSortedLadderQ[factors], {{coef, factors}},
   True,
   With[{fz = iSiteFactorize[factors]},
    If[fz === $Failed,
     Map[{coef #[[1]], #[[2]]} &, iTermList[iNormalOrderLadder[iNCM[factors]]]],
     With[{k = coef fz[[1]]}, Map[{k #[[1]], #[[2]]} &, iCombineGlobal[iLocalTermsAt @@@ fz[[2]]]]]]]];

(*Normal-ordered terms of an expression without symbolic spins and multi-site operators*)
iNormalOrderTerms[e_] := iMergeTerms[Join @@ (iNOTerm /@ iTermList[e])];

(*Terms of expr in normal order (symbolic spins: through NormalOrder)*)
iNormalOrderedTermList[expr_] :=
  With[{e = iExpandMultiSite[expr]},
   If[iSpinVariables[e] === {}, iNormalOrderTerms[e], iTermList[NormalOrder[e]]]];

(*Operators defined on several sites (S.S, ...) are expanded before separating by site*)
iMultiSitePattern[] :=
  Alternatives @@ (Blank /@ Select[Keys[$DefinedOperators], $DefinedOperators[#]["Arity"] >= 2 &]);
iExpandMultiSite[expr_] :=
  With[{pat = iMultiSitePattern[]},
   If[FreeQ[expr, pat], expr, expr /. (op : pat) :> iNonlocalLadder[Head[op], List @@ op]]];

(*The cache of iNormalOrderList depends only on the anticommutation relations. No Return[] inside the memoization*)
ResetNormalOrder[] := (
   Clear[iNormalOrderList];
   iNormalOrderList[{}] = 1;
   iNormalOrderList[{op_}] := op;
   iNormalOrderList[opsList_List] := iNormalOrderList[opsList] =
     Module[{nLen = Length[opsList], pos = None, op1, op2, before, after, fermPair, acVal, idx},
      Do[
       If[(opsList[[idx]] === opsList[[idx + 1]] && MatchQ[opsList[[idx]], $ladderPattern]) ||
          ! iOpLessEqual[opsList[[idx]], opsList[[idx + 1]]],
        pos = idx; Break[]],
       {idx, 1, nLen - 1}];
      Which[
       pos === None, iNCM[opsList],
       opsList[[pos]] === opsList[[pos + 1]] && MatchQ[opsList[[pos]], $ladderPattern], 0,
       True,
       op1 = opsList[[pos]]; op2 = opsList[[pos + 1]];
       before = Take[opsList, pos - 1]; after = Drop[opsList, pos + 1];
       fermPair = MatchQ[op1, $ladderPattern] && MatchQ[op2, $ladderPattern];
       acVal = If[fermPair, iAnticommutator[op1, op2], 0];
       Expand[
        If[acVal === 0, 0, acVal iNormalOrderList[Join[before, after]]] +
         If[fermPair, -1, 1] iNormalOrderList[Join[before, {op2, op1}, after]]]]];
   iResetLocalCaches[];
   iResetNonlocalCache[]);

(*Finds the symbolic spin variables and records them in $SpinSymbols*)
iSpinVariables[expr_] :=
  Module[{vars},

   vars = DeleteDuplicates@If[$genericModes,
      Select[Flatten[Join[Cases[expr, (c | cDag)[a__] :> Pick[{a}, # === "Spin" & /@ iLabelTypes[Length[{a}]]], Infinity],
         Cases[expr, nOp[_, s_] :> s, Infinity]]], (! MemberQ[SpinString, #] && ! NumericQ[#]) &],
      Cases[expr, (c | cDag | nOp)[_, s_] /; (! MemberQ[SpinString, s] && ! NumericQ[s]) :> s, Infinity]];
   If[vars =!= {}, $SpinSymbols = Union[$SpinSymbols, vars]];
   vars];

NormalOrder[expr_] :=
  Module[{e = iExpandMultiSite[expr], spinVars, sVar},
   iScanLabels[e];
   spinVars = iSpinVariables[e];
   If[spinVars =!= {},
    sVar = First[spinVars];
    CleanDeltas[
     KroneckerDelta[sVar, "\[UpArrow]"] NormalOrder[e /. sVar -> "\[UpArrow]"] +
      KroneckerDelta[sVar, "\[DownArrow]"] NormalOrder[e /. sVar -> "\[DownArrow]"]],
    iFromTerms[iNormalOrderTerms[e]]]];

(*--- Hermitian conjugate ---*)
$complexHeads = {\[Chi], \[CapitalDelta], \[Rho], \[Kappa]};
Options[Dagger] = {"Complex" -> {}};
(*Complex conjugate of a scalar: numbers, the amplitudes and the symbols declared complex; the other symbols are real*)
iConjScalar[e_, cplx_] :=
  e /. {Conjugate[x_] :> x, Complex[a_, b_] :> Complex[a, -b],
    x : (h_Symbol[___] /; MemberQ[$complexHeads, h] || MemberQ[cplx, h]) :> Conjugate[x],
    x_Symbol /; MemberQ[cplx, x] :> Conjugate[x]};
iDag[e_Plus, cplx_] := iDag[#, cplx] & /@ e;
iDag[e_NonCommutativeMultiply, cplx_] := NonCommutativeMultiply @@ Reverse[iDag[#, cplx] & /@ List @@ e];
iDag[e_Times, cplx_] := With[{sc = Select[e, iScalarQ], op = Select[e, ! iScalarQ[#] &]}, iConjScalar[sc, cplx] iDag[op, cplx]];
iDag[Power[b_, k_Integer?Positive], cplx_] /; ! iScalarQ[b] := NonCommutativeMultiply @@ ConstantArray[iDag[b, cplx], k];
iDag[c[a___], _] := cDag[a];
iDag[cDag[a___], _] := c[a];
iDag[SpOp[i_], _] := SmOp[i];
iDag[SmOp[i_], _] := SpOp[i];
iDag[B[i_], _] := Bdag[i];
iDag[Bdag[i_], _] := B[i];
iDag[x : (_nOp | _nTotOp | _SzOp | _DOp), _] := x;
(*A defined operator stays compact: itself when Hermitian, itself with the sites permuted when that is its adjoint (a hopping),
  otherwise its adjoint operator, registered by DefineOperator*)
iDag[x_, cplx_] /; KeyExistsQ[$DefinedOperators, Head[x]] :=
  With[{a = Lookup[$DefinedOperators[Head[x]], "Adjoint", None]},
   Switch[a,
    "Hermitian", x,
    {"Permutation", _}, Head[x] @@ (List @@ x)[[a[[2]]]],
    {"Head", _}, a[[2]] @@ x,
    _, If[iDefinedHermitianQ[Head[x]], x, iDag[NormalOrder[$DefinedOperators[Head[x]]["Definition"] @@ x], cplx]]]];
iDag[x_, cplx_] := iConjScalar[x, cplx];
(*A defined operator is Hermitian when its definition is (checked once, on placeholder sites)*)
iDefinedHermitianQ[h_] := iDefinedHermitianQ[h] =
   With[{d = $DefinedOperators[h]["Definition"] @@ $DefinedOperators[h]["Sites"]},
    AllTrue[Values[OperatorCoefficients[iDag[NormalOrder[d], {}] - d]], Simplify[#] === 0 &]];
Dagger[expr_, OptionsPattern[]] := iDag[expr, OptionValue["Complex"]];
PlusHC[expr_, opts : OptionsPattern[Dagger]] := expr + Dagger[expr, opts];
HermitianOperatorQ[expr_, opts : OptionsPattern[Dagger]] := AllTrue[Values[OperatorCoefficients[expr - Dagger[expr, opts]]], Simplify[#] === 0 &];

(* ========================================================================================== *)
(*  3. Deltas                                                                                 *)
(* ========================================================================================== *)

(*Sums over the spin (for any common coefficient x)*)
$spinSumRules = With[{u = "\[UpArrow]", d = "\[DownArrow]"}, {
    KroneckerDelta[u, s_] KroneckerDelta[d, s_] :> 0,
    x_. KroneckerDelta[u, s1_] KroneckerDelta[u, s2_] + x_. KroneckerDelta[d, s1_] KroneckerDelta[d, s2_] :> x KroneckerDelta[s1, s2],
    x_. KroneckerDelta[u, s1_] KroneckerDelta[d, s2_] + x_. KroneckerDelta[d, s1_] KroneckerDelta[u, s2_] :> x - x KroneckerDelta[s1, s2],
    x_. KroneckerDelta[u, s_] + x_. KroneckerDelta[d, s_] :> x}];

(*KroneckerDelta is Orderless, so one version of each rule is enough. extra: rules of the state, applied first*)
CleanDeltas[expr_] := CleanDeltas[expr, {}];
CleanDeltas[expr_, extra_List] :=
  Module[{spinSymbols, isSpin, rules},
   spinSymbols = DeleteDuplicates@Join[
      Cases[expr, (c | cDag | nOp)[_, s_] :> s, Infinity],
      Cases[expr, KroneckerDelta[s_, "\[UpArrow]" | "\[DownArrow]"] :> s, Infinity],
      Cases[expr, s_Symbol /; StringContainsQ[ToString[s], "\[Sigma]"], Infinity],
      $SpinSymbols];
   isSpin[x_] := MemberQ[spinSymbols, x] || MemberQ[SpinString, x] || MemberQ[$LabelSymbols, x];
   rules = Join[extra, {
      KroneckerDelta[a_, b_] /; (a =!= b && ! isSpin[a] && ! isSpin[b]) :> 0,
      KroneckerDelta[a_, b_]^k_Integer?Positive :> KroneckerDelta[a, b]}];
   (*a delta of generic labels substitutes one label for the other in what it multiplies (only when there are generic labels)*)
   If[$LabelSymbols =!= {} && ! FreeQ[expr, KroneckerDelta[_, Alternatives @@ $LabelSymbols]],
    AppendTo[rules, (dd : KroneckerDelta[a_, b_]) x_ /; OrderedQ[{a, b}] && MemberQ[$LabelSymbols, b] && ! FreeQ[x, b] :> dd (x /. b -> a)]];
   If[! FreeQ[expr, KroneckerDelta[_, "\[UpArrow]" | "\[DownArrow]"]], rules = Join[rules, $spinSumRules]];
   expr //. rules];

(* ========================================================================================== *)
(*  4. Gaussian states                                                                         *)
(* ========================================================================================== *)

GaussianState[a_Association][key_String] := Lookup[a, key, Missing["KeyAbsent", key]];
Format[GaussianState[a_Association], StandardForm] :=
  Interpretation[Row[{"GaussianState[", Style[Lookup[a, "Description", "\[Ellipsis]"], Italic], "]"}], GaussianState[a]];

LocalPairQ[x_, y_] := iSiteOf[x] === iSiteOf[y] && ! MissingQ[iSiteOf[x]];
iSiteOf[x : (_c | _cDag)] := With[{q = iSitePos[Length[x]]}, If[q === None, Missing["NoSite"], x[[q]]]];
AnomalousPairQ[x_, y_] := Head[x] === Head[y];

(*Epsilon of the singlet, \[Epsilon]_\[UpArrow]\[DownArrow] = 1*)
iEps[s1_, s2_] := KroneckerDelta[s1, "\[UpArrow]"] KroneckerDelta[s2, "\[DownArrow]"] - KroneckerDelta[s1, "\[DownArrow]"] KroneckerDelta[s2, "\[UpArrow]"];

(*--- from rules ---*)
Options[GaussianState] = {"Rules" -> {}, "Description" -> Automatic, "Hopping" -> \[Chi], "Pairing" -> \[CapitalDelta],
   "Filling" -> Automatic, "Magnetization" -> Automatic, "OnSitePairing" -> Automatic, "Hermitian" -> True};
GaussianState[rules : {(_Rule | _RuleDelayed) ..}, opts : OptionsPattern[]] :=
  GaussianState[<|"Contraction" -> iRuleContraction[rules], "Rules" -> OptionValue[GaussianState, {opts}, "Rules"],
    "Description" -> Replace[OptionValue[GaussianState, {opts}, "Description"], Automatic -> "user-defined"],
    "Symmetries" -> Missing["Unknown"], "Amplitudes" -> {}|>];
iRuleContraction[rules_] := Function[{x, y}, iFromRules[rules, x, y]];
iFromRules[rules_, x_, y_] :=
  Module[{r = Replace[{x, y}, Append[rules, _ -> iNoRule]], n},
   Which[
    r =!= iNoRule, r,
    MatchQ[{x, y}, {_c, _cDag}],
    n = Replace[{y, x}, Append[rules, _ -> iNoRule]];
    If[n === iNoRule, 0, iAnticommutator[x, y] - n],
    True, 0]];

(*--- from symmetries: spin channels ---
  <cDag_{i s} c_{j t}> = G_ij[[s, t]],  G_ij = \[Chi]_ij 1 + Sum_a \[Chi]^a_ij Transpose[\[Sigma]^a]   (\[Chi]^a_ii = <S^a_i>)
  <c_{i s} c_{j t}>    = F_ij[[s, t]],  F_ij = (\[CapitalDelta]_ij 1 + Sum_a \[CapitalDelta]^a_ij \[Sigma]^a) . \[Epsilon]    (\[CapitalDelta]_ji = \[CapitalDelta]_ij, \[CapitalDelta]^a_ji = -\[CapitalDelta]^a_ij)
  Hermiticity: G_ji = G_ij^dagger, i.e. \[Chi]_ji = \[Chi]_ij^*, \[Chi]^a_ji = \[Chi]^a_ij^*. A symmetry is a set of channels that vanish, plus reality
  conditions: SU(2) leaves only the scalars; rotations about z leave the z components; U(1) of charge kills F; time reversal
  T = i\[Sigma]^y K gives G = \[Epsilon] G^* \[Epsilon]^T and F = \[Epsilon] F^* \[Epsilon]^T, i.e. \[Chi], \[CapitalDelta] real and \[Chi]^a, \[CapitalDelta]^a imaginary.*)
$symmetries = <|
   "SpinRotation" -> {{"hop", "x"}, {"hop", "y"}, {"hop", "z"}, {"pair", "x"}, {"pair", "y"}, {"pair", "z"}},
   "SpinRotationZ" -> {{"hop", "x"}, {"hop", "y"}, {"pair", "x"}, {"pair", "y"}},
   "Charge" -> {{"pair", 0}, {"pair", "x"}, {"pair", "y"}, {"pair", "z"}},
   "TimeReversal" -> {}|>;
$symmetryNames = <|"SpinRotation" -> "spin rotations", "SpinRotationZ" -> "spin rotations about z", "Charge" -> "charge",
   "TimeReversal" -> "time reversal"|>;
$pauli = <|"x" -> {{0, 1}, {1, 0}}, "y" -> {{0, -I}, {I, 0}}, "z" -> {{1, 0}, {0, -1}}|>;
$eps = {{0, 1}, {-1, 0}};
$spinIdx = <|"\[UpArrow]" -> 1, "\[DownArrow]" -> 2|>;

GaussianState[] := GaussianState[{}];
GaussianState[sym_String, opts : OptionsPattern[]] := GaussianState[{sym}, opts];
GaussianState[syms : {___String}, opts : OptionsPattern[]] :=
  Module[{bad = Complement[syms, Keys[$symmetries]], cfg, h, d, rules},
   If[bad =!= {}, Message[GaussianState::sym, First[bad], Keys[$symmetries]]; Return[$Failed]];
   h = OptionValue[GaussianState, {opts}, "Hopping"]; d = OptionValue[GaussianState, {opts}, "Pairing"];
   cfg = <|"Hopping" -> h, "Pairing" -> d, "Zero" -> Union @@ Lookup[$symmetries, syms, {}],
     "Hermitian" -> TrueQ[OptionValue[GaussianState, {opts}, "Hermitian"]],
     "Filling" -> OptionValue[GaussianState, {opts}, "Filling"],
     "Magnetization" -> OptionValue[GaussianState, {opts}, "Magnetization"],
     "OnSitePairing" -> OptionValue[GaussianState, {opts}, "OnSitePairing"]|>;
   rules = Join[
     If[cfg["Hermitian"], {Conjugate[h[a_, a_]] :> h[a, a], Conjugate[h[a_, a_, b_String]] :> h[a, a, b]}, {}],
     If[MemberQ[syms, "TimeReversal"],
      {h[a_, a_, _String] :> 0, Conjugate[h[a_, b_]] :> h[a, b], Conjugate[h[a_, b_, k_String]] :> -h[a, b, k],
       Conjugate[d[a_, b_]] :> d[a, b], Conjugate[d[a_, b_, k_String]] :> -d[a, b, k]}, {}],
     If[cfg["Filling"] =!= Automatic, {h[a_, a_] :> cfg["Filling"], KroneckerDelta[a_, b_] h[a_, b_] :> KroneckerDelta[a, b] cfg["Filling"]}, {}],
     OptionValue[GaussianState, {opts}, "Rules"]];
   GaussianState[<|
     "Contraction" -> If[MemberQ[syms, "SpinRotation"],
       With[{hh = h, nu = cfg["Filling"], herm = cfg["Hermitian"], dd = If[MemberQ[syms, "Charge"], None, d], d0 = cfg["OnSitePairing"]},
        Function[{x, y}, iSinglet[hh, nu, herm, dd, d0, x, y]]],
       With[{c0 = cfg}, Function[{x, y}, iChannelContraction[c0, x, y]]]],
     "Rules" -> rules, "Symmetries" -> Sort[syms],
     "Description" -> Replace[OptionValue[GaussianState, {opts}, "Description"],
       Automatic -> If[syms === {}, "no symmetry", "symmetric under " <> StringRiffle[Lookup[$symmetryNames, Sort[syms]], ", "]]],
     "Hopping" -> h, "Amplitudes" -> If[MemberQ[syms, "Charge"], {h}, {h, d}]|>]];

(*Amplitude of a channel {kind, component} on the ordered pair of sites (i, j)*)
iAmp[cfg_, kc : {k_, a_}, i_, j_] := If[MemberQ[cfg["Zero"], kc], 0, If[k === "hop", iHopAmp[cfg, a, i, j], iPairAmp[cfg, a, i, j]]];
iHopSym[h_, 0, i_, j_] := h[i, j];
iHopSym[h_, a_String, i_, j_] := h[i, j, a];
iHopAmp[cfg_, a_, i_, j_] :=
  Which[
   i === j, If[a === 0, Replace[cfg["Filling"], Automatic :> cfg["Hopping"][i, i]], iMagnetization[cfg, a, i]],
   ! cfg["Hermitian"] || OrderedQ[{i, j}], iHopSym[cfg["Hopping"], a, i, j],
   True, Conjugate[iHopSym[cfg["Hopping"], a, j, i]]];
iMagnetization[cfg_, a_, i_] :=
  With[{m = cfg["Magnetization"]},
   Which[
    m === Automatic, cfg["Hopping"][i, i, a],
    Head[m] === Function, iMagComponent[m[i], a],
    True, iMagComponent[m, a]]];
iMagComponent[v_List, a_] := v[[Lookup[<|"x" -> 1, "y" -> 2, "z" -> 3|>, a]]];
iMagComponent[v_, a_] := If[a === "z", v, 0];
iPairAmp[cfg_, a_, i_, j_] :=
  Which[
   i === j, If[a === 0, Replace[cfg["OnSitePairing"], Automatic :> cfg["Pairing"][i, i]], 0],
   a === 0, cfg["Pairing"] @@ Sort[{i, j}],
   OrderedQ[{i, j}], cfg["Pairing"][i, j, a],
   True, -cfg["Pairing"][j, i, a]];

iHopMatrix[cfg_, i_, j_] :=
  iAmp[cfg, {"hop", 0}, i, j] IdentityMatrix[2] + Sum[iAmp[cfg, {"hop", a}, i, j] Transpose[$pauli[a]], {a, {"x", "y", "z"}}];
iPairMatrix[cfg_, i_, j_] :=
  (iAmp[cfg, {"pair", 0}, i, j] IdentityMatrix[2] + Sum[iAmp[cfg, {"pair", a}, i, j] $pauli[a], {a, {"x", "y", "z"}}]) . $eps;
(*Matrix element for explicit spins, or written with KroneckerDelta for symbolic ones*)
iElem[m_, s_, t_] :=
  If[KeyExistsQ[$spinIdx, s] && KeyExistsQ[$spinIdx, t], m[[$spinIdx[s], $spinIdx[t]]],
   Sum[KroneckerDelta[s, SpinString[[u]]] KroneckerDelta[t, SpinString[[v]]] m[[u, v]], {u, 2}, {v, 2}]];
(*Complex conjugate, term by term (Conjugate does not distribute over sums and symbolic products)*)
iConjExpr[e_] := Total[Map[If[Head[#] === Times, Conjugate /@ #, Conjugate[#]] &, If[Head[#] === Plus, List @@ #, {#}] &[Expand[e]]]];

iChannelContraction[cfg_, cDag[i_, s_], c[j_, t_]] := iElem[iHopMatrix[cfg, i, j], s, t];
iChannelContraction[cfg_, c[i_, s_], cDag[j_, t_]] := If[i === j, KroneckerDelta[s, t], 0] - iElem[iHopMatrix[cfg, j, i], t, s];
iChannelContraction[cfg_, c[i_, s_], c[j_, t_]] := iElem[iPairMatrix[cfg, i, j], s, t];
iChannelContraction[cfg_, cDag[i_, s_], cDag[j_, t_]] := iElem[Map[iConjExpr, iPairMatrix[cfg, j, i], {2}], t, s];
iChannelContraction[___] := 0;

(*Fast path for states invariant under spin rotations: G = \[Chi] 1, F = \[CapitalDelta] \[Epsilon], parameters fixed when the state is made*)
iHop[h_, nu_, herm_, i_, j_] := Which[i === j, If[nu === Automatic, h[i, i], nu], ! herm || OrderedQ[{i, j}], h[i, j], True, Conjugate[h[j, i]]];
iPair[d_, d0_, i_, j_] := If[i === j, If[d0 === Automatic, d[i, i], d0], d @@ Sort[{i, j}]];
iSinglet[h_, nu_, herm_, _, _, cDag[i_, s_], c[j_, t_]] := iHop[h, nu, herm, i, j] KroneckerDelta[s, t];
iSinglet[h_, nu_, herm_, _, _, c[i_, s_], cDag[j_, t_]] := If[i === j, 1 - iHop[h, nu, herm, i, i], -iHop[h, nu, herm, j, i]] KroneckerDelta[s, t];
iSinglet[_, _, _, d_, d0_, c[i_, s_], c[j_, t_]] /; d =!= None := iEps[s, t] iPair[d, d0, i, j];
iSinglet[_, _, _, d_, d0_, cDag[i_, s_], cDag[j_, t_]] /; d =!= None := iEps[t, s] Conjugate[iPair[d, d0, j, i]];
iSinglet[___] := 0;

(*--- invariant under spin rotations, the usual shortcut ---*)
Options[SpinSingletState] = {"Hopping" -> \[Chi], "Filling" -> 1/2, "Pairing" -> None, "OnSitePairing" -> 0, "Hermitian" -> True};
SpinSingletState[OptionsPattern[]] :=
  With[{d = OptionValue["Pairing"]},
   GaussianState[If[d === None, {"SpinRotation", "Charge"}, {"SpinRotation"}],
    "Hopping" -> OptionValue["Hopping"], "Pairing" -> Replace[d, None -> \[CapitalDelta]], "Filling" -> OptionValue["Filling"],
    "OnSitePairing" -> OptionValue["OnSitePairing"], "Hermitian" -> OptionValue["Hermitian"],
    "Description" -> If[d === None, "spin singlet", "spin singlet with pairing"]]];

(*--- general, by spin components ---*)
Options[GeneralState] = {"Normal" -> \[Rho], "Anomalous" -> None, "Hermitian" -> True};
GeneralState[OptionsPattern[]] :=
  With[{r = OptionValue["Normal"], k = OptionValue["Anomalous"], herm = TrueQ[OptionValue["Hermitian"]]},
   GaussianState[<|"Contraction" -> Function[{x, y}, iGeneral[r, k, herm, x, y]],
     "Rules" -> If[herm, {Conjugate[r[x__]] /; EvenQ[Length[{x}]] && Take[{x}, Length[{x}]/2] === Drop[{x}, Length[{x}]/2] :> r[x]}, {}],
     "Symmetries" -> {},
     "Description" -> If[k === None, "general, normal", "general, normal and anomalous"],
     "Amplitudes" -> DeleteCases[{r, k}, None]|>]];
iRho[r_, herm_, x_, y_] := If[! herm || OrderedQ[{x, y}], r @@ Join[x, y], Conjugate[r @@ Join[y, x]]];
iGeneral[r_, _, herm_, cDag[a___], c[b___]] := iRho[r, herm, {a}, {b}];
iGeneral[r_, _, herm_, c[a___], cDag[b___]] := iModeDelta[{a}, {b}] - iRho[r, herm, {b}, {a}];
iGeneral[_, k_, _, c[a___], c[b___]] /; k =!= None := iKappa[k, {a}, {b}];
iGeneral[_, k_, _, cDag[a___], cDag[b___]] /; k =!= None := Conjugate[iKappa[k, {b}, {a}]];
iGeneral[___] := 0;
iKappa[k_, x_, y_] := Which[x === y, 0, OrderedQ[{x, y}], k @@ Join[x, y], True, -k @@ Join[y, x]];

(*--- selection ---*)
StateSelect[GaussianState[a_Association], crit_] :=
  GaussianState[Join[a, <|
    "Contraction" -> With[{f = a["Contraction"]}, Function[{x, y}, If[TrueQ[crit[x, y]], f[x, y], 0]]],
    "Description" -> Lookup[a, "Description", ""] <> ", selected"|>]];

$State = SpinSingletState[];

Contraction[x_, y_] := Contraction[$State, x, y];
Contraction[GaussianState[a_Association], x_, y_] := CleanDeltas[a["Contraction"][x, y], a["Rules"]];

ContractionTable[state_GaussianState, sites_List, part_String] :=
  With[{modes = Tuples[{sites, SpinString}]},
   Switch[part,
    "Normal", Table[Contraction[state, cDag @@ x, c @@ y], {x, modes}, {y, modes}],
    "Anomalous", Table[Contraction[state, c @@ x, c @@ y], {x, modes}, {y, modes}]]];
ContractionTable[state_GaussianState, sites_List] :=
  With[{modes = Tuples[{sites, SpinString}]},
   Grid[{{Style["\[LeftAngleBracket]\!\(\*SuperscriptBox[\(c\), \(\[Dagger]\)]\)c\[RightAngleBracket]", Bold], Style["\[LeftAngleBracket]cc\[RightAngleBracket]", Bold]},
     {TableForm[ContractionTable[state, sites, "Normal"], TableHeadings -> {cDag @@@ modes, c @@@ modes}],
      TableForm[ContractionTable[state, sites, "Anomalous"], TableHeadings -> {c @@@ modes, c @@@ modes}]}},
    Spacings -> {4, 1}, Alignment -> Top]];

(*Amplitudes on screen: \[Chi]_ij, \[Chi]^a_ij, ...*)
Format[\[Chi][i_, j_], StandardForm] := Interpretation[Subscript["\[Chi]", Row[{i, j}]], \[Chi][i, j]];
Format[\[Chi][i_, j_, a_String], StandardForm] := Interpretation[Subsuperscript["\[Chi]", Row[{i, j}], a], \[Chi][i, j, a]];
Format[\[CapitalDelta][i_, j_], StandardForm] := Interpretation[Subscript["\[CapitalDelta]", Row[{i, j}]], \[CapitalDelta][i, j]];
Format[\[CapitalDelta][i_, j_, a_String], StandardForm] := Interpretation[Subsuperscript["\[CapitalDelta]", Row[{i, j}], a], \[CapitalDelta][i, j, a]];

(*Complex conjugates and squared moduli on screen, as ToTeX writes them: Conjugate[x] as x^*, and x^* x as |x|^2, for the
  amplitudes ($complexHeads) and for symbols. Only the display changes: the boxes interpret as the expression.*)
conjugateShownQ[x_] := MatchQ[Unevaluated[x], _Symbol | (h_[__] /; MemberQ[$complexHeads, h])];
MakeBoxes[Conjugate[x_?conjugateShownQ], StandardForm] :=
  With[{b = MakeBoxes[x, StandardForm]}, InterpretationBox[SuperscriptBox[b, "*"], Conjugate[x]]];
MakeBoxes[modulusSquared[x_, k_], StandardForm] :=
  With[{b = MakeBoxes[x, StandardForm], e = ToString[2 k]}, InterpretationBox[SuperscriptBox[RowBox[{"\[LeftBracketingBar]", b, "\[RightBracketingBar]"}], e], (Conjugate[x] x)^k]];
(*the factors of a product with their exponents; a pair x^m Conjugate[x]^n (m, n positive integers) gives |x|^(2 min(m, n))*)
factorExponents[t_Times] := Merge[Replace[List @@ t, {Power[b_, e_] :> (b -> e), b_ :> (b -> 1)}, {1}], Total];
pairedQ[t_Times] := With[{f = factorExponents[t]},
   AnyTrue[Keys[f], conjugateShownQ[#] && IntegerQ[f[#]] && f[#] > 0 && IntegerQ[Lookup[f, Key[Conjugate[#]], 0]] &&
      Lookup[f, Key[Conjugate[#]], 0] > 0 &]];
pairUp[t_Times] := Module[{f = factorExponents[t], mods = {}, k},
   Do[If[conjugateShownQ[x] && IntegerQ[f[x]] && f[x] > 0 && IntegerQ[Lookup[f, Key[Conjugate[x]], 0]] && Lookup[f, Key[Conjugate[x]], 0] > 0,
      k = Min[f[x], f[Conjugate[x]]]; f[x] -= k; f[Conjugate[x]] -= k; AppendTo[mods, modulusSquared[x, k]]], {x, Keys[f]}];
   Times @@ Join[KeyValueMap[#1^#2 &, f], mods]];
MakeBoxes[t_Times?pairedQ, StandardForm] := MakeBoxes[#, StandardForm] &[pairUp[t]];

(* ========================================================================================== *)
(*  5. Wick's theorem and mean field                                                          *)
(* ========================================================================================== *)

(*iWick[ops, budget]: sum over the pairings, memoized for one call (the contraction iContract is set by the caller).
  budget = number of contractions between different sites still allowed*)
iResetWick[f_] := (
   Clear[iWick];
   iContract = f;
   iWick[{}, _] = 1;
   iWick[opsList_List] := iWick[opsList, Infinity];
   iWick[opsList_List, budget_] := iWick[opsList, budget] =
     If[OddQ[Length[opsList]], 0,
      Module[{firstOp = First[opsList], sumRes = 0, val, subRes, idx, cross},
       Do[
        cross = opsList[[idx, 1]] =!= firstOp[[1]];
        If[! cross || budget > 0,
         val = iContract[firstOp, opsList[[idx]]];
         If[val =!= 0,
          subRes = iWick[Delete[opsList, {{1}, {idx}}], If[cross, budget - 1, budget]];
          If[subRes =!= 0, sumRes += (-1)^idx val subRes]]],
        {idx, 2, Length[opsList]}];
       sumRes]]);

iStateQ[s_] := MatchQ[s, GaussianState[_Association]] || (Message[WickVEV::nostate, s]; False);

(*"MaxChiOrder": old name of "MaxOrder"*)
Options[WickVEV] = {"NormalOrder" -> True, "MaxOrder" -> Infinity, "MaxChiOrder" -> Automatic};
WickVEV[expr_, opts : OptionsPattern[]] := WickVEV[expr, $State, opts];
WickVEV[expr_, state_GaussianState, OptionsPattern[]] /; iStateQ[state] :=
  Module[{terms, maxOrder = Replace[OptionValue["MaxChiOrder"], Automatic :> OptionValue["MaxOrder"]]},
   iResetWick[state["Contraction"]];
   terms = If[TrueQ[OptionValue["NormalOrder"]], iNormalOrderedTermList[expr],
     With[{lad = ToLadderBasis[expr]}, iSpinVariables[lad]; iTermList[lad]]];
   CleanDeltas[Expand[Total[Map[If[#[[2]] === {}, #[[1]], #[[1]] iWick[#[[2]], maxOrder]] &, terms]]], state["Rules"]]];

VerifyNormalOrder[expr_] := VerifyNormalOrder[expr, $State];
VerifyNormalOrder[expr_, state_GaussianState] :=
  Module[{diff, spins, cases},
   diff = WickVEV[expr, state] - WickVEV[expr, state, "NormalOrder" -> False];
   spins = iSpinVariables[ToLadderBasis[expr]];
   cases = If[spins === {}, {diff}, Map[diff /. Thread[spins -> #] &, Tuples[SpinString, Length[spins]]]];
   If[AllTrue[cases, Simplify[Expand[#]] === 0 &], 0, iSimplify[diff]]];

iSimplify[expr_] := FullSimplify[CleanDeltas[expr], TransformationFunctions -> {Automatic, CleanDeltas}];

(*Mean field: in every normal-ordered product two operators stay out, in every way (in their order, which stays normal),
  and the rest is contracted. Odd terms have zero expectation value and zero variation in Gaussian states: 0*)
iMeanFieldList[ops_List] :=
  Module[{m = Length[ops], p, q},
   Which[
    OddQ[m], 0,
    m <= 2, iNCM[ops],
    True,
    Sum[(-1)^(p + q + 1) iWick[Delete[ops, {{p}, {q}}]] ops[[p]] ** ops[[q]], {p, 1, m - 1}, {q, p + 1, m}] -
     (m/2 - 1) iWick[ops]]];

MeanField[expr_] := MeanField[expr, $State];
MeanField[expr_, state_GaussianState] /; iStateQ[state] :=
  (iResetWick[state["Contraction"]];
   CleanDeltas[Expand[Total[Map[If[#[[2]] === {}, #[[1]], #[[1]] iMeanFieldList[#[[2]]]] &, iNormalOrderedTermList[expr]]]],
    state["Rules"]]);

SetOptions[FullSimplify, ExcludedForms -> _NonCommutativeMultiply];

(* ========================================================================================== *)
(*  6. Local operators and the spin-charge basis                                              *)
(* ========================================================================================== *)

SpinOperator[site_, 1] := (SpOp[site] + SmOp[site])/2;
SpinOperator[site_, 2] := (SpOp[site] - SmOp[site])/(2 I);
SpinOperator[site_, 3] := SzOp[site];

iSpinWeight["\[UpArrow]"] = 1;
iSpinWeight["\[DownArrow]"] = 2;
iSpinWeight[_] = 3;

(*Local order: site, then cDag before c, then spin (cDag: \[UpArrow]\[DownArrow]; c: \[DownArrow]\[UpArrow]); a total order also with symbolic spins*)
iLocalOpLessEqual[op1_, op2_] :=
  Which[
   ! MatchQ[op1, $ladderPattern] || ! MatchQ[op2, $ladderPattern], OrderedQ[{op1, op2}],
   op1[[1]] =!= op2[[1]], OrderedQ[{op1[[1]], op2[[1]]}],
   Head[op1] =!= Head[op2], Head[op1] === cDag,
   Length[op1] =!= 2 || Length[op2] =!= 2, OrderedQ[{op1, op2}],
   True,
   With[{w1 = iSpinWeight[op1[[2]]], w2 = iSpinWeight[op2[[2]]], isDag = Head[op1] === cDag},
    Which[
     w1 < w2, isDag,
     w1 > w2, ! isDag,
     True, If[isDag, OrderedQ[{op1[[2]], op2[[2]]}], OrderedQ[{op2[[2]], op1[[2]]}]]]]];

(*Bubble sort with the fermion sign*)
iFermionSortSign[opsList_List] :=
  Module[{list = opsList, sign = 1, swapped = True, isZero = False, op1, op2, idx},
   While[swapped && ! isZero,
    swapped = False;
    Do[
     op1 = list[[idx]]; op2 = list[[idx + 1]];
     If[op1 === op2 && MatchQ[op1, $ladderPattern], isZero = True; Break[]];
     If[! iLocalOpLessEqual[op1, op2],
      list[[{idx, idx + 1}]] = {op2, op1};
      If[MatchQ[op1, $ladderPattern] && MatchQ[op2, $ladderPattern], sign = -sign];
      swapped = True],
     {idx, 1, Length[list] - 1}]];
   If[isZero, {0, {}}, {sign, list}]];

(*Ladder products -> composite operators*)
$localDensityRules = {
   NonCommutativeMultiply[a___, cDag[s_, "\[UpArrow]"], cDag[s_, "\[DownArrow]"], c[s_, "\[DownArrow]"], c[s_, "\[UpArrow]"], b___] :> NonCommutativeMultiply[a, DOp[s], b],
   NonCommutativeMultiply[a___, cDag[s_, "\[UpArrow]"], cDag[s_, "\[DownArrow]"], b___] :> NonCommutativeMultiply[a, Bdag[s], b],
   NonCommutativeMultiply[a___, c[s_, "\[DownArrow]"], c[s_, "\[UpArrow]"], b___] :> NonCommutativeMultiply[a, B[s], b],
   NonCommutativeMultiply[a___, cDag[s_, "\[UpArrow]"], c[s_, "\[DownArrow]"], b___] :> NonCommutativeMultiply[a, SpOp[s], b],
   NonCommutativeMultiply[a___, cDag[s_, "\[DownArrow]"], c[s_, "\[UpArrow]"], b___] :> NonCommutativeMultiply[a, SmOp[s], b],
   NonCommutativeMultiply[a___, cDag[s_, sp_], c[s_, sp_], b___] :> NonCommutativeMultiply[a, nOp[s, sp], b],
   NonCommutativeMultiply[a___, NonCommutativeMultiply[x__], b___] :> NonCommutativeMultiply[a, x, b],
   NonCommutativeMultiply[x_] :> x,
   NonCommutativeMultiply[] :> 1};

(*The site-by-site reordering is exact only on normal-ordered input: guaranteed here*)
GroupLocalOperators[expr_] :=
  Expand[iMapOpsTerms[Function[ops, With[{res = iFermionSortSign[ops]}, res[[1]] iNCM[res[[2]]]]], NormalOrder[expr]]];

(*From a normal-ordered polynomial in c, cDag to the spin-charge basis*)
iSCFromLadder[lad_] :=
  Module[{sub},
   sub = Expand[iMapOpsTerms[Function[ops, With[{r = iFermionSortSign[ops]}, r[[1]] iNCM[r[[2]]]]], lad]] //. $localDensityRules;
   sub = sub /. {
      nOp[site_, "\[UpArrow]"] :> (nTotOp[site]/2 + SzOp[site]),
      nOp[site_, "\[DownArrow]"] :> (nTotOp[site]/2 - SzOp[site])};
   iNCExpand[sub] //. {
     NonCommutativeMultiply[pre___, NonCommutativeMultiply[mid__], post___] :> NonCommutativeMultiply[pre, mid, post],
     NonCommutativeMultiply[x_] :> x}];

iSCTerm[{coef_, factors_List}] :=
  If[factors === {}, coef,
   With[{fz = iSiteFactorize[factors]},
    If[fz === $Failed,
     coef iSCFromLadder[NormalOrder[iNCM[factors]]],
     coef fz[[1]] iNCM[iLocalSCAt @@@ fz[[2]]]]]];

iToSC[expr_] :=
  Module[{e = iExpandMultiSite[expr], spinVars, sVar},
   spinVars = iSpinVariables[e];
   If[spinVars =!= {},
    sVar = First[spinVars];
    CleanDeltas[Expand[
      KroneckerDelta[sVar, "\[UpArrow]"] iToSC[e /. sVar -> "\[UpArrow]"] +
       KroneckerDelta[sVar, "\[DownArrow]"] iToSC[e /. sVar -> "\[DownArrow]"]]],
    iNCExpand[Total[iSCTerm /@ iTermList[e]]] //. {
      NonCommutativeMultiply[pre___, NonCommutativeMultiply[mid__], post___] :> NonCommutativeMultiply[pre, mid, post],
      NonCommutativeMultiply[x_] :> x}]];

Options[ToSpinChargeBasis] = {"Recognize" -> True};
ToSpinChargeBasis[expr_, OptionsPattern[]] :=
  With[{res = iToSC[expr]}, If[TrueQ[OptionValue["Recognize"]], iRecognizeNonlocal[res], res]];

$ladderRules = {
   DOp[s_] :> cDag[s, "\[UpArrow]"] ** cDag[s, "\[DownArrow]"] ** c[s, "\[DownArrow]"] ** c[s, "\[UpArrow]"],
   Bdag[s_] :> cDag[s, "\[UpArrow]"] ** cDag[s, "\[DownArrow]"],
   B[s_] :> c[s, "\[DownArrow]"] ** c[s, "\[UpArrow]"],
   SpOp[s_] :> cDag[s, "\[UpArrow]"] ** c[s, "\[DownArrow]"],
   SmOp[s_] :> cDag[s, "\[DownArrow]"] ** c[s, "\[UpArrow]"],
   SzOp[s_] :> 1/2 (cDag[s, "\[UpArrow]"] ** c[s, "\[UpArrow]"] - cDag[s, "\[DownArrow]"] ** c[s, "\[DownArrow]"]),
   nTotOp[s_] :> cDag[s, "\[UpArrow]"] ** c[s, "\[UpArrow]"] + cDag[s, "\[DownArrow]"] ** c[s, "\[DownArrow]"],
   nOp[s_, sp_] :> cDag[s, sp] ** c[s, sp],
   nOp[s_] :> cDag[s, "\[UpArrow]"] ** c[s, "\[UpArrow]"] + cDag[s, "\[DownArrow]"] ** c[s, "\[DownArrow]"]};

ToLadderBasis[expr_] :=
  Module[{e = expr},
   If[Length[$DefinedOperators] > 0 && ! FreeQ[e, iNonlocalPattern[]],
    e = e /. (op : iNonlocalPattern[]) :> iNonlocalLadder[Head[op], List @@ op]];
   iNCExpand[e /. $ladderRules] //. {
     NonCommutativeMultiply[pre___, NonCommutativeMultiply[mid__], post___] :> NonCommutativeMultiply[pre, mid, post],
     NonCommutativeMultiply[x_] :> x}];

(* ========================================================================================== *)
(*  7. Defined operators                                                                      *)
(* ========================================================================================== *)

(*A defined operator: head, number of sites, definition (function of the sites, written with known operators) and display.
  Recognition: the off-diagonal terms of the definition (with S+, S-, B, Bdag, c, cDag) must all be present with the same
  common coefficient, possibly times operators on other sites; the diagonal part (Sz, n, D) can leave a remainder:
  e.g. XXZ -> J S.S + (Jz - J) Sz Sz. With "Required" -> All the whole operator must be present.*)

iNonlocalPattern[] := Alternatives @@ (Blank /@ Keys[$DefinedOperators]);

iZeroQ[x_] := x === 0 || Simplify[x] === 0;

iSitesOf[f_] := If[KeyExistsQ[$DefinedOperators, Head[f]], List @@ f, {First[f]}];

iResetNonlocalCache[] := (
   Clear[iNonlocalLadder];
   iNonlocalLadder[h_, sites_List] := iNonlocalLadder[h, sites] = NormalOrder[$DefinedOperators[h]["Definition"] @@ sites];
   iResetLocalCaches[]);

Options[DefineOperator] = {"Recognize" -> Automatic, "Required" -> Automatic};
DefineOperator[head : (_Symbol | AdjointOperator[_Symbol]), arity_Integer?Positive, definition_, display_, OptionsPattern[]] := (
   If[! MemberQ[$compositePattern, Blank[head]],
    $compositePattern = Append[$compositePattern, Blank[head]];
    $opPattern = Append[$opPattern, Blank[head]]];
   $DefinedOperators[head] = <|"Arity" -> arity, "Definition" -> definition,
     "Sites" -> Array[Unique[opSite] &, arity],
     "Recognize" -> TrueQ[Replace[OptionValue["Recognize"], Automatic -> (arity >= 2)]],
     "RequiredOption" -> OptionValue["Required"], "Display" -> display|>;
   iResetNonlocalCache[]; Quiet[Unset[iDefinedHermitianQ[head]]];
   With[{h = head, dsp = display},
    Format[h[args__], StandardForm] := With[{d = dsp[args]}, Interpretation[d, h[args]]];
    If[arity == 1,
     h[sites_List] := h @@ sites;
     h[s1_, s2__] := NonCommutativeMultiply @@ (h /@ {s1, s2})]];
   With[{adj = iAdjointOf[head, arity, definition, display, OptionValue["Recognize"], OptionValue["Required"]]},
    $DefinedOperators[head] = Append[$DefinedOperators[head], "Adjoint" -> adj]];
   head);

(*The adjoint of a defined operator, found once on placeholder sites: "Hermitian", {"Permutation", p} when it is the same operator
  with the sites permuted, otherwise {"Head", AdjointOperator[head]}, a compact operator defined here with a dagger on the display;
  the adjoint of AdjointOperator[head] is head*)
iZeroOperatorQ[e_] := AllTrue[Values[OperatorCoefficients[e]], iZeroQ];
iAdjointOf[AdjointOperator[h_], ___] := {"Head", h};
iAdjointOf[head_, arity_, definition_, display_, recog_, req_] :=
  Module[{s = $DefinedOperators[head]["Sites"], d, dag, perm},
   d = NormalOrder[definition @@ s]; dag = NormalOrder[iDag[d, {}]];
   perm = If[arity > 4, None, SelectFirst[Permutations[Range[arity]], iZeroOperatorQ[dag - NormalOrder[definition @@ s[[#]]]] &, None]];
   Which[
    perm === Range[arity], "Hermitian",
    perm =!= None, {"Permutation", perm},
    True,
    With[{def = definition, dsp = display},
     DefineOperator[AdjointOperator[head], arity, Function[Null, Dagger[def[##]]], Function[Null, Superscript[dsp[##], "\[Dagger]"]],
      "Recognize" -> recog, "Required" -> req]];
    {"Head", AdjointOperator[head]}]];

(*Model of the operator in the spin-charge basis, computed at the first use (on placeholder sites)*)
iNonlocalData[h_] :=
  Module[{d = $DefinedOperators[h], ph, tmpl, req, key},
   If[KeyExistsQ[d, "Template"], d,
    ph = d["Sites"];
    tmpl = iTermList[ToSpinChargeBasis[d["Definition"] @@ ph, "Recognize" -> False]];
    req = If[d["RequiredOption"] === All, tmpl,
      With[{off = Select[tmpl, Function[m, ! AllTrue[m[[2]], MatchQ[#, $diagonalPattern] &]]]},
       If[off === {}, tmpl, off]]];
    key = SelectFirst[SortBy[req, Length[#[[2]]] &], Function[m, SubsetQ[First /@ m[[2]], ph]], None];
    If[key === None, Message[DefineOperator::nokey, h]];
    d = Join[d, <|"Template" -> tmpl, "Required" -> req, "Key" -> key|>];
    $DefinedOperators[h] = d;
    d]];

(*All the ways of finding the factors of the key term in factors: {placeholder rules -> sites, positions used}*)
iMatchKey[keyFactors_, factors_] :=
  Module[{cands, out = {}, rules, combo},
   cands = Map[
     Function[kf, Select[Range[Length[factors]],
       Function[q, Head[factors[[q]]] === Head[kf] && Length[factors[[q]]] === Length[kf] &&
         Rest[List @@ factors[[q]]] === Rest[List @@ kf]]]],
     keyFactors];
   If[FreeQ[cands, {}],
    Do[
     If[DuplicateFreeQ[combo],
      rules = DeleteDuplicates[Thread[(First /@ keyFactors) -> (First /@ factors[[combo]])]];
      If[DuplicateFreeQ[Keys[rules]] && DuplicateFreeQ[Values[rules]], AppendTo[out, {rules, combo}]]],
     {combo, Tuples[cands]}]];
   out];

(*Finds the term with the factors target (on the sites of the model) and rest (elsewhere); checks the coefficient if expected =!= None*)
(*The coefficients are compared in the same order of the factors: the sign of the permutation of the odd ones
  (fermion operators) between the product target ** rest and the stored term*)
iOddQ[f_] := MatchQ[f, $ladderPattern] || With[{p = iFactorSiteParity[f]}, ListQ[p] && p[[2]] == 1];
iOrderSign[from_List, to_List] :=
  With[{a = Select[from, iOddQ], b = Select[to, iOddQ]}, If[Length[a] < 2, 1, Signature[Flatten[Position[a, #, {1}, 1] & /@ b]]]];
iFindTerm[terms_, index_, target_, rest_, siteSet_, expected_] :=
  Module[{cands},
   cands = Select[Lookup[index, Key[Sort[Join[target, rest]]], {}],
     Function[q, Select[terms[[q, 2]], Function[f, ! IntersectingQ[iSitesOf[f], siteSet]]] === rest]];
   Which[
    cands === {}, Missing[],
    expected === None, First[cands],
    iZeroQ[terms[[First[cands], 1]] - expected iOrderSign[Join[target, rest], terms[[First[cands], 2]]]], First[cands],
    True, Missing[]]];

(*Subtracts alpha * model * rest and adds alpha * defined operator * rest*)
iApplyExtraction[terms_, index_, template_, sigma_, rest_, siteSet_, alpha_, newOp_] :=
  Module[{t = terms, touched = {}, target, q, m},
   Do[
    target = m[[2]] /. sigma;
    q = iFindTerm[terms, index, target, rest, siteSet, None];
    If[MissingQ[q],
     AppendTo[t, {-alpha m[[1]], Join[target, rest]}],
     t[[q, 1]] = t[[q, 1]] - alpha m[[1]] iOrderSign[Join[target, rest], terms[[q, 2]]]; AppendTo[touched, q]],
    {m, template}];
   AppendTo[t, {alpha, Prepend[rest, newOp]}];
   Delete[t, List /@ Select[touched, iZeroQ[t[[#, 1]]] &]]];

(*One extraction of the operator h (or $Failed). orderedOnly: only sites in order*)
iExtractOnce[terms_List, h_, orderedOnly_] :=
  Module[{data = iNonlocalData[h], ph, key, index, result = $Failed, gamma, factors, sigma, sites, siteSet, rest, alpha, pos, mt},
   If[data["Key"] === None, $Failed,
    ph = data["Sites"]; key = data["Key"];
    index = GroupBy[Range[Length[terms]], Sort[terms[[#, 2]]] &];
    Do[
     {gamma, factors} = terms[[pos]];
     Do[
      sigma = mt[[1]];
      sites = ph /. sigma;
      siteSet = Values[sigma];
      rest = Delete[factors, List /@ mt[[2]]];
      If[(! orderedOnly || OrderedQ[sites]) && ! AnyTrue[rest, IntersectingQ[iSitesOf[#], siteSet] &],
       alpha = gamma/(key[[1]] iOrderSign[Join[key[[2]] /. sigma, rest], factors]);
       If[AllTrue[data["Required"], ! MissingQ[iFindTerm[terms, index, #[[2]] /. sigma, rest, siteSet, alpha #[[1]]]] &],
        result = iApplyExtraction[terms, index, data["Template"], sigma, rest, siteSet, alpha, h @@ sites];
        Break[]]],
      {mt, iMatchKey[key[[2]], factors]}];
     If[result =!= $Failed, Break[]],
     {pos, Length[terms]}];
    result]];

(*Recognizes repeatedly: first the operators with more sites; first the choices with sites in order (S_i.S_j, not S_j.S_i)*)
iRecognizeNonlocal[expr_] :=
  With[{cand = Select[Keys[$DefinedOperators], $DefinedOperators[#]["Recognize"] &]},
   If[cand === {}, expr,
    Module[{terms = iTermList[expr], heads, changed = True, res, ord, h},
     heads = Flatten[SortBy[GatherBy[Reverse[cand], $DefinedOperators[#]["Arity"] &], -$DefinedOperators[First[#]]["Arity"] &]];
     While[changed,
      changed = False;
      Do[
       res = iExtractOnce[terms, h, ord];
       If[res =!= $Failed, terms = res; changed = True; Break[]],
       {ord, {True, False}}, {h, heads}]];
     Total[#[[1]] iNCM[#[[2]]] & /@ terms]]]];

ResetNormalOrder[];

(* ========================================================================================== *)
(*  8. Peierls phases and fluxes                                                              *)
(* ========================================================================================== *)

(*A constraint on the phases is a DEFINITION (one phase eliminated per loop), attached to the symbols (UpValues)*)
DefinePeierlsHopping[amp_Symbol, phase_Symbol] := (
   amp /: Subscript[amp, a_, b_] := amp Exp[I Subscript[phase, a, b]];
   phase /: Subscript[phase, a_, a_] := 0;
   phase /: Subscript[phase, a_, b_] /; ! OrderedQ[{a, b}] := -Subscript[phase, b, a];
   {amp, phase});

DefineFlux[phase_Symbol, loop_List, flux_] /; Length[loop] >= 3 :=
  Module[{links, sum, cand, coef, val},
   links = Partition[loop, 2, 1, 1];
   sum = Expand[Total[Apply[Subscript[phase, #1, #2] &, links, {1}]]];
   cand = SelectFirst[Reverse[Sort /@ links],
     Function[pr, With[{a = pr[[1]], b = pr[[2]]},
       Hold[Evaluate[Subscript[phase, a, b]]] === Hold[Subscript[phase, a, b]] && Coefficient[sum, Subscript[phase, a, b]] =!= 0]],
     None];
   If[cand === None,
    Message[DefineFlux::dep, loop, sum]; $Failed,
    coef = Coefficient[sum, Subscript[phase, cand[[1]], cand[[2]]]];
    val = Expand[(flux - (sum - coef Subscript[phase, cand[[1]], cand[[2]]]))/coef];
    With[{a = cand[[1]], b = cand[[2]], v = val},
     phase /: Subscript[phase, a, b] = v;
     HoldForm[Subscript[phase, a, b]] -> v]]];

(* ========================================================================================== *)
(*  9. LaTeX                                                                                  *)
(* ========================================================================================== *)

iTeXArgs[a_String] := With[{x = StringDelete[a, " "]}, If[StringCount[x, ","] >= 2, x, StringReplace[x, "," -> ""]]];
$texOps = {"cDag" -> "c^\\dagger", "nTotOp" -> "n", "nOp" -> "n", "SzOp" -> "S^z", "SpOp" -> "S^+", "SmOp" -> "S^-",
   "DOp" -> "D", "Bdag" -> "B^\\dagger", "EmptyProjector" -> "P^0", "SingleOccupationProjector" -> "P^1",
   "GutzwillerProjector" -> "P^G"};
(*Operators defined by the user (and adjoints) are written with their display; those of the package have rules below*)
iUserDisplayQ[AdjointOperator[_]] := True;
iUserDisplayQ[h_Symbol] := ! MemberQ[Join[Keys[$texOps], {"SdotS", "SdotScrossS"}], SymbolName[h]];
iUserDisplayQ[_] := False;
iTeXDisplay[expr_] := If[Length[$DefinedOperators] == 0, expr,
   expr /. op : iNonlocalPattern[] /; iUserDisplayQ[Head[op]] :> $DefinedOperators[Head[op]]["Display"] @@ op];
ToTeX[expr_] :=
  Fold[StringReplace[#1, #2] &, ToString[TeXForm[iTeXDisplay[expr]]], {
    Map[With[{name = #[[1]], tex = #[[2]]},
       RegularExpression["\\\\text\\{" <> name <> "\\}\\(([^()]*)\\)"] :> tex <> "_{" <> iTeXArgs["$1"] <> "}"] &, $texOps],
    RegularExpression["(?<![A-Za-z\\\\])c\\(([^()]*)\\)"] :> "c_{" <> iTeXArgs["$1"] <> "}",
    RegularExpression["(?<![A-Za-z\\\\])B\\(([^()]*)\\)"] :> "B_{" <> iTeXArgs["$1"] <> "}",
    RegularExpression["\\\\text\\{SdotScrossS\\}\\(([^(),]*),([^(),]*),([^(),]*)\\)"] :>
      "\\mathbf S_{$1}\\cdot(\\mathbf S_{$2}\\times\\mathbf S_{$3})",
    RegularExpression["\\\\text\\{SdotS\\}\\(([^(),]*),([^(),]*)\\)"] :> "\\mathbf S_{$1}\\cdot\\mathbf S_{$2}",
    RegularExpression["\\\\(chi|Delta) \\(([^(),]*),([^(),]*),\\\\text\\{([xyz])\\}\\)"] :>
      "{\\" <> "$1" <> "^{" <> "$4" <> "}_{" <> iTeXArgs["$2"] <> iTeXArgs["$3"] <> "}}",
    RegularExpression["\\\\(chi|Delta|rho|kappa) \\(([^()]*)\\)"] :> "\\" <> "$1" <> "_{" <> iTeXArgs["$2"] <> "}",
    {"\\uparrow " -> "\\uparrow", "\\downarrow " -> "\\downarrow"}}];

iTeXSpin["\[UpArrow]"] = "\\uparrow"; iTeXSpin["\[DownArrow]"] = "\\downarrow"; iTeXSpin[s_] := ToString[TeXForm[s]];
iTeXLabels[l_List] := StringRiffle[iTeXSpin /@ l, If[Length[l] > 2, ",", ""]];
iOpTeX[cDag[a___]] := "c^\\dagger_{" <> iTeXLabels[{a}] <> "}";
iOpTeX[c[a___]] := "c_{" <> iTeXLabels[{a}] <> "}";

(*All the pairings of a list, with the sign of the permutation (the convention of iWick)*)
iMatchings[{}] = {{1, {}}};
iMatchings[l_List] := If[OddQ[Length[l]], {},
   Flatten[Table[Map[{(-1)^k #[[1]], Prepend[#[[2]], {l[[1]], l[[k]]}]} &, iMatchings[Delete[l, {{1}, {k}}]]], {k, 2, Length[l]}], 1]];

iMonomial[expr_] :=
  With[{t = iTermList[ToLadderBasis[expr]]},
   If[Length[t] == 1 && t[[1, 2]] =!= {} && AllTrue[t[[1, 2]], MatchQ[#, $ladderPattern] &], t[[1]], $Failed]];

WickPairings[expr_] := WickPairings[expr, $State];
WickPairings[expr_, state_GaussianState] :=
  With[{m = iMonomial[expr], f = state["Contraction"]},
   If[m === $Failed, Message[WickTeX::mono, expr]; $Failed,
    Select[
     Map[<|"Sign" -> #[[1]], "Pairs" -> #[[2]],
        "Value" -> CleanDeltas[Times @@ (f[m[[2, #[[1]]]], m[[2, #[[2]]]]] & /@ #[[2]]), state["Rules"]]|> &,
      iMatchings[Range[Length[m[[2]]]]]],
     #["Value"] =!= 0 &]]];

(*Heights of the brackets: a contraction is drawn above every shorter one it overlaps*)
iWickLevels[pairs_] :=
  Module[{lev = <||>},
   Scan[Function[pq, lev[pq] = 1 + Max[0, Lookup[lev, Select[Keys[lev], ! (#[[2]] < pq[[1]] || #[[1]] > pq[[2]]) &]]]],
    SortBy[pairs, {#[[2]] - #[[1]], #[[1]]} &]];
   lev];
iWickTeX[ops_, pairs_] :=
  Module[{lev = iWickLevels[pairs], mark = <||>},
   If[Max[Values[lev]] > 9, Message[WickTeX::levels]];
   Scan[(mark[#[[1]]] = lev[#]; mark[#[[2]]] = lev[#]) &, pairs];
   "\\wick{" <> StringRiffle[MapIndexed[If[KeyExistsQ[mark, First[#2]], "\\c" <> ToString[mark[First[#2]]] <> " ", ""] <> iOpTeX[#1] &, ops], " "] <> "}"];

iJoinTeX[l_List] := StringReplace[StringRiffle[l, " + "], "+ -" -> "- "];

Options[WickTeX] = {"Values" -> True, "Result" -> True};
WickTeX[expr_, opts : OptionsPattern[]] := WickTeX[expr, $State, opts];
WickTeX[expr_, state_GaussianState, OptionsPattern[]] :=
  Module[{m = iMonomial[expr], prs, coef, ops, wrap, parts},
   If[m === $Failed, Message[WickTeX::mono, expr]; Return[$Failed]];
   {coef, ops} = m;
   prs = WickPairings[expr, state];
   wrap[s_] := If[coef === 1, s, ToTeX[coef] <> "\\left(" <> s <> "\\right)"];
   parts = {If[coef === 1, "", ToTeX[coef]] <> "\\langle " <> StringRiffle[iOpTeX /@ ops, " "] <> "\\rangle",
     If[prs === {}, "0", wrap[iJoinTeX[(If[#["Sign"] == 1, "", "-"] <> iWickTeX[ops, #["Pairs"]]) & /@ prs]]]};
   If[TrueQ[OptionValue["Values"]] && prs =!= {}, AppendTo[parts, wrap[iJoinTeX[ToTeX[#["Sign"] #["Value"]] & /@ prs]]]];
   If[TrueQ[OptionValue["Result"]] && prs =!= {}, AppendTo[parts, ToTeX[Expand[coef Total[#["Sign"] #["Value"] & /@ prs]]]]];
   StringRiffle[DeleteDuplicatesBy[Reverse[parts], StringDelete[#, " "] &] // Reverse, " = "]];

(* ========================================================================================== *)
(*  10. Library of operators                                                                  *)
(* ========================================================================================== *)

DOp[sites_List] := DOp @@ sites;
DOp[s1_, s2__] := NonCommutativeMultiply @@ (DOp /@ {s1, s2});

DefineOperator[EmptyProjector, 1, (1 - nTotOp[#] + DOp[#]) &, Subsuperscript["P", #, "0"] &];
DefineOperator[SingleOccupationProjector, 1, (nTotOp[#] - 2 DOp[#]) &, Subsuperscript["P", #, "1"] &];

(*The order of the factors is site1, site2, site3: valid also for coinciding sites*)
DefineOperator[SdotS, 2,
  (SzOp[#1] ** SzOp[#2] + (SpOp[#1] ** SmOp[#2] + SmOp[#1] ** SpOp[#2])/2) &,
  Row[{Subscript[Style["S", Bold], #1], "\[CenterDot]", Subscript[Style["S", Bold], #2]}] &];

DefineOperator[SdotScrossS, 3,
  (I/2 (SzOp[#1] ** SpOp[#2] ** SmOp[#3] - SzOp[#1] ** SmOp[#2] ** SpOp[#3] +
        SpOp[#1] ** SmOp[#2] ** SzOp[#3] - SpOp[#1] ** SzOp[#2] ** SmOp[#3] -
        SmOp[#1] ** SpOp[#2] ** SzOp[#3] + SmOp[#1] ** SzOp[#2] ** SpOp[#3])) &,
  Row[{Subscript[Style["S", Bold], #1], "\[CenterDot](", Subscript[Style["S", Bold], #2], "\[Times]", Subscript[Style["S", Bold], #3], ")"}] &];

(* ========================================================================================== *)
(*  11. The banner and the update check                                                       *)
(* ========================================================================================== *)

(* the banner shown when the package is loaded in a notebook (as FunctionalDiagrams does): the logo (the small badge of
   Logo/FermiWick_logo_small.svg, the outer contraction on the gradient of the package, y upwards), the attribution
   "FermiWick, by Nazzareno Africani" word for word (term (b) of NOTICE), the version, license, how to cite, with links to
   the repository and the website; printed by the end of this file only when $Notebooks *)
$fwStart = RGBColor[100/255, 42/255, 140/255]; $fwEnd = RGBColor[176/255, 37/255, 115/255];
$fwAccent = RGBColor[138/255, 40/255, 128/255]; $fwTint = RGBColor[247/255, 181/255, 218/255];
(* a colour for light notebooks and one for dark notebooks (LightDarkSwitched, Mathematica 14.3 or later) *)
lightDarkColor[light_, dark_] := If[Names["System`LightDarkSwitched"] =!= {}, Symbol["LightDarkSwitched"][light, dark], light];
badgePoints[] := With[{r = 28.}, Join @@ MapThread[Function[{c, a0}, Table[c + r {Cos[a], Sin[a]}, {a, a0, a0 + Pi/2, Pi/16}]],
     {{{128 - r, r}, {128 - r, 128 - r}, {r, 128 - r}, {r, r}}, {-Pi/2, 0, Pi/2, Pi}}]];
quarterArc[c_, r_, a0_, a1_] := Table[c + r {Cos[a], Sin[a]}, {a, a0, a1, (a1 - a0)/12}];
packageLogo[size_: 30] := Module[{pts = badgePoints[], x0 = 26, x1 = 102, y = 36., top = 100., rr = 18.},
   Graphics[{
     Polygon[pts, VertexColors -> (Blend[{$fwStart, $fwEnd}, (#[[1]] + 128 - #[[2]])/256] & /@ pts)],
     {EdgeForm[{White, Opacity[0.16], Thickness[2/128]}], FaceForm[None], Polygon[badgePoints[]]},
     White, CapForm["Butt"], JoinForm["Round"], Thickness[11/128],
     Line[Join[{{x0, y + 9.5}}, quarterArc[{x0 + rr, top - rr}, rr, Pi, Pi/2], quarterArc[{x1 - rr, top - rr}, rr, Pi/2, 0],
       {{x1, y + 10.5}}]],
     EdgeForm[{White, Thickness[2/128], JoinForm["Round"]}], Polygon[{{54, top}, {70, top + 10}, {70, top - 10}}],
     EdgeForm[], Disk[{x0, y}, 11], Thickness[6/128], Circle[{x1, y}, 8.5]},
    PlotRange -> {{0, 128}, {0, 128}}, PlotRangePadding -> None, ImagePadding -> None, ImageSize -> size,
    BaselinePosition -> Scaled[0.2]]];
$repositoryURL = "https://github.com/n4sd4q/FermiWick";
$websiteURL = "https://n4sd4q.github.io/FermiWick/";
(* the name in the accent of the logo (its tint in a dark notebook), the rest in grey; U+2010 (a hyphen, not a minus
   sign), U+00B7 and U+00A0 by character code, so that the file stays plain ASCII. Each group is unbreakable and ends with
   its separator, so that in a narrow window the lines wrap after a separator, never before one *)
bannerDot[] := With[{nbsp = FromCharacterCode[160]}, nbsp <> nbsp <> FromCharacterCode[183] <> nbsp <> " "];
bannerGroup[items__] := Style[Row[{items}], LineBreakWithin -> False];
bannerAccent[] := lightDarkColor[$fwAccent, $fwTint];
bannerGrey[] := lightDarkColor[GrayLevel[0.45], GrayLevel[0.6]];
loadBanner[release_ : None] :=
  Module[{dot = bannerDot[], group = bannerGroup, accent = bannerAccent[], grey = bannerGrey[], file = $repositoryURL <> "/blob/main/"},
   Style[Grid[{{packageLogo[30],
       Column[{
         Row[{group[Style["FermiWick", Bold, 14, accent], Style[", by Nazzareno Africani", 14, accent], Style[dot, 12, grey]],
           Style["fermion operators, as symbols and as matrices", 12, grey]}],
         Style[Row[{group["version " <> $FermiWickVersion, dot],
             group[StringReplace["GPL-3.0-or-later", "-" -> FromCharacterCode[8208]], dot],
             group["please cite it: ", Hyperlink["CITATION.cff", file <> "CITATION.cff"], dot],
             group[Hyperlink["GitHub", $repositoryURL], dot],
             Hyperlink["Website", $websiteURL]}], 11, grey],
         If[AssociationQ[release], updateNotice[release], Nothing]},
        Spacings -> 0.1]}},
      Alignment -> {Left, Center}, Spacings -> {0.7, 0}], FontFamily -> "Source Sans Pro", LineIndent -> 0]];

(* the update check. In a notebook (never in a script), at most once a day, the package asks GitHub for its latest release
   in the background (URLSubmit: loading never waits for the network) and, if it is newer than the version loaded, the
   banner gets a third line with the PacletInstall line of its paclet. Any failure (offline, a proxy, a time limit, the rate
   limit of GitHub, an answer without the fields) is silent. The handler of the answer runs after the evaluation that
   loaded the package: it writes the banner again in its own cell, found by a tag given to it at load, without moving the
   selection. FermiWickUpdateCheck[] asks by hand, and waits for the answer *)
$latestReleaseURL = "https://api.github.com/repos/n4sd4q/FermiWick/releases/latest";
latestReleaseRequest[] := HTTPRequest[$latestReleaseURL, <|"Headers" -> {"Accept" -> "application/vnd.github+json"}|>];
$updateCheckDirectory := FileNameJoin[{$UserBaseDirectory, "ApplicationData", "FermiWick"}];
updateCheckDate[] := DateString[{"Year", "-", "Month", "-", "Day"}];
networkOptions[f_, opts_List] := FilterRules[opts, Options[f]];

(* the numbers of a version "0.10.0" or of a tag "v0.10.0"; $Failed for anything else (a pre-release "v0.3.0-rc1") *)
versionNumbers[s_String] := If[StringMatchQ[s, RegularExpression["v?[0-9]+\\.[0-9]+\\.[0-9]+"]],
   FromDigits /@ StringSplit[StringReplace[s, StartOfString ~~ "v" -> ""], "."], $Failed];
versionNumbers[_] := $Failed;
(* whether version a is newer than version b, number by number (0.10.0 is newer than 0.9.0) *)
newerVersionQ[a_, b_] := With[{x = versionNumbers[a], y = versionNumbers[b]}, ListQ[x] && ListQ[y] && Order[y, x] === 1];

(* the latest release in the answer of GitHub (JSON): <|"Version" -> "0.2.0", "Paclet" -> url of its .paclet, "Page" -> url
   of its page|>, or $Failed for a draft or a pre-release, a malformed tag or no paclet of this repository. ImportString
   reads bytes: the text is encoded in UTF-8 first. The paclet's address must be one of this repository's downloads, made
   of plain characters, since it is shown inside a line of code to copy *)
latestReleaseOf[body_String] := Module[{json, tag, version, paclets},
   json = Quiet[ImportString[FromCharacterCode[ToCharacterCode[body, "UTF8"]], "RawJSON"]];
   If[! AssociationQ[json], Return[$Failed]];
   tag = Lookup[json, "tag_name"];
   If[TrueQ[Lookup[json, "draft"]] || TrueQ[Lookup[json, "prerelease"]] || ! ListQ[versionNumbers[tag]], Return[$Failed]];
   version = StringReplace[tag, StartOfString ~~ "v" -> ""];
   paclets = Cases[Lookup[json, "assets", {}], KeyValuePattern["browser_download_url" -> url_String] /;
       StringMatchQ[url, RegularExpression[StringReplace[$repositoryURL, "." -> "\\."] <> "/releases/download/" <>
          StringReplace[tag, "." -> "\\."] <> "/[A-Za-z0-9._-]+\\.paclet"]] :> url];
   If[paclets === {}, $Failed,
    <|"Version" -> version, "Paclet" -> First[paclets], "Page" -> $repositoryURL <> "/releases/tag/" <> tag|>]];
latestReleaseOf[_] := $Failed;
installLine[release_Association] := "PacletInstall[\"" <> release["Paclet"] <> "\"]";

(* the third line of the banner: the newer version, a link to its page and the line that installs it, as code to copy *)
updateNotice[release_Association] := With[{dot = bannerDot[], group = bannerGroup},
   Style[Row[{group[Style["version " <> release["Version"] <> " is available", Bold, bannerAccent[]], dot],
      group[Hyperlink["what is new", release["Page"]], dot],
      With[{url = release["Paclet"]}, Style[Defer[PacletInstall[url]], ShowStringCharacters -> True,
        PrivateFontOptions -> {"OperatorSubstitution" -> False}]]}], 11, bannerGrey[]]];
updateNoticeText[release_Association] :=
  "FermiWick " <> release["Version"] <> " is available (" <> $FermiWickVersion <> " is loaded): " <> installLine[release];

(* whether the check at load is due today: the file LastUpdateCheck.txt of dir holds the date of the last check; the date of
   today is written into it before asking. Any error reading or writing it: the check is due *)
updateCheckDueQ[dir_String, today_String] := Module[{file = FileNameJoin[{dir, "LastUpdateCheck.txt"}], last, stream},
   last = Quiet[ReadString[file]];
   If[StringQ[last] && StringTrim[last] === today, Return[False]];
   Quiet[If[! DirectoryQ[dir], CreateDirectory[dir, CreateIntermediateDirectories -> True]];
    stream = OpenWrite[file];
    If[Head[stream] === OutputStream, WriteString[stream, today <> "\n"]; Close[stream]]];
   True];

$updateCheckTask = None;
submitUpdateCheck[nb_, tag_String] := Quiet[
   $updateCheckTask = URLSubmit[latestReleaseRequest[],
     HandlerFunctions -> <|"TaskFinished" -> (showUpdateNotice[nb, tag, #] &)|>,
     HandlerFunctionsKeys -> {"StatusCode", "Body"}, Sequence @@ networkOptions[URLSubmit, {TimeConstraint -> 30}]]];
showUpdateNotice[nb_, tag_String, answer_] := Quiet[Module[{release, cells},
    release = If[Lookup[answer, "StatusCode"] === 200, latestReleaseOf[Lookup[answer, "Body"]], $Failed];
    If[AssociationQ[release] && newerVersionQ[release["Version"], $FermiWickVersion],
     cells = Cells[nb, CellTags -> tag];
     If[ListQ[cells] && Length[cells] == 1,
      NotebookWrite[First[cells], ReplacePart[NotebookRead[First[cells]], 1 -> BoxData[ToBoxes[loadBanner[release]]]]]]]]];

FermiWickUpdateCheck::latest = "The latest release could not be read from GitHub (`1`).";
FermiWickUpdateCheck[] := Module[{response, status, release},
   response = Quiet[URLRead[latestReleaseRequest[], Sequence @@ networkOptions[URLRead, {TimeConstraint -> 10}]]];
   status = If[Head[response] === HTTPResponse, response["StatusCode"], None];
   release = If[status === 200, latestReleaseOf[response["Body"]], $Failed];
   If[! AssociationQ[release],
    Message[FermiWickUpdateCheck::latest, Which[status === None, "no answer", status =!= 200, "status " <> ToString[status],
       True, "no release with a paclet in the answer"]];
    $Failed,
    If[newerVersionQ[release["Version"], $FermiWickVersion],
     Print[If[TrueQ[$Notebooks], Style[updateNotice[release], FontFamily -> "Source Sans Pro"], updateNoticeText[release]]]];
    <|"InstalledVersion" -> $FermiWickVersion, "LatestVersion" -> release["Version"], "InstallLine" -> installLine[release]|>]];

(* the check at load is on unless $FermiWickUpdateCheck was set before loading (kept by the start of the package) *)
If[! ValueQ[$FermiWickUpdateCheck], $FermiWickUpdateCheck = True];

(* in a notebook (not in scripts, which stay quiet), the banner, and the update check when it is due *)
If[TrueQ[$Notebooks],
  If[$FermiWickUpdateCheck =!= False && $AllowInternet =!= False && updateCheckDueQ[$updateCheckDirectory, updateCheckDate[]],
   With[{tag = "FermiWickBanner" <> ToString[Round[10^6 AbsoluteTime[]]]},
    CellPrint[ExpressionCell[loadBanner[], "Print", CellTags -> tag]];
    submitUpdateCheck[EvaluationNotebook[], tag]],
   Print[loadBanner[]]]];

End[];
EndPackage[];
