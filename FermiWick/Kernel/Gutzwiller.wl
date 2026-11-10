(* ::Package:: *)

(* FermiWick`Gutzwiller`: Gutzwiller projectors, projected operators and their expectation values to first order in the
   contractions between different sites. *)

BeginPackage["FermiWick`Gutzwiller`", {"FermiWick`"}];

GutzwillerProjector::usage = "GutzwillerProjector[i] = Sqrt[2] (Cos[\[Theta]] P1_i + Sin[\[Theta]] (P0_i + D_i)), the Gutzwiller operator at half filling \
with the variational angle \[Theta] (\[Theta] = 0: projection on single occupation); with several sites it is the product.";
\[Theta]::usage = "\[Theta] is the variational angle of GutzwillerProjector.";
LocalSandwich::usage = "LocalSandwich[P, {i, j, ...}, X] = P_i\[CenterDot]P_j\[CenterDot]...\[CenterDot]X\[CenterDot]P_i\[CenterDot]P_j\[CenterDot]...; LocalSandwich[PL, sites, X, PR] uses PL on the left \
and PR on the right.";
GutzwillerVEV::usage = "GutzwillerVEV[expr, state] is WickVEV[expr, state] to first order in the contractions between different sites \
(the hopping amplitudes of the state), simplified. GutzwillerVEV[expr] uses $State.";
GutzwillerSimplify::usage = "GutzwillerSimplify[expr] simplifies expr with the rules of the deltas.";

Begin["`Private`"];

DefineOperator[GutzwillerProjector, 1,
  (Sqrt[2] (Cos[\[Theta]] SingleOccupationProjector[#] + Sin[\[Theta]] (EmptyProjector[#] + DOp[#]))) &,
  Subsuperscript["P", #, "G"] &];

LocalSandwich[proj_, sites_List, x_] := LocalSandwich[proj, sites, x, proj];
LocalSandwich[projL_, sites_List, x_, projR_] := CenterDot @@ Join[projL /@ sites, {x}, projR /@ sites];

GutzwillerSimplify[expr_] := FullSimplify[CleanDeltas[expr], TransformationFunctions -> {Automatic, CleanDeltas}];

GutzwillerVEV[expr_, opts : OptionsPattern[WickVEV]] := GutzwillerVEV[expr, $State, opts];
GutzwillerVEV[expr_, state_GaussianState, opts : OptionsPattern[WickVEV]] :=
  Module[{vev, scaled, s, amps = Alternatives @@ Replace[state["Amplitudes"], _Missing -> {}]},
   vev = WickVEV[expr, state, opts, "MaxOrder" -> 1];
   (*first order in the amplitudes between different sites, also when they appear conjugated*)
   scaled = Normal[Series[vev /. {x : Conjugate[amps[a_, b_, ___]] /; a =!= b :> s x, x : amps[a_, b_, ___] /; a =!= b :> s x},
       {s, 0, 1}]] /. s -> 1;
   GutzwillerSimplify[scaled]];

End[];
EndPackage[];
