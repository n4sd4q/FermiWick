(* Runs the examples of a document and writes, for each, its code and its output.
     SRC=<examples.wl> OUT=<outdir> math -script Development/run_examples.wl
   The examples file is a sequence of blocks, each starting with a line  (* ::Example:: name *), or  (* ::Example:: name text *)
   for an output shown as text, or  (* ::Example:: name tex *)  for an output that is LaTeX (shown as text and typeset);
   the blocks are evaluated in order in one session, in the context Global`. For each block: outdir/name.wl (the code,
   verbatim) and outdir/name.tex (the output in LaTeX) or outdir/name.txt (the output as text). The names of this script
   live in the context DocsRunner`, so that the examples can use any symbol. *)

Begin["DocsRunner`"];
{src, out} = {Environment["SRC"], Environment["OUT"]};
root = ParentDirectory[DirectoryName[$InputFileName]];
PacletDirectoryLoad[FileNameJoin[{root, "FermiWick"}]];
Needs["FermiWick`"]; Needs["FermiWick`Gutzwiller`"]; Needs["FermiWick`Fock`"];   (* the contexts before the code is parsed *)
SetOptions[$Output, PageWidth -> 110];

(* LaTeX: ToTeX of the package, plus the names used in the documents *)
docPost[t_String] := StringReplace[t, {
    RegularExpression["\\\\text\\{Hop\\}\\(([^()]*)\\)"] :> "T_{" <> StringReplace[StringDelete["$1", " "], "," -> ""] <> "}",
    RegularExpression["(?<![A-Za-z\\\\])m\\(([^()]*)\\)"] :> "m_{$1}",
    RegularExpression["\\\\text\\{nk\\}\\(([^()]*)\\)"] :> "n_{$1}",
    "\\text{eps0}" -> "\\varepsilon_0", "\\text{Jc}" -> "J_\\chi"}];
(* a ket of sites with one orbital each (StateToKet) written |0,↑↓⟩, without the braces of the sites, its arrows as
   symbols (no space of a relation after the bar) *)
ketPost[t_String] := StringReplace[t, RegularExpression["\\|([^|]*?)\\\\rangle"] :>
    "|" <> StringReplace["$1", {"\\uparrow" -> "{\\uparrow}", "\\downarrow" -> "{\\downarrow}"}] <> "\\rangle"];
toTeX[e_] := docPost[ketPost[ToTeX[e /. Ket[k : {{_} ..}] :> Ket[Flatten[k]]]]];
toText[s_String] := s;
toText[e_] := StringReplace[ToString[e /. FermiWick`GaussianState[a_Association] :> stateShown[Lookup[a, "Description", "..."]], InputForm,
     PageWidth -> 110], "DocsRunner`stateShown" | "stateShown" -> "GaussianState"];
shown[code_] := StringReplace[code, Thread[("\\[" <> # <> "]" & /@ {"UpArrow", "DownArrow", "CapitalDelta", "Theta", "Nu", "Eta", "Chi",
       "CenterDot", "Rho", "Kappa", "CurlyPhi", "CapitalPhi", "Sigma", "Mu", "CapitalOmega"}) ->
    {"\[UpArrow]", "\[DownArrow]", "\[CapitalDelta]", "\[Theta]", "\[Nu]", "\[Eta]", "\[Chi]", "\[CenterDot]", "\[Rho]", "\[Kappa]",
     "\[CurlyPhi]", "\[CapitalPhi]", "\[Sigma]", "\[Mu]", "\[CapitalOmega]"}]];

text = Import[src, "Text"];
blocks = StringCases[text, RegularExpression["\\(\\* ::Example:: (\\w+)( text| tex)? \\*\\)\\n([\\s\\S]*?)(?=\\(\\* ::Example::|\\z)"] :>
    {"$1", StringTrim["$2"], StringTrim["$3"]}];
If[! DuplicateFreeQ[blocks[[All, 1]]],
  Print["duplicate example names: ", Select[Tally[blocks[[All, 1]]], #[[2]] > 1 &][[All, 1]]]; Exit[1]];
If[! DirectoryQ[out], CreateDirectory[out]];
End[];

Do[
  With[{name = DocsRunner`blk[[1]], asText = DocsRunner`blk[[2]] === "text", asTeX = DocsRunner`blk[[2]] === "tex",
    code = DocsRunner`blk[[3]], out = DocsRunner`out},
   Export[FileNameJoin[{out, name <> ".wl"}], DocsRunner`shown[code], "Text", CharacterEncoding -> "UTF8"];
   DocsRunner`res = ToExpression[code];
   Which[
    DocsRunner`res === Null, Null,
    asTeX, Export[FileNameJoin[{out, name <> ".txt"}], DocsRunner`res, "Text", CharacterEncoding -> "UTF8"];
     Export[FileNameJoin[{out, name <> "_typeset.tex"}], DocsRunner`docPost[DocsRunner`res], "Text"],
    asText || AssociationQ[DocsRunner`res], Export[FileNameJoin[{out, name <> ".txt"}], DocsRunner`toText[DocsRunner`res], "Text", CharacterEncoding -> "UTF8"],
    True, Export[FileNameJoin[{out, name <> ".tex"}], DocsRunner`toTeX[DocsRunner`res], "Text"]];
   Print[name, ": ", StringTake[DocsRunner`toText[DocsRunner`res], UpTo[90]]]],
  {DocsRunner`blk, DocsRunner`blocks}];
