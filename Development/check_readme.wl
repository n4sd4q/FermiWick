(* Checks the first example of README.md: its code blocks are run in order, in one fresh kernel with the paclet of this
   repository, and the LaTeX of what each block gives must be the formula printed after it (as GitHub receives it:
   ^{\ast} for ^*, \\{ for \{).     math -script Development/check_readme.wl *)

root = ParentDirectory[DirectoryName[$InputFileName]];
PacletDirectoryLoad[FileNameJoin[{root, "FermiWick"}]];
Needs["FermiWick`"]; Needs["FermiWick`Fock`"];   (* the contexts before the code is parsed *)
readme = Import[FileNameJoin[{root, "README.md"}], "Text", CharacterEncoding -> "UTF8"];
section = First[StringCases[readme, "## A first example" ~~ Shortest[x__] ~~ "\n## " :> x]];
(*the code blocks, each with the formula that follows it before the next block (or None)*)
pieces = StringCases[section, RegularExpression["```mathematica\\n([\\s\\S]*?)```"] :> "$1"];
formulas = Map[First[StringCases[#, RegularExpression["^```[\\s\\S]*?```([\\s\\S]*?)(?:```|$)"] :> "$1"], ""] &,
   StringCases[section, RegularExpression["```mathematica[\\s\\S]*?(?=```mathematica|$)"]]];
formulaOf[s_String] := First[StringCases[s, RegularExpression["\\$\\$([^$]+)\\$\\$"] :> "$1"], None];
fromGitHub[f_String] := StringReplace[f, {"^{\\ast}" -> "^*", "\\\\{" -> "\\{", "\\\\}" -> "\\}", "\\thinspace " -> "\\,"}];
docPost[t_String] := StringReplace[t, {"\\text{Jc}" -> "J_\\chi", "\\left\\{" -> "\\left\\{", ",\\frac" -> ",\\frac"}];
fail = 0;
Do[
  With[{code = pieces[[n]], expected = formulaOf[formulas[[n]]]},
   Module[{res, got},
    res = Check[ToExpression[code], $Failed];
    If[res === $Failed, fail++; Print["FAIL block ", n, ": messages"],
     If[expected === None, Print["ok   block ", n, " (no formula)"],
      got = docPost[ToTeX[res]];
      If[StringDelete[got, " "] === StringDelete[fromGitHub[expected], " "], Print["ok   block ", n],
       fail++; Print["FAIL block ", n, "\n  README: ", fromGitHub[expected], "\n  gives:  ", got]]]]]],
  {n, Length[pieces]}];
Print[If[fail == 0, "README: all blocks give what is shown", "README: " <> ToString[fail] <> " block(s) differ"]];
Exit[If[fail == 0, 0, 1]];
