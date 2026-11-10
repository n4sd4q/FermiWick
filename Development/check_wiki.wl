(* Checks a page of the wiki (wiki/NAME.md): its ```mathematica blocks are run in order, in one fresh kernel with the
   paclet of this repository, and each must give what the page shows after it, before the next block: a display formula
   $$...$$ is compared with the LaTeX of the result (ToTeX, written as GitHub needs it: ^{\ast} for ^*, \\{ for \{), a
   ```math block (for formulas on several lines, a matrix: GitHub leaves its backslashes alone) with the plain LaTeX, a
   ```text block with the result as a terminal prints it (ToString); a block with neither shows no output. A block after
   the line <!-- not run --> is skipped (installation lines, paths). A block that prints a message fails.
     math -script Development/check_wiki.wl wiki/NAME.md          checks the page
     math -script Development/check_wiki.wl wiki/NAME.md --fill   writes the outputs where the page has "?" in place
                                                                  of the formula or of the text, then checks *)

root = ParentDirectory[DirectoryName[$InputFileName]];
PacletDirectoryLoad[FileNameJoin[{root, "FermiWick"}]];
Needs["FermiWick`"]; Needs["FermiWick`Fock`"]; Needs["FermiWick`Gutzwiller`"];
args = If[Length[$ScriptCommandLine] > 1, Rest[$ScriptCommandLine], $CommandLine];   (* wolfram -script fills only $CommandLine *)
fill = MemberQ[args, "--fill"];
file = First[Select[args, StringEndsQ[#, ".md"] &]];
page = Import[file, "Text", CharacterEncoding -> "UTF8"];

toGitHub[t_String] := StringReplace[t, {"^*" -> "^{\\ast}", "\\{" -> "\\\\{", "\\}" -> "\\\\}", "\\," -> "\\thinspace "}];
fromGitHub[t_String] := StringReplace[t, {"^{\\ast}" -> "^*", "\\\\{" -> "\\{", "\\\\}" -> "\\}", "\\thinspace " -> "\\,"}];
asText[r_] := StringTrim[ToString[r]];

(* the pieces of the page: the text before the first block, then each block with what follows it *)
pieces = StringSplit[page, "```mathematica\n" -> "\[Bullet]BLOCK\[Bullet]"];
blocks = {};
Do[If[pieces[[k]] === "\[Bullet]BLOCK\[Bullet]",
   With[{body = pieces[[k + 1]], before = If[k > 1, pieces[[k - 1]], ""]},
    With[{code = First[StringSplit[body, "```", 2]], after = Last[StringSplit[body, "```", 2]]},
     AppendTo[blocks, <|"Code" -> code, "Skip" -> StringEndsQ[StringTrim[before], "<!-- not run -->"],
       "Formula" -> First[StringCases[after, StartOfString ~~ Whitespace ~~ "$$" ~~ f : Shortest[__] ~~ "$$" :> f], None],
       "Text" -> First[StringCases[after, StartOfString ~~ Whitespace ~~ "```text\n" ~~ t : Shortest[___] ~~ "```" :> t], None],
       "Math" -> First[StringCases[after, StartOfString ~~ Whitespace ~~ "```math\n" ~~ t : Shortest[___] ~~ "```" :> t], None]|>]]]],
  {k, Length[pieces]}];

fail = 0; filled = page;
Do[With[{b = blocks[[n]]},
   If[b["Skip"], Continue[]];
   Module[{res = Check[ToExpression[b["Code"]], $Failed], got},
    Which[
     res === $Failed, fail++; Print["FAIL ", FileNameTake[file], " block ", n, ": messages"],
     StringQ[b["Formula"]],
      got = ToTeX[res];
      Which[
       fill && StringTrim[b["Formula"]] === "?",
        filled = StringReplace[filled, b["Code"] <> "```\n\n$$?$$" -> b["Code"] <> "```\n\n$$" <> toGitHub[got] <> "$$", 1],
       StringDelete[got, " "] =!= StringDelete[fromGitHub[b["Formula"]], " "],
        fail++; Print["FAIL ", FileNameTake[file], " block ", n, "\n  page:  ", fromGitHub[b["Formula"]], "\n  gives: ", got]],
     StringQ[b["Text"]],
      got = asText[res];
      Which[
       fill && StringTrim[b["Text"]] === "?",
        filled = StringReplace[filled, b["Code"] <> "```\n\n```text\n?\n```" -> b["Code"] <> "```\n\n```text\n" <> got <> "\n```", 1],
       got =!= StringTrim[b["Text"]],
        fail++; Print["FAIL ", FileNameTake[file], " block ", n, "\n  page:  ", StringTrim[b["Text"]], "\n  gives: ", got]],
     StringQ[b["Math"]],
      got = ToTeX[res];
      Which[
       fill && StringTrim[b["Math"]] === "?",
        filled = StringReplace[filled, b["Code"] <> "```\n\n```math\n?\n```" -> b["Code"] <> "```\n\n```math\n" <> got <> "\n```", 1],
       StringDelete[got, Whitespace] =!= StringDelete[b["Math"], Whitespace],
        fail++; Print["FAIL ", FileNameTake[file], " block ", n, "\n  page:  ", b["Math"], "\n  gives: ", got]],
     True, Null]]],
  {n, Length[blocks]}];

If[fill && filled =!= page, Export[file, filled, "Text", CharacterEncoding -> "UTF8"]; Print["filled ", FileNameTake[file]]];
Print[If[fail == 0, "ok   ", "FAIL "], FileNameTake[file], ": ", Count[blocks, b_ /; ! b["Skip"]], " blocks",
  If[fail == 0, "", ", " <> ToString[fail] <> " differ"]];
Exit[If[fail == 0, 0, 1]];
