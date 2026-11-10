# FermiWick — notes for Claude Code

A Wolfram Language package for fermion ladder operators (see `README.md`). Author: Nazzareno Africani (SISSA). Code,
comments and documents in English. Personal notes of a developer (resume point, local workflow) go in `CLAUDE.local.md`,
which is not versioned.

## Layout
- `FermiWick/` — the paclet: `PacletInfo.wl`, `Kernel/FermiWick.wl` (context ``FermiWick` ``), `Kernel/Gutzwiller.wl`
  (``FermiWick`Gutzwiller` ``), `Kernel/Fock.wl` (``FermiWick`Fock` ``: the operators as sparse matrices on Fock spaces and
  symmetry sectors, exact diagonalisation). `FockMatrix` reads an expression through
  `ToLadderBasis` and the private `iTermList` (no normal ordering): keep those two stable.
- `Tests/FermiWick.wlt`, `Tests/Fock.wlt`, `Tests/RunTests.wls` — `math -script Tests/RunTests.wls`; every test must pass.
  The exact checks of the symbolic engine use `FockMatrix`; the Jordan–Wigner construction is written out once, in the
  first test of `Fock.wlt`.
- `Showcase/` and `Documentation/` — LaTeX sources; the examples (`examples.wl`, blocks `(* ::Example:: name *)`) are run
  by `Development/run_examples.wl`, which writes `examples/NAME.wl` (code) and `NAME.tex`/`NAME.txt` (output).
  `Showcase/build_pdf.sh`, `Documentation/build_manual.sh` compute the examples and typeset. Common style:
  `Development/style.tex` (accent #8E2A7B). Logo: `Logo/` (SVG; the PDF is made with rsvg-convert).
- `Examples.wl` — the guided tour (cells marked as in a package); `FermiWick/Examples.nb` is generated from it by
  `math -script Development/make_notebook.wl` (saved with Visible -> True). Edit the `.wl`, never the `.nb`.
- `Development/check_docs.sh` — runs every documented command again in fresh kernels (showcase, manual, tour, README) and
  compares with what the documents show; it must pass before a commit that touches the kernel or the documents.
- `Development/operators_table.tex` — the table of the composite operators, included by the showcase and the manual.
- License: GPL-3.0-or-later with the additional terms of `NOTICE` (as FunctionalDiagrams).

## Rules
- Commit after each step, with a clear message.
- Speed matters: profile before and after a change of the kernel; results must not change.
- Every command shown in the documents gives exactly the output shown after it (computed by the scripts).
- After a change: tests, then rebuild the showcase and the manual, look at the pages, update README and CHANGELOG.
- README math: GitHub's Markdown runs before MathJax, so no literal `*` (write `^{\ast}`), no `\,` (write
  `\thinspace`), braces as `\\{ \\}`. Check with `gh api markdown -f mode=gfm -f text="$(cat README.md)"`: every
  `<math-renderer>` must contain the formula as written. Pictures need a background that reads on light and dark pages.
