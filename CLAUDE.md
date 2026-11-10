# FermiWick — notes for Claude Code

A Wolfram Language package for fermion ladder operators (see `README.md`). Author: Nazzareno Africani (SISSA). Code,
comments and documents in English. Personal notes of a developer (resume point, local workflow) go in `CLAUDE.local.md`,
which is not versioned.

## Layout
- `FermiWick/` — the paclet: `PacletInfo.wl`, `Kernel/FermiWick.wl` (context ``FermiWick` ``), `Kernel/Gutzwiller.wl`
  (``FermiWick`Gutzwiller` ``), `Kernel/Fock.wl` (``FermiWick`Fock` ``: the operators as sparse matrices on Fock spaces and
  symmetry sectors, exact diagonalisation). `FockMatrix` reads an expression through
  `ToLadderBasis` and the private `iTermList` (no normal ordering): keep those two stable.
- The banner and the update check (as FunctionalDiagrams): the end of `Kernel/FermiWick.wl`; `SECURITY.md` describes the
  only access to the network; the tests read a saved answer of GitHub, `Tests/LatestRelease.json` (the release 0.2.0, trimmed),
  and load the package in fresh kernels.
- `Tests/FermiWick.wlt`, `Tests/Fock.wlt`, `Tests/RunTests.wls` — `math -script Tests/RunTests.wls`; every test must pass.
  The exact checks of the symbolic engine use `FockMatrix`; the Jordan–Wigner construction is written out once, in the
  first test of `Fock.wlt`.
- `Showcase/` and `Documentation/` — LaTeX sources; the examples (`examples.wl`, blocks `(* ::Example:: name *)`) are run
  by `Development/run_examples.wl`, which writes `examples/NAME.wl` (code) and `NAME.tex`/`NAME.txt` (output).
  `Showcase/build_pdf.sh`, `Documentation/build_manual.sh` compute the examples and typeset. Common style:
  `Development/style.tex` (on `n4sd4q.sty`, the Orchid colours of the palette). Logo: `Logo/` (SVG; the PDFs are made with rsvg-convert).
- `Examples.wl` — the guided tour (cells marked as in a package); `FermiWick/Examples.nb` is generated from it by
  `math -script Development/make_notebook.wl` (saved with Visible -> True). Edit the `.wl`, never the `.nb`.
- `Development/check_docs.sh` — runs every documented command again in fresh kernels (showcase, manual, tour, README) and
  compares with what the documents show; it must pass before a commit that touches the kernel or the documents.
- `Media/` — the demo video `demo.mp4` (1920×1080, 60 fps, H.264, about 50 s), its looping preview `demo_preview.webp`
  (960×540, 60 fps, at the top of the README, linking the video) and its poster `demo_poster.png` (the title card, which
  is also the end card, so that the video loops without a seam), made by `python3 Development/make_demo_video.py` (about
  eight minutes, most of it the WebP; `--frames` writes frames to `build/video/frames`, `--poster` only the poster). Its
  cards are one notebook session, listed in `CARDS`: cells of the README and computed examples of the showcase and the
  manual, code and outputs (never typed by hand); `--check` runs the session in a fresh kernel and compares every output.
  Run it again when those examples change, and look at the frames. The website shows the video under its hero.
- `Development/operators_table.tex` — the table of the composite operators, included by the showcase and the manual.
- `wiki/` — the GitHub wiki (practical guides and a FAQ, as FunctionalDiagrams): `Home.md`, `_Sidebar.md` (every page),
  `_Footer.md`, `Page-Name.md` (no `#` title but on Home), links `[text](Page-Name)`. `.github/workflows/wiki.yml` mirrors
  the folder into the wiki on every push to main that touches it: edit the pages here, never on GitHub.
  `Development/check_wiki.wl wiki/Page.md` runs the blocks of a page in a fresh kernel; each must give the `$$formula$$`
  (`ToTeX`, GitHub's escapes), the ```` ```math ```` block (formulas on several lines) or the ```` ```text ```` block
  (`ToString`) after it; `<!-- not run -->` before a block skips it; `--fill` writes the outputs where the page has `?`.
  Wiki math: no `\,` (`\thinspace`), no `<` or `>` (`\lt`, `\gt`: GitHub reads `<q` as a tag), no `)_`, `]_`, `}_`
  or `+_` (GitHub reads them as emphasis: subscript first).
- `index.html` — the website (GitHub Pages serves main: https://n4sd4q.github.io/FermiWick/), **generated** by
  `python3 Development/make_site.py` from the README, `PacletInfo.wl`, `CITATION.cff` and the PDFs (formulas typeset by
  pdflatex, inlined as SVG); never edit it, run the script after changing its sources (`check_docs.sh` checks it). `.nojekyll`
  (empty) lets Pages serve the files as they are.
- The common look of the packages: `Development/n4sd4q.sty` with `n4sd4q-palette.tex` (the style of the PDFs, read through
  `TEXINPUTS` by the build scripts), `Development/family.css` and `Development/palette.json` (the website) and the logos in
  `Logo/` are copies of the files of the hub (github.com/n4sd4q/n4sd4q.github.io, `design/`): change them there, then copy.
- License: GPL-3.0-or-later with the additional terms of `NOTICE` (as FunctionalDiagrams).

## Rules
- Commit after each step, with a clear message.
- Speed matters: profile before and after a change of the kernel; results must not change.
- Every command shown in the documents gives exactly the output shown after it (computed by the scripts).
- After a change: tests, then rebuild the showcase and the manual, look at the pages, update README and CHANGELOG.
- README math: GitHub's Markdown runs before MathJax, so no literal `*` (write `^{\ast}`), no `\,` (write
  `\thinspace`), braces as `\\{ \\}`. Check with `gh api markdown -f mode=gfm -f text="$(cat README.md)"`: every
  `<math-renderer>` must contain the formula as written. Pictures need a background that reads on light and dark pages.
- The README links the PDFs on GitHub Pages (`https://n4sd4q.github.io/FermiWick/Showcase/...pdf`), so that they open in
  the browser, not on GitHub's page of the file (as FunctionalDiagrams); `make_site.py` makes those links relative.

## Releases
- Semantic versioning (as FunctionalDiagrams). The version is the same in `FermiWick/PacletInfo.wl`, `$FermiWickVersion`,
  `CITATION.cff` (with `date-released`), the install line of the README and of the manual (the URL of the paclet), the
  title lines of the showcase and the manual (`\familytitle`), the wiki (`Installation.md`, `_Footer.md`) and `CHANGELOG.md`.
- The paclet: `CreatePacletArchive["FermiWick", "build"]` from a clean checkout (build/ is not versioned); test it unzipped,
  with `PacletDirectoryLoad`, never with `PacletInstall` on the author's machine. Then the tag vX.Y.Z and the GitHub release
  with the paclet and the two PDFs; `python3 Development/make_site.py` after every change of the README or the version.

