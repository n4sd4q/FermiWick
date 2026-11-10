Every result of FermiWick can be written in LaTeX, for a paper or for notes; Wick's theorem can be written out, with the
contractions drawn. Every block on this page runs after

```mathematica
Needs["FermiWick`"]
```

## Formulas

`ToTeX[expr]` writes ladder and local operators, spins and amplitudes as in a paper:

```mathematica
ToTeX[ToSpinChargeBasis[OperatorCommutator[SdotS[1, 2], SdotS[2, 3]]]]
```

```text
-i \mathbf S_{1}\cdot(\mathbf S_{2}\times\mathbf S_{3})
```

```mathematica
ToTeX[WickVEV[SdotS[1, 2], SpinSingletState["Pairing" -> Δ]]]
```

```text
-\frac{3}{2} \Delta_{12} \Delta_{12}^*-\frac{3}{2} \chi_{12} \chi_{12}^*
```

## Wick's theorem written out

`WickTeX[product, state]` writes $\langle O_1\cdots O_{2n}\rangle$ as the sum of its pairings, drawn as contractions
(the LaTeX package [simpler-wick](https://ctan.org/pkg/simpler-wick)), then their values and the result; the options
`"Values" -> False` and `"Result" -> False` stop earlier. Only the pairings that the state allows are drawn:

```mathematica
WickTeX[cDag[1, SpinUp] ** cDag[2, SpinDown] ** c[3, SpinDown] ** c[4, SpinUp], SpinSingletState["Pairing" -> Δ], "Values" -> False, "Result" -> False]
```

```text
\langle c^\dagger_{1\uparrow} c^\dagger_{2\downarrow} c_{3\downarrow} c_{4\uparrow}\rangle = \wick{\c1 c^\dagger_{1\uparrow} \c1 c^\dagger_{2\downarrow} \c1 c_{3\downarrow} \c1 c_{4\uparrow}} + \wick{\c2 c^\dagger_{1\uparrow} \c1 c^\dagger_{2\downarrow} \c1 c_{3\downarrow} \c2 c_{4\uparrow}}
```

`WickPairings[product, state]` gives the same pairings as data: their signs, pairs and values.

## In a document

The formulas need `amsmath` and, for the contractions, `simpler-wick`; compile twice, so that the contractions are
placed:

<!-- not run -->
```mathematica
Export["wick.tex", WickTeX[cDag[1, SpinUp] ** cDag[2, SpinDown] ** c[3, SpinDown] ** c[4, SpinUp], SpinSingletState["Pairing" -> Δ]], "Text"]
```

```latex
\documentclass{article}
\usepackage{amsmath}
\usepackage{simpler-wick}
\begin{document}
\[ \input{wick.tex} \]
\end{document}
```

Long results are best broken by hand, or with `breqn`'s `dmath` environment.