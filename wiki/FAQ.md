Short answers to the usual questions. Every block on this page runs after

```mathematica
Needs["FermiWick`"]
```

## Which function do I need?

| To... | use |
|---|---|
| bring a product to normal order | `NormalOrder` |
| read a result as densities and spins | `ToSpinChargeBasis` (and back: `ToLadderBasis`) |
| take an expectation value in a Gaussian state | `WickVEV` |
| decouple an interaction | `MeanField`, then `OperatorCoefficients` for the fields |
| solve a cluster exactly | `FockMatrix`, `LowestStates`, `FockExpectation` (``FermiWick`Fock` ``) |
| write LaTeX | `ToTeX`, `WickTeX` |

## Why does `c[i, s] ** cDag[j, s]` give no delta?

Different symbols are different sites: `i` and `j` never coincide. Use one symbol for one site, and write the case
`i = j` apart when you need it:

```mathematica
{NormalOrder[c[i, SpinUp] ** cDag[j, SpinUp]], NormalOrder[c[i, SpinUp] ** cDag[i, SpinUp]]}
```

$$\\{-c^\dagger_{j\uparrow}c_{i\uparrow},1-c^\dagger_{i\uparrow}c_{i\uparrow}\\}$$

## How do I put numbers in a result?

Replace the amplitudes; their conjugates follow:

```mathematica
WickVEV[SdotS[1, 2], SpinSingletState[]] /. χ[1, 2] -> 1/3
```

$$-\frac{1}{6}$$

For a complex value give the conjugate too, `{χ[1, 2] -> z, Conjugate[χ[1, 2]] -> Conjugate[z]}`, or write the
amplitude as `r Exp[I φ]` with real `r` and `φ`.

## Why is my expectation value zero?

The state decides what survives: a contraction that its symmetries forbid is zero. The chirality, for example, vanishes
in every state with time reversal; `ContractionTable` shows what a state keeps
([Declaring a state](Declaring-a-State)).

## Can a site have its own filling or magnetization?

Yes, with rules on the ladder operators, which can depend on the site in any way
([Declaring a state](Declaring-a-State#by-rules)).

## Can the spin be a symbol?

Yes: `c[i, s]` with a symbol `s` is split into the two spins when needed, and a sum over the spin is `Sum[..., {s, $Spins}]`:

```mathematica
NormalOrder[Sum[cDag[1, s] ** c[1, s] ** cDag[1, r] ** c[1, r], {s, $Spins}, {r, $Spins}]] - NormalOrder[nTotOp[1] ** nTotOp[1]]
```

$$0$$

## How large a cluster can I diagonalise?

Sectors of up to about a million states: a Hubbard ring of twelve sites at half filling (853 776 states) is written as
a matrix in about four seconds, and its ground state found by Arnoldi. Lattice symmetries (momentum, point group) are
not used.

## How do I cite FermiWick?

As described in [`CITATION.cff`](https://github.com/n4sd4q/FermiWick/blob/main/CITATION.cff); GitHub's "Cite this
repository" button gives the BibTeX.