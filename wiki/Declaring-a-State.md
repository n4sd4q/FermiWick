Wick's theorem needs a Gaussian state: a state fixed by its two-point contractions, $\langle c^\dagger c\rangle$ and
$\langle c\thinspace c\rangle$. You declare it by the symmetries you assume, or by rules; `WickVEV`, `MeanField` and `WickTeX`
take it as their last argument. Every block on this page runs after

```mathematica
Needs["FermiWick`"]
```

## By its symmetries

`GaussianState[symmetries]` is the most general Gaussian state with those symmetries. Its contractions are written in
spin channels: hopping χ and pairing Δ, and their spin parts χ<sup>a</sup> (magnetization) and Δ<sup>a</sup> (triplet
pairing), $a = x, y, z$ ([Conventions](Conventions)). Each symmetry sets channels to zero or makes them real:

| Symmetry | What it does |
|---|---|
| `"SpinRotation"` | only χ and Δ: no magnetization, no triplet pairing |
| `"SpinRotationZ"` | rotations about z: only the z components of χ<sup>a</sup> and Δ<sup>a</sup> |
| `"Charge"` | particle number conserved: no pairing |
| `"TimeReversal"` | χ and Δ real, χ<sup>a</sup> and Δ<sup>a</sup> imaginary |

`ContractionTable` shows what survives, here $\langle c^\dagger_a c_b\rangle$ for the modes of two sites, in the order
1↑, 1↓, 2↑, 2↓:

```mathematica
ContractionTable[GaussianState[{"SpinRotation", "Charge", "TimeReversal"}], {1, 2}, "Normal"]
```

```math
\left(
\begin{array}{cccc}
 \chi_{11} & 0 & \chi_{12} & 0 \\
 0 & \chi_{11} & 0 & \chi_{12} \\
 \chi_{12} & 0 & \chi_{22} & 0 \\
 0 & \chi_{12} & 0 & \chi_{22} \\
\end{array}
\right)
```

The on-site contractions have their own options: `"Filling" -> ν` ($\langle n_{i\sigma}\rangle$), `"Magnetization" -> m`
($\langle S^z_i\rangle$) and `"OnSitePairing" -> η` ($\langle c_{i\uparrow}c_{i\downarrow}\rangle$).

## States invariant under spin rotations

`SpinSingletState[]` is the state invariant under spin rotations that conserves the charge, at half filling;
`"Pairing" -> Δ` adds singlet pairing. In it, hopping and pairing enter the exchange with the same weight:

```mathematica
WickVEV[SdotS[1, 2], SpinSingletState["Pairing" -> Δ]]
```

$$-\frac{3}{2} \Delta_{12} \Delta_{12}^{\ast}-\frac{3}{2} \chi_{12} \chi_{12}^{\ast}$$

## Magnetic states

With rotations about z only, the state can carry a magnetization along z:

```mathematica
magnetic = GaussianState[{"SpinRotationZ", "Charge"}, "Filling" -> 1/2, "Magnetization" -> m];
WickVEV[SzOp[1] ** SzOp[2], magnetic]
```

$$-\frac{1}{2} \chi_{12} \chi_{12}^{\ast}-\frac{1}{2} {\chi^{z}_{12}} {\chi^{z}_{12}}^{\ast}+m^2$$

## By rules

`GaussianState[{rule, ...}]` gives the contractions by rules on the ladder operators, with any spin structure and any
dependence on the sites. Here a hopping g that keeps the spin and a hopping h that flips it:

```mathematica
flip = GaussianState[{{cDag[i_, s_], c[j_, r_]} :> If[i === j, KroneckerDelta[s, r]/2, If[s === r, g, h]]}];
WickVEV[SzOp[1] ** SzOp[2], flip]
```

$$\frac{h^2}{2}-\frac{g^2}{2}$$

and a magnetization that alternates from site to site, on numbered sites:

```mathematica
neel = GaussianState[{{cDag[i_, s_], c[j_, r_]} /; i === j :> KroneckerDelta[s, r] (1/2 + (-1)^i If[s === SpinUp, m, -m])}];
WickVEV[SzOp[1] ** SzOp[2], neel]
```

$$-m^2$$

A contraction that no rule gives is zero. `GeneralState[]` is the state with no assumption at all, written with spin
components instead of channels (`"Anomalous" -> κ` adds pairing).

## Choosing the contractions

`StateSelect[state, crit]` keeps only the contractions for which `crit[op1, op2]` is true: `LocalPairQ` keeps those on one
site, `! AnomalousPairQ[##] &` drops the pairing:

```mathematica
WickVEV[SdotS[1, 2], StateSelect[SpinSingletState["Pairing" -> Δ], ! AnomalousPairQ[##] &]]
```

$$-\frac{3}{2} \chi_{12} \chi_{12}^{\ast}$$

It is a filter, not a symmetry: the result need not be a physical state. `WickVEV[expr, state, "MaxOrder" -> m]` keeps
the terms with at most `m` contractions between different sites; the chirality needs three:

```mathematica
WickVEV[SdotScrossS[1, 2, 3], GaussianState[{"SpinRotation", "Charge"}], "MaxOrder" -> 2]
```

$$0$$

## One contraction, and the default state

`Contraction[state, op1, op2]` gives a single contraction:

```mathematica
Contraction[GaussianState["SpinRotation"], c[1, SpinUp], c[2, SpinDown]]
```

$$\Delta_{12}$$

`$State` is the state used when none is given (`SpinSingletState[]` unless you set it): `WickVEV[expr]` is
`WickVEV[expr, $State]`.