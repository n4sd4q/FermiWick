Wick's theorem needs a Gaussian state: a state fixed by its two-point contractions, $\langle c^\dagger c\rangle$ and
$\langle c\thinspace c\rangle$. The most general one has them all; the one you assume is declared by its symmetries, or by
rules. `WickVEV`, `MeanField` and `WickTeX` take it as their last argument. Every block on this page runs after

```mathematica
Needs["FermiWick`"]
```

## The general state

`GaussianState[]` has every contraction, one amplitude each: the hopping
$\chi^{\sigma\sigma'}_{ij}=\langle c^\dagger_{i\sigma}c_{j\sigma'}\rangle$ and the pairing
$\Delta^{\sigma\sigma'}_{ij}=\langle c_{i\sigma}c_{j\sigma'}\rangle$. `ContractionTable` shows them between the modes of some
sites, here 1↑, 1↓, 2↑, 2↓; first $\langle c^\dagger_a c_b\rangle$:

```mathematica
general = GaussianState[];
ContractionTable[general, {1, 2}, "Normal"]
```

```math
\left(
\begin{array}{cccc}
 {\chi^{\uparrow\uparrow}_{11}} & {\chi^{\uparrow\downarrow}_{11}} & {\chi^{\uparrow\uparrow}_{12}} & {\chi^{\uparrow\downarrow}_{12}} \\
 {\chi^{\uparrow\downarrow}_{11}}^* & {\chi^{\downarrow\downarrow}_{11}} & {\chi^{\downarrow\uparrow}_{12}} & {\chi^{\downarrow\downarrow}_{12}} \\
 {\chi^{\uparrow\uparrow}_{12}}^* & {\chi^{\downarrow\uparrow}_{12}}^* & {\chi^{\uparrow\uparrow}_{22}} & {\chi^{\uparrow\downarrow}_{22}} \\
 {\chi^{\uparrow\downarrow}_{12}}^* & {\chi^{\downarrow\downarrow}_{12}}^* & {\chi^{\uparrow\downarrow}_{22}}^* & {\chi^{\downarrow\downarrow}_{22}} \\
\end{array}
\right)
```

Nothing is fixed: on one site the diagonal entries are the densities, $\langle n_{1\uparrow}\rangle=\chi^{\uparrow\uparrow}_{11}$ and
$\langle n_{1\downarrow}\rangle=\chi^{\downarrow\downarrow}_{11}$, and the other one, $\langle c^\dagger_{1\uparrow}c_{1\downarrow}\rangle=\langle S^+_1\rangle$, a
transverse magnetization. Hermiticity, $\langle c^\dagger_bc_a\rangle=\langle c^\dagger_ac_b\rangle^{\ast}$, always holds: an amplitude is
written with its modes in canonical order, the other order as its conjugate. Then $\langle c_ac_b\rangle$, antisymmetric, so
that on one site only $\langle c_{1\uparrow}c_{1\downarrow}\rangle$ is left:

```mathematica
ContractionTable[general, {1, 2}, "Anomalous"]
```

```math
\left(
\begin{array}{cccc}
 0 & {\Delta^{\uparrow\downarrow}_{11}} & {\Delta^{\uparrow\uparrow}_{12}} & {\Delta^{\uparrow\downarrow}_{12}} \\
 -{\Delta^{\uparrow\downarrow}_{11}} & 0 & {\Delta^{\downarrow\uparrow}_{12}} & {\Delta^{\downarrow\downarrow}_{12}} \\
 -{\Delta^{\uparrow\uparrow}_{12}} & -{\Delta^{\downarrow\uparrow}_{12}} & 0 & {\Delta^{\uparrow\downarrow}_{22}} \\
 -{\Delta^{\uparrow\downarrow}_{12}} & -{\Delta^{\downarrow\downarrow}_{12}} & -{\Delta^{\uparrow\downarrow}_{22}} & 0 \\
\end{array}
\right)
```

## The symmetries, one at a time

A symmetry you assume removes contractions, or the indices that told them apart.

**Charge conserved** (`"Charge"`): no pairing, the anomalous table vanishes:

```mathematica
ContractionTable[GaussianState["Charge"], {1, 2}, "Anomalous"]
```

```math
\left(
\begin{array}{cccc}
 0 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 \\
\end{array}
\right)
```

**$S^z$ conserved** (`"SpinRotationZ"`, rotations about $z$): the hopping keeps the spin, one amplitude for each,
$\chi^\uparrow_{ij}$ and $\chi^\downarrow_{ij}$; the pairing joins opposite spins, $\Delta^{\uparrow\downarrow}_{ij}$ and
$\Delta^{\downarrow\uparrow}_{ij}$:

```mathematica
sz = GaussianState["SpinRotationZ"];
ContractionTable[sz, {1, 2}, "Normal"]
```

```math
\left(
\begin{array}{cccc}
 {\chi^{\uparrow}_{11}} & 0 & {\chi^{\uparrow}_{12}} & 0 \\
 0 & {\chi^{\downarrow}_{11}} & 0 & {\chi^{\downarrow}_{12}} \\
 {\chi^{\uparrow}_{12}}^* & 0 & {\chi^{\uparrow}_{22}} & 0 \\
 0 & {\chi^{\downarrow}_{12}}^* & 0 & {\chi^{\downarrow}_{22}} \\
\end{array}
\right)
```

```mathematica
ContractionTable[sz, {1, 2}, "Anomalous"]
```

```math
\left(
\begin{array}{cccc}
 0 & {\Delta^{\uparrow\downarrow}_{11}} & 0 & {\Delta^{\uparrow\downarrow}_{12}} \\
 -{\Delta^{\uparrow\downarrow}_{11}} & 0 & {\Delta^{\downarrow\uparrow}_{12}} & 0 \\
 0 & -{\Delta^{\downarrow\uparrow}_{12}} & 0 & {\Delta^{\uparrow\downarrow}_{22}} \\
 -{\Delta^{\uparrow\downarrow}_{12}} & 0 & -{\Delta^{\uparrow\downarrow}_{22}} & 0 \\
\end{array}
\right)
```

On one site $\chi^\uparrow_{ii}$ and $\chi^\downarrow_{ii}$ are the densities of each spin, and half their difference the
magnetization. In the exchange $\mathbf S_1\cdot\mathbf S_2$ of such a state the hopping of each spin enters, the spin
flip $\chi^\uparrow\chi^{\downarrow\ast}$, and the product of the magnetizations:

```mathematica
WickVEV[SdotS[1, 2], GaussianState[{"SpinRotationZ", "Charge"}]]
```

$$-\frac{1}{4} {\chi^{\downarrow}_{12}} {\chi^{\downarrow}_{12}}^{\ast}-\frac{1}{2} {\chi^{\downarrow}_{12}} {\chi^{\uparrow}_{12}}^{\ast}-\frac{1}{2} {\chi^{\uparrow}_{12}} {\chi^{\downarrow}_{12}}^{\ast}-\frac{1}{4} {\chi^{\uparrow}_{12}} {\chi^{\uparrow}_{12}}^{\ast}+\frac{1}{4} {\chi^{\downarrow}_{11}} {\chi^{\downarrow}_{22}}-\frac{1}{4} {\chi^{\downarrow}_{22}} {\chi^{\uparrow}_{11}}-\frac{1}{4} {\chi^{\downarrow}_{11}} {\chi^{\uparrow}_{22}}+\frac{1}{4} {\chi^{\uparrow}_{11}} {\chi^{\uparrow}_{22}}$$

**Spin rotations** (`"SpinRotation"`, SU(2)): the hopping is the same for both spins, $\chi_{ij}$; the pairing is a
singlet, $\Delta_{ij}=\langle c_{i\uparrow}c_{j\downarrow}\rangle=-\langle c_{i\downarrow}c_{j\uparrow}\rangle$, symmetric in the sites; no
magnetization:

```mathematica
su2 = GaussianState["SpinRotation"];
ContractionTable[su2, {1, 2}, "Normal"]
```

```math
\left(
\begin{array}{cccc}
 \chi_{11} & 0 & \chi_{12} & 0 \\
 0 & \chi_{11} & 0 & \chi_{12} \\
 \chi_{12}^* & 0 & \chi_{22} & 0 \\
 0 & \chi_{12}^* & 0 & \chi_{22} \\
\end{array}
\right)
```

```mathematica
ContractionTable[su2, {1, 2}, "Anomalous"]
```

```math
\left(
\begin{array}{cccc}
 0 & \Delta_{11} & 0 & \Delta_{12} \\
 -\Delta_{11} & 0 & -\Delta_{12} & 0 \\
 0 & \Delta_{12} & 0 & \Delta_{22} \\
 -\Delta_{12} & 0 & -\Delta_{22} & 0 \\
\end{array}
\right)
```

This is the state of most calculations, so it has a short name: `SpinSingletState[]` is
`GaussianState[{"SpinRotation", "Charge"}, "Filling" -> 1/2]`, and `SpinSingletState["Pairing" -> Δ]` is
`GaussianState["SpinRotation", "Filling" -> 1/2, "OnSitePairing" -> 0]` (the on-site values fixed, as below). In it,
hopping and pairing enter the exchange with the same weight:

```mathematica
WickVEV[SdotS[1, 2], SpinSingletState["Pairing" -> Δ]]
```

$$-\frac{3}{2} \Delta_{12} \Delta_{12}^{\ast}-\frac{3}{2} \chi_{12} \chi_{12}^{\ast}$$

**Time reversal** (`"TimeReversal"`): the components of opposite spins become conjugates of each other, and with SU(2)
the amplitudes are real. The scalar chirality, odd under time reversal, then vanishes:

```mathematica
ContractionTable[GaussianState[{"SpinRotation", "TimeReversal"}], {1, 2}, "Normal"]
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

```mathematica
WickVEV[SdotScrossS[1, 2, 3], GaussianState[{"SpinRotation", "Charge", "TimeReversal"}]]
```

$$0$$

## On one site

The contractions on one site are the densities, the magnetization and the on-site pairing, and the options fix them:
`"Filling" -> ν` sets $\langle n_{i\sigma}\rangle=\nu$; `"Magnetization" -> m` makes them $\nu\pm m$ (in the general state a
vector `{mx, my, mz}`, and a function of the site, `m[#] &`, for a pattern); `"OnSitePairing" -> η` sets
$\langle c_{i\uparrow}c_{i\downarrow}\rangle=\eta$. A state with a magnetization along $z$, and the correlation of two spins in it:

```mathematica
magnetic = GaussianState[{"SpinRotationZ", "Charge"}, "Filling" -> 1/2, "Magnetization" -> m];
ContractionTable[magnetic, {1}, "Normal"]
```

```math
\left(
\begin{array}{cc}
 m+\frac{1}{2} & 0 \\
 0 & \frac{1}{2}-m \\
\end{array}
\right)
```

```mathematica
WickVEV[SzOp[1] ** SzOp[2], magnetic]
```

$$-\frac{1}{4} {\chi^{\downarrow}_{12}} {\chi^{\downarrow}_{12}}^{\ast}-\frac{1}{4} {\chi^{\uparrow}_{12}} {\chi^{\uparrow}_{12}}^{\ast}+m^2$$

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

A contraction that no rule gives is zero.

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