`MeanField[expr, state]` is the Hartree–Fock–Bogoliubov Hamiltonian of an operator in a Gaussian state: in every
normal-ordered product two operators are kept, in every way, and the rest is replaced by its expectation value, with
the constant that keeps the energy right. The result is quadratic, and has the same expectation value and the same
first variation as the operator. Every block on this page runs after

```mathematica
Needs["FermiWick`"]
```

## The mean field of an interaction

The Hubbard interaction of a site, in a state with magnetization $m$ along z (see
[Declaring a state](Declaring-a-State)):

```mathematica
magnetic = GaussianState[{"SpinRotationZ", "Charge"}, "Filling" -> 1/2, "Magnetization" -> m];
ToSpinChargeBasis[MeanField[U DOp[1], magnetic]]
```

$$m^2 U-2 m S^z_{1} U+\frac{1}{2} n_{1} U-\frac{U}{4}$$

and in a state with on-site pairing η, the Hartree term and the pairing fields:

```mathematica
pairing = GaussianState["SpinRotation", "Filling" -> ν, "OnSitePairing" -> η];
Collect[ToSpinChargeBasis[MeanField[U DOp[1], pairing]], U]
```

$$U \left(-B_{1} \eta ^{\ast}-B^\dagger_{1} \eta -\eta  \eta ^{\ast}-\nu ^2+\nu  n_{1}\right)$$

## The fields

The fields are the coefficients of the result (`OperatorCoefficients`). The exchange of a bond, in the state invariant
under spin rotations, gives a hopping field proportional to the bond's own amplitude:

```mathematica
fields = OperatorCoefficients[MeanField[J SdotS[i, j], SpinSingletState[]]];
fields[cDag[i, SpinUp] ** c[j, SpinUp]]
```

$$-\frac{3}{4} J \chi_{ij}^{\ast}$$

## Any order

An operator of six fermions is no harder. The hopping field that the chirality $J\thinspace \mathbf S_i\cdot(\mathbf S_j\times\mathbf S_k)$
produces on the bond $i\to j$ is the product of the amplitudes along the other two sides of the triangle:

```mathematica
fields = OperatorCoefficients[MeanField[J SdotScrossS[i, j, k], SpinSingletState[]]];
fields[cDag[i, SpinUp] ** c[j, SpinUp]]
```

$$-\frac{3}{4} i J \chi_{jk} \chi_{ik}^{\ast}$$

## The variational property

The mean-field Hamiltonian has the expectation value of the operator, in the same state:

```mathematica
Simplify[WickVEV[MeanField[SdotS[1, 2], SpinSingletState[]], SpinSingletState[]] - WickVEV[SdotS[1, 2], SpinSingletState[]]]
```

$$0$$

## A self-consistent loop

The antiferromagnet of the Hubbard model on a ring of six sites, at half filling, with $U=4t$: a state with the
magnetization $\pm m$ alternating from site to site (declared by a rule), its mean-field Hamiltonian, and the
single-particle matrix of each spin, read from the coefficients.

```mathematica
L = 6;
neel = GaussianState[{{cDag[i_, s_], c[j_, r_]} /; i === j :> KroneckerDelta[s, r] (1/2 + (-1)^i If[s === SpinUp, m, -m])}];
hMF = -Sum[PlusHC[cDag[i, s] ** c[Mod[i, L] + 1, s]], {i, L}, {s, $Spins}] + 4 MeanField[Sum[DOp[i], {i, L}], neel];
h[s_] := h[s] = Table[OperatorCoefficients[hMF, cDag[a, s] ** c[b, s]], {a, L}, {b, L}];
h[SpinUp] // MatrixForm
```

```math
\left(
\begin{array}{cccccc}
 4 m+2 & -1 & 0 & 0 & 0 & -1 \\
 -1 & 2-4 m & -1 & 0 & 0 & 0 \\
 0 & -1 & 4 m+2 & -1 & 0 & 0 \\
 0 & 0 & -1 & 2-4 m & -1 & 0 \\
 0 & 0 & 0 & -1 & 4 m+2 & -1 \\
 -1 & 0 & 0 & 0 & -1 & 2-4 m \\
\end{array}
\right)
```

Each step fills the three lowest levels of each spin and reads the new $m$ from the density of site 2:

```mathematica
step[m0_] := Module[{e, v}, {e, v} = Eigensystem[N[h[SpinUp] /. m -> m0]];
   Total[Abs[v[[Ordering[e, 3]]]]^2][[2]] - 1/2];
FixedPoint[step, 0.3, 500, SameTest -> (Abs[#1 - #2] < 10^-12 &)]
```

```text
0.378884
```