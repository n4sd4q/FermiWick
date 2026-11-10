``FermiWick`Gutzwiller` `` handles projected operators: sandwiches $P_i\cdots X\cdots P_i$ of an operator between local
operators of your choice, and their expectation values in the Gutzwiller approximation, to first order in the
contractions between different sites. The sandwich is cheap: every product is split site by site, so each projector
meets only the part of $X$ on its own site, and each site is reduced once in its local algebra of four states.
`GutzwillerProjector` is one example; any local operator will do ([below](#projectors-of-your-own)). Every block on
this page runs after

```mathematica
Needs["FermiWick`"]; Needs["FermiWick`Gutzwiller`"]
```

## The projector

`GutzwillerProjector[i]` is the variational Gutzwiller operator at half filling,
$P_i=\sqrt2\thinspace [\cos\theta\thinspace P_i^{(1)}+\sin\theta\thinspace (P_i^{(0)}+D_i)]$, with $P^{(1)}$ and $P^{(0)}$ the projectors on single
occupation and on the empty site: $\theta=0$ projects out the empty and doubly occupied sites, $\theta=\pi/4$ is no
projection at all. `LocalSandwich[P, {i, j, ...}, X]` is $P_iP_j\cdots X\cdots P_jP_i$.

## Expectation values

`GutzwillerVEV[expr, state]` is `WickVEV` to first order in the contractions between different sites. In a state with
real hopping at half filling, the double occupancy of a projected site:

```mathematica
real = GaussianState[{"SpinRotation", "Charge", "TimeReversal"}, "Filling" -> 1/2];
d = GutzwillerVEV[LocalSandwich[GutzwillerProjector, {1}, DOp[1]], real]
```

$$\frac{\sin ^2(\theta )}{2}$$

and the hopping between two projected sites, reduced by a factor:

```mathematica
hop = Sum[PlusHC[cDag[1, s] ** c[2, s]], {s, $Spins}];
q = Simplify[GutzwillerVEV[LocalSandwich[GutzwillerProjector, {1, 2}, hop], real]/WickVEV[hop, real]]
```

$$\sin ^2(2 \theta )$$

that is $q=8d(1-2d)$.

## The Brinkman–Rice transition

The energy per site of the half-filled Hubbard model is then $q\thinspace \varepsilon_0+U d$, with $\varepsilon_0$ the (negative) kinetic
energy of the uncorrelated state. Written with $x=\cos2\theta$, its minimum:

```mathematica
energy = Simplify[TrigReduce[q Subscript[ε, 0] + U d] /. θ -> ArcCos[x]/2, -1 <= x <= 1];
Solve[D[energy, x] == 0, x]
```

$$\left\\{\left\\{x\to -\frac{U}{8 \varepsilon _0}\right\\}\right\\}$$

so $d=\frac14(1-U/U_c)$ with $U_c=8|\varepsilon_0|$: the double occupancy vanishes at $U_c$, the Brinkman–Rice
transition. The [showcase](https://n4sd4q.github.io/FermiWick/Showcase/FermiWick_Showcase.pdf) goes through it.

## Projectors of your own

A projector is any function of the site that gives a local operator. The general Gutzwiller operator weighs the empty,
singly and doubly occupied states with their own factors:

```mathematica
P[i_] := Subscript[λ, 0] EmptyProjector[i] + Subscript[λ, 1] SingleOccupationProjector[i] + Subscript[λ, 2] DOp[i];
GutzwillerVEV[LocalSandwich[P, {1}, 1], real]
```

$$\frac{1}{4} \left(\lambda _0^2+2 \lambda _1^2+\lambda _2^2\right)$$

is its normalization in the state above, and the hopping between two sites is reduced by

```mathematica
Simplify[GutzwillerVEV[LocalSandwich[P, {1, 2}, hop], real]/WickVEV[hop, real]]
```

$$\frac{1}{4} \lambda _1^2 \left(\lambda _0+\lambda _2\right){}^2$$

Larger sandwiches stay cheap. The chirality between three projected sites, a product of six fermions between six
projectors, takes under a second, and keeps only the singly occupied part of each site; expanding the same product into
ladder operators before normal-ordering it takes minutes:

```mathematica
ToSpinChargeBasis[LocalSandwich[P, {1, 2, 3}, SdotScrossS[1, 2, 3]]]
```

$$\lambda _1^6 \mathbf S_{1}\cdot(\mathbf S_{2}\times\mathbf S_{3})$$

`LocalSandwich[PL, sites, X, PR]` takes different operators on the left and on the right.
