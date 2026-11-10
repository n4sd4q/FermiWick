``FermiWick`Gutzwiller` `` handles Gutzwiller-projected operators and their expectation values in the Gutzwiller
approximation, to first order in the contractions between different sites. Every block on this page runs after

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