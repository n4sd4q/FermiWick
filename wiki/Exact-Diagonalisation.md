``FermiWick`Fock` `` writes any expression of the package as a sparse matrix, on the Fock space of a few sites or on a
sector with fixed quantum numbers: a model written as on paper is solved exactly, in closed form when it is small, and a
symbolic result can be checked. Every block on this page runs after

```mathematica
Needs["FermiWick`"]; Needs["FermiWick`Fock`"]
```

## Spaces and sectors

`FockSpace[sites]` has the modes `c[i, s]` of the sites, numbered site by site with the spin fastest;
`FockSector[space, "N" -> n, "Sz" -> sz]` is the part with `n` particles and spin `sz`:

```mathematica
space = FockSpace[Range[4]];
{space["Dimension"], FockSector[space, "N" -> 4, "Sz" -> 0]["Dimension"]}
```

```text
{256, 36}
```

`"Counts" -> {{modes, n}, ...}` fixes the number of particles in any set of modes. `FockSpace[sites, "Orbitals" -> {...}]`
has the modes `c[i, o, s]`, and `"Spins" -> None` gives spinless modes `c[i]`.

## The Hubbard dimer, analytically

`FockMatrix[expr, sector]` is the matrix of an expression, symbolic coefficients included:

```mathematica
dimer = -t Sum[PlusHC[cDag[1, s] ** c[2, s]], {s, $Spins}] + U Sum[DOp[i], {i, 2}];
sector = FockSector[FockSpace[{1, 2}], "N" -> 2, "Sz" -> 0];
Eigenvalues[Normal[FockMatrix[dimer, sector]]]
```

$$\left\\{0,U,\frac{1}{2} \left(U-\sqrt{16 t^2+U^2}\right),\frac{1}{2} \left(\sqrt{16 t^2+U^2}+U\right)\right\\}$$

The eigenvectors are exact too, and `StateToKet` writes a state with the occupations of the sites. The ground state is
the lowest of the four energies for every $t, U \gt 0$; their values at one point, with `N`, only rank them (`Ordering`
of the exact values would sort them in the canonical order of expressions, not by size):

```mathematica
{vals, vecs} = Eigensystem[Normal[FockMatrix[dimer, sector]]];
ground = Normalize[vecs[[First[Ordering[N[vals /. {t -> 1, U -> 1}], 1]]]]];
FullSimplify[StateToKet[sector, ground], t > 0 && U > 0]
```

```math
\frac{4 t (|\{0\},\{\uparrow\downarrow\}\rangle +|\{\uparrow\downarrow\},\{0\}\rangle )-\left(\sqrt{16 t^2+U^2}+U\right) (|\{\downarrow\},\{\uparrow\}\rangle -|\{\uparrow\},\{\downarrow\}\rangle )}{2 \sqrt{U \left(\sqrt{16 t^2+U^2}+U\right)+16 t^2}}
```

A singlet, whose weight on the doubly occupied sites falls as $U$ grows.

### Thermodynamics

The whole Fock space, sixteen states, at the chemical potential $U/2$ of half filling: the free energy
$F = -T \log \mathrm{Tr}\thinspace e^{-H/T}$,

```mathematica
h = dimer - U/2 Sum[nTotOp[i], {i, 2}];
levels = Eigenvalues[Normal[FockMatrix[h, FockSpace[{1, 2}]]]];
F = FullSimplify[-T Log[Total[Exp[-levels/T]]], T > 0]
```

$$-T \log \left(2 e^{\frac{U}{2 T}} \left(\cosh \left(\frac{\sqrt{16 t^2+U^2}}{2 T}\right)+4 \cosh \left(\frac{t}{T}\right)+3 \cosh \left(\frac{U}{2 T}\right)\right)\right)$$

and the entropy $S = -\partial F/\partial T$, from zero (the ground state is a singlet; the limit at $U = 4t$) to
$\log 16$, the sixteen states equally likely:

```mathematica
S = -D[F, T];
{Limit[S /. U -> 4 t, T -> 0, Assumptions -> t > 0], Limit[S, T -> Infinity]}
```

$$\\{0,\log (16)\\}$$

### The Green's function

`LehmannGreen` takes numbers (below). Analytically, the Lehmann sum is written with the matrices of `c` and `cDag`
between sectors, `FockMatrix[op, {fromSector, toSector}]`, and the eigenstates of the sectors with three and with one
particle:

```mathematica
e0 = ground . Normal[FockMatrix[h, sector]] . ground;
part[op_, to_, sign_] := Module[{e, v, w},
  {e, v} = Eigensystem[Normal[FockMatrix[h, to]]];
  w = Normal[FockMatrix[op, {sector, to}]] . ground;
  Sum[(Normalize[v[[n]]] . w)^2/(z - sign (e[[n]] - e0)), {n, Length[e]}]];
FullSimplify[part[cDag[1, SpinUp], FockSector[FockSpace[{1, 2}], "N" -> 3, "Sz" -> 1/2], 1] +
  part[c[1, SpinUp], FockSector[FockSpace[{1, 2}], "N" -> 1, "Sz" -> -1/2], -1], t > 0 && U > 0]
```

$$-\frac{4 z \left(36 t^2+U^2-4 z^2\right)}{-8 z^2 \left(20 t^2+U^2\right)+\left(12 t^2+U^2\right)^2+16 z^4}$$

Its four poles, at $\pm(\frac{1}{2}\sqrt{16t^2+U^2} \pm t)$, are the energies of adding and of removing an electron.

## With numbers

`LowestStates[H, k]` gives the `k` lowest energies and states of a matrix with numbers, dense or (on large sectors) by
Arnoldi, and `StateToKet` writes them:

```mathematica
{energies, states} = LowestStates[FockMatrix[dimer /. {t -> 1, U -> 4}, sector], 1];
Chop[StateToKet[sector, First[states]]]
```

```text
-0.270598 Ket[{{0}, {↑↓}}] + 0.653281 Ket[{{↓}, {↑}}] - 0.653281 Ket[{{↑}, {↓}}] - 0.270598 Ket[{{↑↓}, {0}}]
```

## A Hubbard ring

Ten sites at half filling, a sector of 63 504 states; the spin correlations from site 1 alternate in sign:

```mathematica
ring[L_] := -Sum[PlusHC[cDag[i, s] ** c[Mod[i, L] + 1, s]], {i, L}, {s, $Spins}] + 4 Sum[DOp[i], {i, L}];
ringSector = FockSector[FockSpace[Range[10]], "N" -> 10, "Sz" -> 0];
{energies, states} = LowestStates[FockMatrix[ring[10], ringSector], 1];
Table[FockExpectation[SdotS[1, 1 + d], First[states], ringSector], {d, 5}]
```

```text
{-0.322566, 0.103488, -0.107313, 0.0724331, -0.0859561}
```

A ring of twelve sites (853 776 states) is written as a matrix in about four seconds; above `$FockDenseLimit` the
lowest states are found by Arnoldi, with degenerate multiplets complete.

## A spin model, in its own space

With one fermion on every site (`"Counts"`), the sector is the space of the spins, and `SdotS` is the Heisenberg
exchange. The ring of ten spins, whose energy per site tends to $\frac14-\ln 2$ for long rings:

```mathematica
spins = FockSector[FockSpace[Range[10]], "Sz" -> 0, "Counts" -> Table[{{c[i, SpinUp], c[i, SpinDown]}, 1}, {i, 10}]];
{spins["Dimension"], First[First[LowestStates[FockMatrix[Sum[SdotS[i, Mod[i, 10] + 1], {i, 10}], spins], 1]]]}
```

```text
{252, -4.51545}
```

## The Green's function, with numbers

`LehmannGreen[H, space, modes]` is the single-particle Green's function at zero temperature, in the Lehmann
representation, of a Hamiltonian with numbers; its poles are the excitation energies. For the dimer at half filling:

```mathematica
g = LehmannGreen[FockMatrix[dimer /. {t -> 1, U -> 4}, FockSpace[{1, 2}]], FockSpace[{1, 2}], {c[1, SpinUp]}, "N" -> {2}];
Sort[g["Poles"]]
```

```text
{-1.82843, 0.171573, 3.82843, 5.82843}
```

## Checking a symbolic result

The ground state of free fermions is Gaussian, so Wick's theorem must be exact there: the hopping amplitude of the free
ring of six sites, read from its exact ground state, and the exchange of a bond by Wick's theorem with it, against the
exact value.

```mathematica
freeSector = FockSector[FockSpace[Range[6]], "N" -> 6, "Sz" -> 0];
free = First[Last[LowestStates[FockMatrix[ring[6] /. DOp[_] -> 0, freeSector], 1]]];
chi = FockExpectation[cDag[1, SpinUp] ** c[2, SpinUp], free, freeSector];
{WickVEV[SdotS[1, 2], SpinSingletState[]] /. χ[1, 2] -> chi, FockExpectation[SdotS[1, 2], free, freeSector]}
```

```text
{-0.166667, -0.166667}
```