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

`FockMatrix[expr, sector]` is the matrix of an expression, symbolic coefficients included. The dimer, with the
interaction in the form symmetric about half filling, $\frac U2\sum_i(n_i-1)^2$:

```mathematica
dimer = -t Sum[PlusHC[cDag[1, s] ** c[2, s]], {s, $Spins}] + U/2 Sum[(nTotOp[i] - 1)^2, {i, 2}];
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

Everything at a finite temperature follows from the grand-canonical partition function,
$Z=\mathrm{Tr}\thinspace e^{-(H-\mu N)/T}$, the trace over the whole Fock space, here sixteen states.
`PartitionFunction[H, space, T, μ]` gives it exactly, as a sum of exponentials over the levels of every sector:

```mathematica
Z = PartitionFunction[dimer, FockSpace[{1, 2}], T, μ]
```

$$e^{\frac{2 \mu +\frac{1}{2} \left(-\sqrt{16 t^2+U^2}-U\right)}{T}}+e^{\frac{2 \mu +\frac{1}{2} \left(\sqrt{16 t^2+U^2}-U\right)}{T}}+2 e^{\frac{\mu +\frac{1}{2} (2 t-U)}{T}}+2 e^{\frac{\mu -t-\frac{U}{2}}{T}}+2 e^{\frac{3 \mu +\frac{1}{2} (2 t-U)}{T}}+2 e^{\frac{3 \mu -t-\frac{U}{2}}{T}}+3 e^{\frac{2 \mu }{T}}+e^{\frac{2 \mu -U}{T}}+e^{\frac{4 \mu -U}{T}}+e^{-\frac{U}{T}}$$

Its logarithm is the grand potential $\Omega=-T\log Z$, and the physical quantities are its derivatives. The particle
number $-\partial\Omega/\partial\mu$ is two at $\mu=0$, half filling, by the symmetry of the interaction; its derivative is
the compressibility:

```mathematica
Ω = -T Log[Z];
FullSimplify[{-D[Ω, μ], -D[Ω, {μ, 2}]} /. μ -> 0, T > 0]
```

$$\left\\{2,\frac{4 \left(\cosh \left(\frac{t}{T}\right)+e^{-\frac{U}{2 T}}\right)}{T \left(\cosh \left(\frac{\sqrt{16 t^2+U^2}}{2 T}\right)+4 \cosh \left(\frac{t}{T}\right)+3 \cosh \left(\frac{U}{2 T}\right)\right)}\right\\}$$

The free energy at half filling,

```mathematica
F = FullSimplify[Ω /. μ -> 0, T > 0]
```

$$-T \log \left(2 e^{-\frac{U}{2 T}} \left(\cosh \left(\frac{\sqrt{16 t^2+U^2}}{2 T}\right)+4 \cosh \left(\frac{t}{T}\right)+3 \cosh \left(\frac{U}{2 T}\right)\right)\right)$$

and the entropy $S = -\partial F/\partial T$, from zero (the ground state is a singlet; the limit at $U = 4t$) to
$\log 16$, the sixteen states equally likely:

```mathematica
S = -D[F, T];
{Limit[S /. U -> 4 t, T -> 0, Assumptions -> t > 0], Limit[S, T -> Infinity]}
```

$$\\{0,\log (16)\\}$$

With numbers, `GrandPotential[H, space, T, μ]` is $-T\log Z$ computed from the lowest level, without overflow at a low
temperature.

### The Green's function

`LehmannGreen[H, space, modes]` is the single-particle Green's function as a sum over the same levels, the Lehmann
representation, $G_{ab}(z)=\sum_{nm}(w_n+w_m)\thinspace \langle n|c_a|m\rangle\langle m|c^\dagger_b|n\rangle/(z-E_m+E_n)$, with $w_n$ the
weight of a state; it returns the poles and the residues, and `GreenMatrix` evaluates it. At zero temperature, in
closed form:

```mathematica
g = LehmannGreen[dimer, FockSpace[{1, 2}], {c[1, SpinUp]}];
FullSimplify[GreenMatrix[g, z][[1, 1]], t > 0 && U > 0]
```

$$-\frac{4 z \left(36 t^2+U^2-4 z^2\right)}{-8 z^2 \left(20 t^2+U^2\right)+\left(12 t^2+U^2\right)^2+16 z^4}$$

Its four poles, at $\pm(\frac{1}{2}\sqrt{16t^2+U^2} \pm t)$, are the energies of adding and of removing an electron. At
a temperature, `"Temperature" -> T` gives every state the weight $e^{-E_n/T}/Z$, and the transitions between excited
states add four poles:

```mathematica
gT = LehmannGreen[dimer, FockSpace[{1, 2}], {c[1, SpinUp]}, "Temperature" -> T];
Complement[gT["Poles"], g["Poles"]]
```

$$\left\\{-t-\frac{U}{2},t-\frac{U}{2},t+\frac{U}{2},\frac{1}{2} (U-2 t)\right\\}$$

Each of the new poles has the same weight, which vanishes as $T\to0$:

```mathematica
weights = AssociationThread[gT["Poles"], gT["Residues"][[All, 1, 1]]];
FullSimplify[weights[t + U/2], T > 0]
```

$$\frac{3 \left(\cosh \left(\frac{t}{T}\right)+\cosh \left(\frac{U}{2 T}\right)\right)}{4 \left(\cosh \left(\frac{\sqrt{16 t^2+U^2}}{2 T}\right)+4 \cosh \left(\frac{t}{T}\right)+3 \cosh \left(\frac{U}{2 T}\right)\right)}$$

With symbols at zero temperature, the ground level must be the lowest for every value that `"Assumptions"` allows
(none are needed here); `"ChemicalPotential" -> μ` uses $H-\mu N$.

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

## Thermodynamics and the Green's function, with numbers

With numbers the same functions work on larger clusters; at a temperature every sector is diagonalised in full. The
ring of four sites above (256 states) at $T = t/2$, with the chemical potential $U/2 = 2$ of half filling for its form
of the interaction: the grand potential, and the number of poles of the local Green's function:

```mathematica
g4 = LehmannGreen[ring[4], FockSpace[Range[4]], {c[1, SpinUp]}, "Temperature" -> 1/2, "ChemicalPotential" -> 2];
{g4["GrandPotential"], GrandPotential[ring[4], FockSpace[Range[4]], 1/2, 2], Length[g4["Poles"]]}
```

```text
{-10.8025, -10.8025, 165}
```

At zero temperature only the ground multiplet has weight, and the poles are the excitation energies:

```mathematica
Sort[LehmannGreen[dimer /. {t -> 1, U -> 4}, FockSpace[{1, 2}], {c[1, SpinUp]}]["Poles"]]
```

```text
{-3.82843, -1.82843, 1.82843, 3.82843}
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