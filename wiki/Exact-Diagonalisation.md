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
$Z(T,\mu)=\mathrm{Tr}\thinspace e^{-(H-\mu N)/T}$, the trace over the whole Fock space (here sixteen states), with the
temperature $T$ and the chemical potential $\mu$ as variables. `PartitionFunction[H, space, T, μ]` gives it exactly, as a
sum over the levels of every sector:

```mathematica
space = FockSpace[{1, 2}];
Z = PartitionFunction[dimer, space, T, μ]
```

$$e^{\frac{2 \mu +\frac{1}{2} \left(-\sqrt{16 t^2+U^2}-U\right)}{T}}+e^{\frac{2 \mu +\frac{1}{2} \left(\sqrt{16 t^2+U^2}-U\right)}{T}}+2 e^{\frac{\mu +\frac{1}{2} (2 t-U)}{T}}+2 e^{\frac{\mu -t-\frac{U}{2}}{T}}+2 e^{\frac{3 \mu +\frac{1}{2} (2 t-U)}{T}}+2 e^{\frac{3 \mu -t-\frac{U}{2}}{T}}+3 e^{\frac{2 \mu }{T}}+e^{\frac{2 \mu -U}{T}}+e^{\frac{4 \mu -U}{T}}+e^{-\frac{U}{T}}$$

Its logarithm is the grand potential $\Omega(T,\mu)=-T\log Z$, and the physical quantities are its derivatives, taken before
any value is given to $\mu$. The particle number $N=-\partial\Omega/\partial\mu$ is two at $\mu=0$: with the interaction
written as $\frac U2(n_i-1)^2$, $\mu=0$ is half filling, by particle–hole symmetry. Its derivative $\partial N/\partial\mu$
is the compressibility:

```mathematica
Ω = -T Log[Z];
FullSimplify[{-D[Ω, μ], -D[Ω, {μ, 2}]} /. μ -> 0, T > 0]
```

$$\left\\{2,\frac{4 \left(\cosh \left(\frac{t}{T}\right)+e^{-\frac{U}{2 T}}\right)}{T \left(\cosh \left(\frac{\sqrt{16 t^2+U^2}}{2 T}\right)+4 \cosh \left(\frac{t}{T}\right)+3 \cosh \left(\frac{U}{2 T}\right)\right)}\right\\}$$

The entropy $S = -\partial\Omega/\partial T$ at half filling goes from zero (the ground state is a singlet; the limit at
$U = 4t$) to $\log 16$, the sixteen states equally likely:

```mathematica
S = -D[Ω, T];
{Limit[S /. {μ -> 0, U -> 4 t}, T -> 0, Assumptions -> t > 0], Limit[S /. μ -> 0, T -> Infinity]}
```

$$\\{0,\log (16)\\}$$

With numbers, `GrandPotential[H, space, T, μ]` is $-T\log Z$ computed from the lowest level, without overflow at a low
temperature.

### The Green's function

`LehmannGreen[H, space, modes]` is the single-particle Green's function as a sum over the same levels, the Lehmann
representation, $G_{ab}(z)=\sum_{nm}(w_n+w_m)\thinspace \langle n|c_a|m\rangle\langle m|c^\dagger_b|n\rangle/(z-E_m+E_n)$, with $w_n$
the weight of a state and the energies those of $H-\mu N$ (`"ChemicalPotential" -> μ`); `"Temperature" -> T` sets the
temperature. It works in the grand-canonical ensemble on the whole Fock space: no sector is chosen by hand. At zero
temperature only the ground state of $H-\mu N$ has weight; it depends on $\mu$, and `"Assumptions"` say where $\mu$ is,
here inside the gap, where the ground state is the singlet with two electrons and $S^z=0$ (the sector $\\{N,S^z\\}$):

```mathematica
g = LehmannGreen[dimer, space, {c[1, SpinUp], c[2, SpinUp]}, "ChemicalPotential" -> μ,
  "Assumptions" -> t > 0 && U > 0 && Abs[μ] < Sqrt[16 t^2 + U^2]/2 - t];
First[g["GroundSectors"]]
```

$$\\{2,0\\}$$

`GreenMatrix[g, z]` is the matrix $G_{ab}(z)$ between the modes, here those of spin up (spin down gives the same block). In
the basis of the bonding and antibonding orbitals, $(c_1\pm c_2)/\sqrt2$, it is diagonal:

```mathematica
P = {{1, 1}, {1, -1}}/Sqrt[2];
MatrixForm[Simplify[P . GreenMatrix[g, z] . P, t > 0 && U > 0]]
```

```math
\left(
\begin{array}{cc}
 \frac{12 t-4 (\mu +z)}{12 t^2+8 t (\mu +z)+U^2-4 (\mu +z)^2} & 0 \\
 0 & -\frac{4 (\mu +3 t+z)}{12 t^2-8 t (\mu +z)+U^2-4 (\mu +z)^2} \\
\end{array}
\right)
```

The frequency enters as $z+\mu$: the poles, at $\pm(\frac{1}{2}\sqrt{16t^2+U^2} \pm t)-\mu$, are the energies of adding and
of removing an electron, measured from $\mu$. At a temperature every state of the Fock space has the weight
$e^{-(E_n-\mu N_n)/T}/Z$, and the transitions between excited states add poles. The matrix is still diagonal, each orbital
with $G(z)=\frac1Z\sum_p B_p/(z-\varepsilon_p)$, $B_p$ a sum of Boltzmann factors. The four poles of the bonding orbital,
with their $B_p$:

```mathematica
gT = LehmannGreen[dimer, space, {c[1, SpinUp], c[2, SpinUp]}, "Temperature" -> T, "ChemicalPotential" -> μ];
bonding = Select[Transpose[{gT["Poles"], Z (P . # . P)[[1, 1]] & /@ gT["Residues"]}], ! PossibleZeroQ[Last[#]] &];
MatrixForm[FullSimplify[bonding, T > 0 && t > 0 && U > 0]]
```

```math
\left(
\begin{array}{cc}
 -\mu -t-\frac{U}{2} & \frac{1}{2} e^{-\frac{t+U}{T}} \left(e^{\frac{2 \mu +t}{T}}+2 e^{t/T}+3 e^{\frac{2 \mu +U}{2 T}}\right) \left(e^{\frac{2 \mu +2 t+U}{2 T}}+1\right) \\
 -\mu -\frac{1}{2} \sqrt{16 t^2+U^2}+t & \frac{\left(\sqrt{16 t^2+U^2}+4 t\right) e^{-\frac{U-4 \mu }{2 T}} \left(\cosh \left(\frac{\sqrt{16 t^2+U^2}}{2 T}\right)+\cosh \left(\frac{t-\mu }{T}\right)\right)}{\sqrt{16 t^2+U^2}} \\
 -\mu +\frac{1}{2} \sqrt{16 t^2+U^2}+t & \frac{\left(\sqrt{16 t^2+U^2}-4 t\right) e^{-\frac{U-4 \mu }{2 T}} \left(\cosh \left(\frac{\sqrt{16 t^2+U^2}}{2 T}\right)+\cosh \left(\frac{t-\mu }{T}\right)\right)}{\sqrt{16 t^2+U^2}} \\
 \frac{1}{2} (-2 \mu -2 t+U) & \frac{1}{2} e^{-\frac{-\mu +t+U}{T}} \left(e^{\frac{\mu +t}{T}}+e^{\frac{U}{2 T}}\right) \left(3 e^{\frac{2 \mu +2 t+U}{2 T}}+2 e^{\frac{2 \mu }{T}}+1\right) \\
\end{array}
\right)
```

The antibonding orbital has the other four, those of the bonding orbital with $t\to-t$ (changing the sign of $c_2$ maps one
orbital onto the other).

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