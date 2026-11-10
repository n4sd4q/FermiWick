``FermiWick`Fock` `` writes any expression of the package as a sparse matrix, on the Fock space of a few sites or on a
sector with fixed quantum numbers: a model written as on paper is solved exactly, and a symbolic result can be checked.
Every block on this page runs after

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

## The Hubbard dimer

`FockMatrix[expr, sector]` is the matrix of an expression, symbolic coefficients included:

```mathematica
dimer = -t Sum[PlusHC[cDag[1, s] ** c[2, s]], {s, $Spins}] + U Sum[DOp[i], {i, 2}];
sector = FockSector[FockSpace[{1, 2}], "N" -> 2, "Sz" -> 0];
Eigenvalues[Normal[FockMatrix[dimer, sector]]]
```

$$\left\\{0,U,\frac{1}{2} \left(U-\sqrt{16 t^2+U^2}\right),\frac{1}{2} \left(\sqrt{16 t^2+U^2}+U\right)\right\\}$$

With numbers, `LowestStates[H, k]` gives the `k` lowest energies and states, and `StateToKet` writes a state with the
occupations of the sites:

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

## The Green's function

`LehmannGreen[H, space, modes]` is the single-particle Green's function at zero temperature, in the Lehmann
representation; its poles are the excitation energies. For the dimer at half filling:

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