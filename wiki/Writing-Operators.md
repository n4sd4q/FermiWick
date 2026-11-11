Operators are written as on paper: ladder operators, multiplied and summed. The package keeps them compact and reads
them back as local operators and spins. Every block on this page runs after

```mathematica
Needs["FermiWick`"]
```

## Ladder operators

`c[i, s]` annihilates a fermion on site `i` with spin `s`, and `cDag[i, s]` creates it. The spin is `SpinUp`, `SpinDown`
(shown as ↑ and ↓) or a symbol; the site is a number or a symbol. `**` multiplies in the order written, and `NormalOrder`
brings every creator to the left:

```mathematica
NormalOrder[c[1, SpinUp] ** cDag[1, SpinUp]]
```

$$1-c^\dagger_{1\uparrow}c_{1\uparrow}$$

`a · b` (typed `\[CenterDot]`, or `Esc . Esc`) is the product already in normal order.

## Sums, and sites as symbols

Sums are ordinary sums; `$Spins` is the list of the two spins. Different symbols are different sites, so a result
written with symbols holds on any lattice:

```mathematica
NormalOrder[Sum[c[i, s] ** cDag[j, s], {s, $Spins}]]
```

$$-c^\dagger_{j\downarrow}c_{i\downarrow}-c^\dagger_{j\uparrow}c_{i\uparrow}$$

To include the case `i = j`, write it apart, with the same symbol for both.

## Local operators

Products on one site are local operators, which stay compact: `nOp[i, s]` and `nTotOp[i]` (densities), `SzOp[i]`,
`SpOp[i]`, `SmOp[i]` (spin), `DOp[i]` (double occupancy), `B[i]` and `Bdag[i]` (doublons), `EmptyProjector[i]` and
`SingleOccupationProjector[i]`. `ToSpinChargeBasis` writes any expression with them:

```mathematica
ToSpinChargeBasis[SpOp[1] ** SmOp[1]]
```

$$-D_{1}+\frac{n_{1}}{2}+S^z_{1}$$

and `ToLadderBasis` writes them back with ladder operators:

```mathematica
ToLadderBasis[DOp[1]]
```

$$c^\dagger_{1\uparrow}c^\dagger_{1\downarrow}c_{1\downarrow}c_{1\uparrow}$$

## Spins on several sites

`SdotS[i, j]` is $\mathbf S_i\cdot\mathbf S_j$ and `SdotScrossS[i, j, k]` the scalar chirality
$\mathbf S_i\cdot(\mathbf S_j\times\mathbf S_k)$. `ToSpinChargeBasis` recognizes them in a result, also from their
components:

```mathematica
ToSpinChargeBasis[SpOp[1] ** SmOp[2] + SmOp[1] ** SpOp[2] + 2 SzOp[1] ** SzOp[2]]
```

$$2 \mathbf S_{1}\cdot\mathbf S_{2}$$

## The Hermitian conjugate

`Dagger` reverses products, exchanges creators and annihilators and conjugates the coefficients: the amplitudes χ, Δ
are complex, other symbols are taken as real (the option `"Complex" -> {...}` lists more complex symbols). `PlusHC`
adds the "+ h.c.", and `HermitianOperatorQ` checks a Hamiltonian:

```mathematica
hubbard = -t Sum[PlusHC[cDag[1, s] ** c[2, s]], {s, $Spins}] + U (DOp[1] + DOp[2]);
HermitianOperatorQ[hubbard]
```

```text
True
```

A complex hopping is declared complex, and its conjugate appears in the "+ h.c.":

```mathematica
PlusHC[τ cDag[1, SpinUp] ** c[2, SpinUp], "Complex" -> {τ}]
```

$$\tau ^{\ast} c^\dagger_{2\uparrow}c_{1\uparrow}+\tau  c^\dagger_{1\uparrow}c_{2\uparrow}$$

## Coefficients

`OperatorCoefficients` gives the coefficient of every monomial, in normal order; the constant has the key `1`:

```mathematica
OperatorCoefficients[ToLadderBasis[hubbard], cDag[2, SpinDown] ** c[1, SpinDown]]
```

$$-t$$

## Modes with any labels

`c[l1, l2, ...]` can carry any labels: orbitals, momenta, labels of your own. `DeclareModeLabels` says what each one is:
`"Site"` (used for locality), `"Spin"` (a spin ½), or any other name, compared with `KroneckerDelta`:

```mathematica
DeclareModeLabels[{"Site", "Orbital", "Spin"}];
NormalOrder[c[1, a, SpinUp] ** cDag[1, b, SpinUp]]
```

$$\delta _{a,b}-c^\dagger_{1,b,\uparrow}c_{1,a,\uparrow}$$

The defaults are `{"Site"}` for one label and `{"Site", "Spin"}` for two.

## Operators of your own

`DefineOperator[head, n, definition, display]` defines an operator on `n` sites from known ones; it stays compact, and
`ToSpinChargeBasis` recognizes it in a result. Here the operator $b_{ij}$ that annihilates a singlet on two sites, and
the number of singlets on a bond, $b^\dagger_{12}b_{12}$:

```mathematica
DefineOperator[pair, 2, Function[{i, j}, (c[i, SpinUp] ** c[j, SpinDown] - c[i, SpinDown] ** c[j, SpinUp])/Sqrt[2]],
  Function[{i, j}, Subscript[b, Row[{i, j}]]]];
ToSpinChargeBasis[Dagger[pair[1, 2]] ** pair[1, 2]]
```

$$\frac{1}{4} n_{1}n_{2}-\mathbf S_{1}\cdot\mathbf S_{2}$$

On two singly occupied sites ($n_1=n_2=1$) this is $\frac14-\mathbf S_1\cdot\mathbf S_2$, the
projector on the singlet.

The adjoint of a defined operator stays compact too: `Dagger` gives the operator itself when it is Hermitian, the
operator with its sites permuted when that is its adjoint (a hopping), and otherwise `AdjointOperator[head][sites]`,
shown with a dagger, expanded and recognized like the operator itself. Here $b^\dagger_{12}$ and its expectation value
in a state with the pairing Δ:

```mathematica
{Dagger[pair[1, 2]], WickVEV[Dagger[pair[1, 2]], SpinSingletState["Pairing" -> Δ]]}
```

$$\left\\{b_{12}{}^{\dagger },\sqrt{2} \Delta_{12}^{\ast}\right\\}$$