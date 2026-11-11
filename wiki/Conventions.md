The conventions of FermiWick, in one place. The [manual](https://n4sd4q.github.io/FermiWick/Documentation/FermiWick_Manual.pdf)
gives them with every function.

## Operators

- `c[i, s]` annihilates a fermion on site `i` with spin `s`, `cDag[i, s]` creates it; $\\{c_{i\sigma},c^\dagger_{j\sigma'}\\}=\delta_{ij}\delta_{\sigma\sigma'}$.
  The spin is `SpinUp` (↑), `SpinDown` (↓) or a symbol. Different symbolic sites are different sites.
- `a ** b` is the product in the order written; `a · b` (`\[CenterDot]`) the product in normal order, every `cDag` to the
  left of every `c`.
- Local operators: $n_{i\sigma}=c^\dagger_{i\sigma}c_{i\sigma}$, $n_i=n_{i\uparrow}+n_{i\downarrow}$,
  $\mathbf S_i=\frac12\sum_{\sigma\sigma'}c^\dagger_{i\sigma}\boldsymbol\sigma_{\sigma\sigma'}c_{i\sigma'}$,
  $S_i^+=c^\dagger_{i\uparrow}c_{i\downarrow}$, $D_i=n_{i\uparrow}n_{i\downarrow}$, $B_i=c_{i\downarrow}c_{i\uparrow}$,
  $B^\dagger_i=c^\dagger_{i\uparrow}c^\dagger_{i\downarrow}$.
- `Dagger` takes the amplitudes χ, Δ as complex and every other symbol as real, unless listed with `"Complex"`.

## Contractions

A Gaussian state has one amplitude for each contraction: the hopping
$\chi^{\sigma\sigma'}_{ij}=\langle c^\dagger_{i\sigma}c_{j\sigma'}\rangle$ (`χ[i, j, s, t]`) and the pairing
$\Delta^{\sigma\sigma'}_{ij}=\langle c_{i\sigma}c_{j\sigma'}\rangle$ (`Δ[i, j, s, t]`). A symmetry removes indices: with $S^z$ conserved
$\chi^\sigma_{ij}$ (`χ[i, j, s]`) and $\Delta^{\uparrow\downarrow}_{ij}$, $\Delta^{\downarrow\uparrow}_{ij}$; with SU(2) $\chi_{ij}$ and the
singlet $\Delta_{ij}=\langle c_{i\uparrow}c_{j\downarrow}\rangle$, symmetric in the sites. Hermiticity gives
$\langle c^\dagger_bc_a\rangle=\langle c^\dagger_ac_b\rangle^{\ast}$ and antisymmetry $\langle c_bc_a\rangle=-\langle c_ac_b\rangle$: an amplitude
is written with its modes in canonical order (the sites, then the spins, ↑ first), the other order as its conjugate or its
opposite. On a site, $\chi^{\sigma\sigma}_{ii}$ are the densities (`"Filling"`, `"Magnetization"`),
$\chi^{\uparrow\downarrow}_{ii}=\langle S^+_i\rangle$ and $\Delta^{\uparrow\downarrow}_{ii}=\langle c_{i\uparrow}c_{i\downarrow}\rangle$ (`"OnSitePairing"`); with
SU(2), $\chi_{ii}$ is the density per spin. Modes with more labels group them by type, the spins last:
`χ[i, j, o, o', s, t]`. Time reversal, $T=i\sigma^yK$, makes the components of opposite spins conjugates of each other:
$\chi^\downarrow_{ij}=\chi^{\uparrow\ast}_{ij}$, $\chi^{\downarrow\uparrow}_{ij}=-\chi^{\uparrow\downarrow\ast}_{ij}$,
$\Delta^{\downarrow\uparrow}_{ij}=-\Delta^{\uparrow\downarrow\ast}_{ij}$, $\Delta^{\downarrow\downarrow}_{ij}=\Delta^{\uparrow\uparrow\ast}_{ij}$; with SU(2),
$\chi$ and $\Delta$ are real.

## Wick's theorem and the mean field

`WickVEV` brings the operator to normal order and applies Wick's theorem to each product, with the sign of the
pairing. `MeanField` replaces each normal-ordered product $O_1\cdots O_{2n}$ by

$$\sum_{p\lt q}(-1)^{p+q+1}\thinspace O_pO_q\thinspace \langle O_1\cdots\widehat{O_p}\cdots\widehat{O_q}\cdots O_{2n}\rangle-(n-1)\thinspace \langle O_1\cdots O_{2n}\rangle,$$

the hats marking the operators left out: the result has the expectation value and the first variation of the operator.

## Matrices

The modes of a `FockSpace` are numbered site by site, then orbital, with the spin fastest (↑ before ↓). The basis
states are the occupation vectors, mode 1 the most significant bit, and the ladder operators are the Jordan–Wigner
matrices $c_m=Z\otimes\cdots\otimes Z\otimes\sigma^-\otimes1\otimes\cdots$, with $Z$ on the modes before $m$.
`LehmannGreen` gives

$$G_{ab}(z)=\langle c_a\thinspace (z-H+E_0)^{-1}\thinspace c^\dagger_b\rangle+\langle c^\dagger_b\thinspace (z+H-E_0)^{-1}\thinspace c_a\rangle,$$

averaged over a degenerate ground multiplet.

## Gutzwiller

`GutzwillerProjector[i]` is $\sqrt2\thinspace [\cos\theta\thinspace P_i^{(1)}+\sin\theta\thinspace (P_i^{(0)}+D_i)]$ at half filling, and
`GutzwillerVEV` is `WickVEV` to first order in the contractions between different sites.
