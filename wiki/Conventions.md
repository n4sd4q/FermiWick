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
- `Dagger` takes the amplitudes χ, Δ, ρ, κ as complex and every other symbol as real, unless listed with `"Complex"`.

## Contractions

A Gaussian state is written in spin channels, with $\sigma^a$ the Pauli matrices and $\epsilon_{\uparrow\downarrow}=1$:

$$\langle c^\dagger_{i\sigma}c_{j\sigma'}\rangle=\chi_{ij}\delta_{\sigma\sigma'}+\sum_a\chi^a_{ij}\sigma^a_{\sigma'\sigma},\qquad
\langle c_{i\sigma}c_{j\sigma'}\rangle=\sum_\tau\Big(\Delta_{ij}\delta_{\sigma\tau}+\sum_a\Delta^a_{ij}\sigma^a_{\sigma\tau}\Big)\epsilon_{\tau\sigma'}.$$

Hermiticity gives $\chi_{ji}=\chi_{ij}^{\ast}$ and $\chi^a_{ji}={\chi^a_{ij}}^{\ast}$; antisymmetry $\Delta_{ji}=\Delta_{ij}$ and
$\Delta^a_{ji}=-\Delta^a_{ij}$. On a site, $\chi_{ii}$ is the density per spin (`"Filling"`), $\chi^a_{ii}=\langle S^a_i\rangle$
(`"Magnetization"`) and $\Delta_{ii}=\langle c_{i\uparrow}c_{i\downarrow}\rangle$ (`"OnSitePairing"`). In a result,
$\chi_{ji}$ is written $\chi_{ij}^{\ast}$.

## Wick's theorem and the mean field

`WickVEV` brings the operator to normal order and applies Wick's theorem to each product, with the sign of the
pairing. `MeanField` replaces each normal-ordered product $O_1\cdots O_{2n}$ by

$$\sum_{p<q}(-1)^{p+q+1}\thinspace O_pO_q\thinspace \langle O_1\cdots\widehat{O_p}\cdots\widehat{O_q}\cdots O_{2n}\rangle-(n-1)\thinspace \langle O_1\cdots O_{2n}\rangle,$$

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
