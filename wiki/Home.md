# FermiWick

Practical guides and a FAQ for FermiWick, a Wolfram Language package that does the algebra of fermion operators: you
write operators as on paper, and it brings them to normal order, evaluates them by Wick's theorem in the state you
assume, gives their mean field, and solves small clusters exactly. The
[website](https://n4sd4q.github.io/FermiWick/) gives the overview and a first example, the
[guided tour](https://github.com/n4sd4q/FermiWick/blob/main/Examples.wl) builds everything up step by step, the
[showcase](https://n4sd4q.github.io/FermiWick/Showcase/FermiWick_Showcase.pdf) goes through the physics and the
[manual](https://n4sd4q.github.io/FermiWick/Documentation/FermiWick_Manual.pdf) describes every function.

## Pages

- Start
  - [Installation](Installation): installing, loading, updating, and what to do when loading fails
  - [Writing operators](Writing-Operators): ladder operators, sums, local operators, spins, "+ h.c.", any labels
  - [FAQ](FAQ): short answers to the usual questions
- Recipes
  - [Declaring a state](Declaring-a-State): Gaussian states by their symmetries or by rules, and what they contract
  - [Mean field](Mean-Field): the Hartree–Fock–Bogoliubov Hamiltonian of an operator, its fields, a self-consistent loop
  - [Exact diagonalisation](Exact-Diagonalisation): matrices on Fock spaces and sectors, spectra, correlations, Green's
    functions, spin models
  - [Gutzwiller projectors](Gutzwiller-Projectors): projected operators and the Gutzwiller approximation
  - [LaTeX output](LaTeX-Output): formulas for a paper, and Wick's theorem written out
  - [Working in scripts](Working-in-Scripts): running the package without a notebook
- Reference
  - [Conventions](Conventions): operators, contractions and signs, in one place

The code on these pages is run, and gives what is shown after it: a formula as the package writes it in LaTeX
(`ToTeX`), or a number as a terminal prints it.
