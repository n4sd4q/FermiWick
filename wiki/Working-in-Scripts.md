FermiWick runs the same in a script as in a notebook, which is the way for long calculations, loops over parameters and
batch jobs.

## Running a script

A file `calc.wl`:

<!-- not run -->
```mathematica
Needs["FermiWick`"];
result = WickVEV[SdotS[1, 2] ** SdotS[2, 3], SpinSingletState[]];
Print[InputForm[result]];
Put[result, "result.m"];
```

runs with

<!-- not run -->
```mathematica
wolfram -script calc.wl
```

(`math -script` and `wolframscript -file` do the same where they are installed). To use a clone of the repository
instead of the installed paclet, `PacletDirectoryLoad["/path/to/FermiWick/FermiWick"]` before `Needs`.

In a script the package prints no banner and never looks for a newer release.

## Printing results

A script has no notebook to typeset the output: print it as input, which can be read back, or as LaTeX:

```mathematica
Needs["FermiWick`"];
ToString[InputForm[WickVEV[SdotS[1, 2], SpinSingletState[]]]]
```

```text
(-3*Conjugate[χ[1, 2]]*χ[1, 2])/2
```

`Put[expr, "file.m"]` saves a result and `Get["file.m"]` reads it back, in a later session with the package loaded.

## Speed

Most of the time of a large calculation goes to normal ordering and to Wick's theorem; both work on lists of terms, so
it pays to expand a sum into terms you need and to set numbers early (`/. t -> 1`) for matrices. `FockMatrix` keeps the
matrices of the operators it meets in a cache, which `ClearFockCache[]` empties when memory runs short.