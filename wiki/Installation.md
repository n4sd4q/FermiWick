## Installing

FermiWick needs Mathematica 13.0 or later and nothing else. Install the paclet of the latest
[release](https://github.com/n4sd4q/FermiWick/releases) once:

<!-- not run -->
```mathematica
PacletInstall["https://github.com/n4sd4q/FermiWick/releases/download/v0.3.0/FermiWick-0.3.0.paclet"]
```

Without a connection, download the `.paclet` file from the release page and give `PacletInstall` its path.

## Loading

In every session:

```mathematica
Needs["FermiWick`"]
```

and, when needed, the two modules:

```mathematica
Needs["FermiWick`Fock`"];         (* matrices and exact diagonalisation *)
Needs["FermiWick`Gutzwiller`"];   (* Gutzwiller projectors *)
$FermiWickVersion
```

```text
0.3.0
```

The guided tour, a notebook to run cell by cell, comes with the paclet:

<!-- not run -->
```mathematica
NotebookOpen[PacletObject["FermiWick"]["AssetLocation", "Examples"]]
```

## From a clone of the repository

To use a clone instead of the installed paclet (to try the latest changes, or to change the package):

<!-- not run -->
```mathematica
PacletDirectoryLoad["/path/to/FermiWick/FermiWick"];
Needs["FermiWick`"]
```

## Updating

In a notebook the package looks for a newer release on GitHub at most once a day, in the background, and the banner then
shows the line that installs it. `FermiWickUpdateCheck[]` asks by hand; `$FermiWickUpdateCheck = False` before `Needs`
turns the check off ([`SECURITY.md`](https://github.com/n4sd4q/FermiWick/blob/main/SECURITY.md) says what it sends).
Installing a newer release replaces the older one; `PacletUninstall["FermiWick"]` removes the package.

## When loading fails

- `Needs::nocont` or "cannot open": the paclet is not installed, or not found. `PacletFind["FermiWick"]` lists what
  Mathematica sees; install again, or load the clone with `PacletDirectoryLoad`.
- An older version loads: `PacletFind["FermiWick"]` lists more than one; `PacletUninstall` the older ones.
- Shadowing messages (`::shdw`): a symbol of the package was used before loading it. Restart the kernel, then load the
  package first.