gT = LehmannGreen[H, space, {c[1, "↑"], c[2, "↑"]}, "Temperature" -> T, "ChemicalPotential" -> μ];
bonding = Select[Transpose[{gT["Poles"], Z (P . # . P)[[1, 1]] & /@ gT["Residues"]}], ! PossibleZeroQ[Last[#]] &];
MatrixForm[FullSimplify[bonding, T > 0 && t > 0 && U > 0]]