gT = LehmannGreen[dimer, space, {c[1, "↑"]}, "Temperature" -> 1/2];
{N[Round[10^4 gT["GrandPotential"]]/10^4], Length[gT["Poles"]], Max[Abs[Total[gT["Residues"]] - IdentityMatrix[1]]] < 10^-12}