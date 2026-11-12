ring4 = -Sum[PlusHC[cDag[i, s] ** c[Mod[i, 4] + 1, s]], {i, 4}, {s, $Spins}] + 2 Sum[(nTotOp[i] - 1)^2, {i, 4}];
g4 = LehmannGreen[ring4, FockSpace[Range[4]], {c[1, "↑"]}, "Temperature" -> 1/2, "ChemicalPotential" -> 0];
{N[Round[10^4 g4["GrandPotential"]]/10^4], Length[g4["Poles"]], N[Round[10^4 Im[GreenMatrix[g4, I Pi/2][[1, 1]]]]/10^4]}