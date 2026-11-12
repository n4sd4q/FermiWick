space = FockSpace[{1, 2}];
g = LehmannGreen[dimer, space, {c[1, "↑"]}, "ChemicalPotential" -> 0];
{N[Round[10^4 g["GroundEnergy"]]/10^4], g["Degeneracy"], Max[Abs[Total[g["Residues"]] - IdentityMatrix[1]]] < 10^-12}