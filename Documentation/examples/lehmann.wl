space = FockSpace[{1, 2}];
g = LehmannGreen[FockMatrix[dimer, space], space, {c[1, "↑"]}, "N" -> 2];
{N[Round[10^4 g["GroundEnergy"]]/10^4], g["Degeneracy"], Max[Abs[Total[g["Residues"]] - IdentityMatrix[1]]] < 10^-12}