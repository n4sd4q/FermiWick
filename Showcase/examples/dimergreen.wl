g = LehmannGreen[H, space, {c[1, "↑"], c[2, "↑"]}, "ChemicalPotential" -> μ,
  "Assumptions" -> t > 0 && U > 0 && Abs[μ] < Sqrt[16 t^2 + U^2]/2 - t];
First[g["GroundSectors"]]