dimer = FockMatrix[-hop + 4 (DOp[1] + DOp[2]), space];
g = LehmannGreen[dimer, space, {c[1, "↑"], c[2, "↑"]}, "N" -> 2];
N[Round[10^4 Sort[g["Poles"]]]/10^4]