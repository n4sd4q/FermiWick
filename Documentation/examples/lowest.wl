dimer = -Sum[cDag[1, s] ** c[2, s] + cDag[2, s] ** c[1, s], {s, $Spins}] + 4 (DOp[1] + DOp[2]);
N[Round[10^4 First[LowestStates[FockMatrix[dimer, sec], 2]]]/10^4]