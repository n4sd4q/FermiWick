dimer = -Sum[PlusHC[cDag[1, s] ** c[2, s]], {s, $Spins}] + 2 Sum[(nTotOp[i] - 1)^2, {i, 2}];
N[Round[10^4 First[LowestStates[FockMatrix[dimer, sec], 2]]]/10^4]