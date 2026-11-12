H = -t Sum[PlusHC[cDag[1, s] ** c[2, s]], {s, $Spins}] + U/2 Sum[(nTotOp[i] - 1)^2, {i, 2}];
Eigenvalues[Normal[FockMatrix[H, sector]]]