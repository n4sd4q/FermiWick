H = -t Sum[PlusHC[cDag[1, s] ** c[2, s]], {s, $Spins}] +
  U Sum[nOp[i, SpinUp] ** nOp[i, SpinDown], {i, 2}];
Eigenvalues[Normal[FockMatrix[H, sector]]]