h = H - U/2 Sum[nTotOp[i], {i, 2}];
levels = Eigenvalues[Normal[FockMatrix[h, space]]];
F = FullSimplify[-T Log[Total[Exp[-levels/T]]], T > 0]