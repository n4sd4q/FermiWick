atom = U/2 (nTotOp[1] - 1)^2;
FullSimplify[GreenMatrix[LehmannGreen[atom, FockSpace[{1}], {c[1, "↑"]}, "Temperature" -> T], z][[1, 1]]]