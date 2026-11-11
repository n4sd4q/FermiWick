{vals, vecs} = Eigensystem[Normal[FockMatrix[H, sector]]];
ground = Normalize[vecs[[First[Ordering[N[vals /. {t -> 1, U -> 1}], 1]]]]];
FullSimplify[StateToKet[sector, ground], t > 0 && U > 0]