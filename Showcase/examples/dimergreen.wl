e0 = ground . Normal[FockMatrix[h, sector]] . ground;
part[op_, to_, sign_] := Module[{e, v, w},
  {e, v} = Eigensystem[Normal[FockMatrix[h, to]]];
  w = Normal[FockMatrix[op, {sector, to}]] . ground;
  Sum[(Normalize[v[[n]]] . w)^2/(z - sign (e[[n]] - e0)), {n, Length[e]}]];
FullSimplify[part[cDag[1, "↑"], FockSector[space, "N" -> 3, "Sz" -> 1/2], 1] +
  part[c[1, "↑"], FockSector[space, "N" -> 1, "Sz" -> -1/2], -1], t > 0 && U > 0]