six = FockSector[FockSpace[Range[8]], "N" -> 6, "Sz" -> 0];
v0 = First[Last[LowestStates[FockMatrix[hop8, six], 1]]];
pairs = Flatten[Table[{cDag[i, s], c[j, s]}, {i, 8}, {j, 8}, {s, $Spins}], 2];
exact = GaussianState[Thread[pairs -> (FockExpectation[#1 ** #2, v0, six] & @@@ pairs)]];
op = SdotS[1, 2] + DOp[1];
Chop[WickVEV[op, exact] - FockExpectation[op, v0, six]]