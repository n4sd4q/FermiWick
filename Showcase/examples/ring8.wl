nb[i_] := Mod[i - 1, 8] + 1;
hop8 = -Sum[cDag[i, s] ** c[nb[i + 1], s] + cDag[nb[i + 1], s] ** c[i, s], {i, 8}, {s, $Spins}];
half = FockSector[FockSpace[Range[8]], "N" -> 8, "Sz" -> 0];
{e, v} = LowestStates[FockMatrix[hop8 + 4 Sum[DOp[i], {i, 8}], half], 1];
N[Round[10^4 {First[e]/8, FockExpectation[DOp[1], First[v], half],
    FockExpectation[SdotS[1, 2], First[v], half]}]/10^4]