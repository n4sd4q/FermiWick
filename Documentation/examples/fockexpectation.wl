sec = FockSector[FockSpace[{1, 2}], "N" -> 2, "Sz" -> 0];
singlet = {0, 1, -1, 0}/Sqrt[2];
{FockExpectation[SdotS[1, 2], singlet, sec], FockExpectation[DOp[1], singlet, sec]}