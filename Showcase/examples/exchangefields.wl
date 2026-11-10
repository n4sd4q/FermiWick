mf = MeanField[J SdotS[1, 2], state];
fields = OperatorCoefficients[mf];
{fields[cDag[1, "↑"] ** c[2, "↑"]], fields[cDag[1, "↑"] ** cDag[2, "↓"]], fields[1]}