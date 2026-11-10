DeclareModeLabels[{"Momentum"}];
WickVEV[cDag[k] ** c[k] ** cDag[q] ** c[q], GaussianState[{{cDag[p1_], c[p2_]} :> KroneckerDelta[p1, p2] nk[p1]}]]