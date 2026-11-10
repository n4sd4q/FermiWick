dimer = GaussianState[{{cDag[i_, s_], c[j_, t_]} :> KroneckerDelta[s, t] If[i === j, 1/2, g]}];
WickVEV[nOp[1, "↑"] ** nOp[2, "↑"], dimer]