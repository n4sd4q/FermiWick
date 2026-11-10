flip = GaussianState[{{cDag[i_, s_], c[j_, t_]} :> If[i === j, KroneckerDelta[s, t]/2, If[s === t, g, h]]}];
WickVEV[SzOp[1] ** SzOp[2], flip]