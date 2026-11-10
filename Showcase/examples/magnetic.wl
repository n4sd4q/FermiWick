sgn["↑"] = 1; sgn["↓"] = -1;
neel = GaussianState[{{cDag[i_, s_], c[i_, t_]} :> KroneckerDelta[s, t] (1/2 + sgn[s] m[i])}]