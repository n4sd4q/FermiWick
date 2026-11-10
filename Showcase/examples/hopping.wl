real = GaussianState[{"SpinRotation", "Charge", "TimeReversal"}, "Filling" -> 1/2];
hop = Sum[cDag[1, s] ** c[2, s] + cDag[2, s] ** c[1, s], {s, $Spins}];
Simplify[GutzwillerVEV[LocalSandwich[GutzwillerProjector, {1, 2}, hop], real]/WickVEV[hop, real]]