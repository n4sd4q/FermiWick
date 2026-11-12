weights = AssociationThread[gT["Poles"], gT["Residues"][[All, 1, 1]]];
FullSimplify[weights[t + U/2], T > 0]