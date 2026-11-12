g = LehmannGreen[H, space, {c[1, "↑"]}];
FullSimplify[GreenMatrix[g, z][[1, 1]], t > 0 && U > 0]