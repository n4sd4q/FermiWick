Ω = -T Log[Z];
FullSimplify[{-D[Ω, μ], -D[Ω, {μ, 2}]} /. μ -> 0, T > 0]