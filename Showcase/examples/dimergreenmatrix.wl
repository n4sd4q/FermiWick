P = {{1, 1}, {1, -1}}/Sqrt[2];
MatrixForm[Simplify[P . GreenMatrix[g, z] . P, t > 0 && U > 0]]