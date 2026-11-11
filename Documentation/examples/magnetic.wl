magnetic = GaussianState[{"SpinRotationZ", "Charge"}, "Filling" -> 1/2, "Magnetization" -> m];
MatrixForm[ContractionTable[magnetic, {1}, "Normal"]]