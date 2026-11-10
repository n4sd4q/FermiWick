DefineOperator[Hop, 2, (cDag[#1, "↑"] ** c[#2, "↑"] + cDag[#1, "↓"] ** c[#2, "↓"]) &,
  Subscript["T", Row[{#1, #2}]] &];
ToSpinChargeBasis[OperatorCommutator[Hop[1, 2], nTotOp[2]]]