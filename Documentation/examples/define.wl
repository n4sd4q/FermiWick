DefineOperator[Hop, 2, (cDag[#1, "↑"] ** c[#2, "↑"] + cDag[#1, "↓"] ** c[#2, "↓"]) &,
  Subscript["T", Row[{#1, #2}]] &];
ToSpinChargeBasis[NormalOrder[2 Hop[1, 2] ** nTotOp[3]]]