DefineOperator[Pair, 2,
  (c[#1, "↑"] ** c[#2, "↓"] - c[#1, "↓"] ** c[#2, "↑"])/Sqrt[2] &,
  Subscript["b", Row[{#1, #2}]] &];
{Dagger[Hop[1, 2]], Dagger[Pair[1, 2]],
  ToSpinChargeBasis[Dagger[Pair[1, 2]] ** Pair[1, 2]]}