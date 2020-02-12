Function My2Fx (X As Double, Y As Double, params() As Double) As Double
  My2Fx = X * X + Y * Y + params(0)
End Function

Function My2Fy (X As Double, Y As Double, params() As Double) As Double
  My2Fy = X * X - Y * Y + params(0)
End Function

Function MyDeriv1 (X As Double, params() As Double) As Double
  MyDeriv1 = Exp(X) - 2 * params(0) * X
End Function

Function MyDeriv2 (X As Double, params() As Double) As Double
  MyDeriv2 = Exp(X) - 6
End Function

Function MyFx (X As Double, params() As Double) As Double
  MyFx = Exp(X) - params(0) * X ^ 2
End Function

Function MySNE (X() As Double, index As Integer) As Double
  Select Case index
    Case 0
      MySNE = X(1) * Exp(X(0)) - 2

    Case 1
      MySNE = X(0) ^ 2 + X(1) - 4

    Case Else
      MySNE = 1
  End Select
End Function


