Function MyArrFx (X As Double, Y() As Double, I As Integer) As Double

  Select Case I
    Case 0
      MyArrFx = Y(1)
      
    Case 1
      MyArrFx = Y(0) + X

    Case Else
      MyArrFx = 1
  End Select
End Function

Function MyFx (X As Double, Y As Double) As Double
  MyFx = -X * Y ^ 2
End Function

