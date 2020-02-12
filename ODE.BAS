
Sub RungeKutta4 (x0 As Double, y0 As Double, deltaX As Double, yArr() As Double, numElem As Integer)

  Dim x As Double, y As Double, halfDelta As Double
  Dim k1 As Double, k2 As Double
  Dim k3 As Double, k4 As Double
  Dim i As Integer

  If numElem < 1 Then Exit Sub

  x = x0
  y = y0
  halfDelta = deltaX / 2
  For i = 0 To numElem - 1
     k1 = deltaX * MyFx(x, y)
     k2 = deltaX * MyFx(x + halfDelta, y + k1 / 2)
     k3 = deltaX * MyFx(x + halfDelta, y + k2 / 2)
     k4 = deltaX * MyFx(x + deltaX, y + k3)
     y = y + (k1 + 2 * (k2 + k3) + k4) / 6
     yArr(i) = y
     x = x + deltaX
  Next i
End Sub

Sub RungeKuttaFehlberg (x0 As Double, y0 As Double, deltaX As Double, yArr() As Double, numElem As Integer)

  Dim x As Double, y As Double
  Dim k1 As Double, k2 As Double
  Dim k3 As Double, k4 As Double, k5 As Double
  Dim i As Integer

  If numElem < 1 Then Exit Sub

  x = x0
  y = y0
  For i = 0 To numElem - 1
     k1 = deltaX * MyFx(x, y)
     k2 = deltaX * MyFx(x + deltaX / 4, y + k1 / 4)
     k3 = deltaX * MyFx(x + 3# / 8 * deltaX, y + 3# / 32 * (k1 + 3 * k2))
     k4 = deltaX * MyFx(x + 12# / 13 * deltaX, y + 1932# / 2197 * k1 - 7200# / 2197 * k2 + 7296# / 2197 * k3)
     k5 = deltaX * MyFx(x + deltaX, y + 439# / 216 * k1 - 8 * k2 + 3680# / 513 * k3 - 845# / 4104 * k4)
     y = y + 25# / 216 * k1 + 1408# / 2565 * k3 + 2197# / 4104 * k4 - k5 / 5
     yArr(i) = y
     x = x + deltaX
  Next i
End Sub

Sub RungeKuttaGill (x0 As Double, y0 As Double, deltaX As Double, yArr() As Double, numElem As Integer)

  Dim x As Double, y As Double, halfDelta As Double
  Dim k1 As Double, k2 As Double
  Dim k3 As Double, k4 As Double
  Dim c1 As Double, c2 As Double
  Dim c3 As Double, c4 As Double
  Dim c5 As Double
  Dim i As Integer

  c1 = 1 / Sqr(2)
  c2 = 1 + c1
  c3 = 1 - c1
  c4 = -.5 + c1
  c5 = -c1

  If numElem < 1 Then Exit Sub

  x = x0
  y = y0
  halfDelta = deltaX / 2
  For i = 0 To numElem - 1
     k1 = deltaX * MyFx(x, y)
     k2 = deltaX * MyFx(x + halfDelta, y + k1 / 2)
     k3 = deltaX * MyFx(x + halfDelta, y + c4 * k1 + c3 * k2)
     k4 = deltaX * MyFx(x + deltaX, y + c5 * k2 + c2 * k3)
     y = y + (k1 + k4) / 6 + (c3 * k2 + c2 * k3) / 3
     yArr(i) = y
     x = x + deltaX
  Next i
End Sub

Sub VectRungeKutta4 (x0 As Double, y0() As Double, deltaX As Double, y1() As Double, numYVar As Integer)

  Dim x As Double, halfDelta As Double
  Dim y() As Double
  Dim k1() As Double, k2() As Double
  Dim k3() As Double, k4() As Double
  Dim i As Integer

  If numYVar < 1 Then Exit Sub

  ' allocate local dynamic arrays
  ReDim k1(numYVar)
  ReDim k2(numYVar)
  ReDim k3(numYVar)
  ReDim k4(numYVar)
  ReDim y(numYVar)

  x = x0
  For i = 0 To numYVar - 1
     y(i) = y0(i)
  Next i
  halfDelta = deltaX / 2
  For i = 0 To numYVar - 1
     k1(i) = deltaX * MyArrFx(x, y(), i)
  Next i
  For i = 0 To numYVar - 1
     y(i) = y(i) + k1(i) / 2
  Next i
  For i = 0 To numYVar - 1
     k2(i) = deltaX * MyArrFx(x + halfDelta, y(), i)
  Next i
  For i = 0 To numYVar - 1
     y(i) = y0(i) + k2(i) / 2
  Next i
  For i = 0 To numYVar - 1
     k3(i) = deltaX * MyArrFx(x + halfDelta, y(), i)
  Next i
  For i = 0 To numYVar - 1
     y(i) = y0(i) + k3(i)
  Next i
  For i = 0 To numYVar - 1
     k4(i) = deltaX * MyArrFx(x + deltaX, y(), i)
  Next i
  For i = 0 To numYVar - 1
     y1(i) = y0(i) + (k1(i) + 2 * (k2(i) + k3(i)) + k4(i)) / 6
  Next i
End Sub

Sub VectRungeKuttaFehlberg (x0 As Double, y0() As Double, deltaX As Double, y1() As Double, numYVar As Integer)

  Dim x As Double, halfDelta As Double
  Dim y() As Double
  Dim k1() As Double, k2() As Double
  Dim k3() As Double, k4() As Double
  Dim k5() As Double
  Dim i As Integer

  If numYVar < 1 Then Exit Sub

  ' allocate local dynamic arrays
  ReDim k1(numYVar)
  ReDim k2(numYVar)
  ReDim k3(numYVar)
  ReDim k4(numYVar)
  ReDim k5(numYVar)
  ReDim y(numYVar)

  x = x0
  For i = 0 To numYVar - 1
     y(i) = y0(i)
  Next i
  For i = 0 To numYVar - 1
     k1(i) = deltaX * MyArrFx(x, y(), i)
  Next i
  For i = 0 To numYVar - 1
     y(i) = y0(i) + k1(i) / 4
  Next i
  For i = 0 To numYVar - 1
     k2(i) = deltaX * MyArrFx(x + deltaX / 4, y(), i)
  Next i
  For i = 0 To numYVar - 1
     y(i) = y0(i) + 3# / 32 * (k1(i) + 3 * k2(i))
  Next i
  For i = 0 To numYVar - 1
     k3(i) = deltaX * MyArrFx(x + 3# / 8 * deltaX, y(), i)
  Next i
  For i = 0 To numYVar - 1
     y(i) = y0(i) + 1932# / 2197 * k1(i) - 7200# / 2197 * k2(i) + 7296# / 2197 * k3(i)
  Next i
  For i = 0 To numYVar - 1
     k4(i) = deltaX * MyArrFx(x + 12# / 13 * deltaX, y(), i)
  Next i
  For i = 0 To numYVar - 1
     y(i) = y0(i) + 439# / 216 * k1(i) - 8 * k2(i) + 3680# / 513 * k3(i) - 845# / 4104 * k4(i)
  Next i
  For i = 0 To numYVar - 1
     k5(i) = deltaX * MyArrFx(x + deltaX, y(), i)
  Next i
  For i = 0 To numYVar - 1
     y1(i) = y0(i) + 25# / 216 * k1(i) + 1408# / 2565 * k3(i) + 2197# / 4104 * k4(i) - k5(i) / 5
  Next i
End Sub

Sub VectRungeKuttaGill (x0 As Double, y0() As Double, deltaX As Double, y1() As Double, numYVar As Integer)

  Dim x As Double, halfDelta As Double
  Dim y() As Double
  Dim k1() As Double, k2() As Double
  Dim k3() As Double, k4() As Double
  Dim c1 As Double, c2 As Double
  Dim c3 As Double, c4 As Double
  Dim c5 As Double
  Dim i As Integer

  c1 = 1 / Sqr(2)
  c2 = 1 + c1
  c3 = 1 - c1
  c4 = -.5 + c1
  c5 = -c1

  If numYVar < 1 Then Exit Sub

  ' allocate local dynamic arrays
  ReDim k1(numYVar)
  ReDim k2(numYVar)
  ReDim k3(numYVar)
  ReDim k4(numYVar)
  ReDim y(numYVar)

  x = x0
  For i = 0 To numYVar - 1
     y(i) = y0(i)
  Next i
  halfDelta = deltaX / 2
  For i = 0 To numYVar - 1
     k1(i) = deltaX * MyArrFx(x, y(), i)
  Next i
  For i = 0 To numYVar - 1
     y(i) = y(i) + k1(i) / 2
  Next i
  For i = 0 To numYVar - 1
     k2(i) = deltaX * MyArrFx(x + halfDelta, y(), i)
  Next i
  For i = 0 To numYVar - 1
     y(i) = y0(i) + c4 * k1(i) + c3 * k2(i)
  Next i
  For i = 0 To numYVar - 1
     k3(i) = deltaX * MyArrFx(x + halfDelta, y(), i)
  Next i
  For i = 0 To numYVar - 1
     y(i) = y0(i) + c5 * k2(i) + c2 * k3(i)
  Next i
  For i = 0 To numYVar - 1
     k4(i) = deltaX * MyArrFx(x + deltaX, y(), i)
  Next i
  For i = 0 To numYVar - 1
     y1(i) = y0(i) + (k1(i) + k4(i)) / 6 + (c3 * k2(i) + c2 * k3(i)) / 3
  Next i
End Sub

