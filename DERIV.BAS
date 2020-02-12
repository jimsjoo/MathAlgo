Global Const DERIV_DEFAULT_INCR# = .01
Global Const DERIV_BAD_RESULT# = -1E+30

Function CDiffDeriv1 (yArr() As Double, incr As Double, n As Integer, y0Idx As Integer) As Double

  ' verify arguments for parameters n and y0Idx
  If n < 3 Or y0Idx < 1 Or y0Idx > (n - 2) Or incr = 0 Then
    CDiffDeriv1 = DERIV_BAD_RESULT
    Exit Function
  End If

  CDiffDeriv1 = (yArr(y0Idx + 1) - yArr(y0Idx - 1)) / 2 / incr
End Function

Function CDiffDeriv2 (yArr() As Double, incr As Double, n As Integer, y0Idx As Integer) As Double

  ' verify arguments for parameters n and y0Idx
  If n < 3 Or y0Idx < 1 Or y0Idx > (n - 2) Or incr = 0 Then
    CDiffDeriv2 = DERIV_BAD_RESULT
    Exit Function
  End If

  CDiffDeriv2 = (yArr(y0Idx - 1) - 2 * yArr(y0Idx) + yArr(y0Idx + 1)) / incr ^ 2
End Function

Function CDiffDeriv3 (yArr() As Double, incr As Double, n As Integer, y0Idx As Integer) As Double

  ' verify arguments for parameters n and y0Idx
  If n < 5 Or y0Idx < 2 Or y0Idx > (n - 3) Or incr = 0 Then
    CDiffDeriv3 = DERIV_BAD_RESULT
    Exit Function
  End If

  CDiffDeriv3 = (-1 * yArr(y0Idx - 2) + 2 * yArr(y0Idx - 1) - 2 * yArr(y0Idx + 1) + yArr(y0Idx + 2)) / 2 / incr ^ 3
End Function

Function CDiffDeriv4 (yArr() As Double, incr As Double, n As Integer, y0Idx As Integer) As Double

  ' verify arguments for parameters n and y0Idx
  If n < 5 Or y0Idx < 2 Or y0Idx > (n - 3) Or incr = 0 Then
    CDiffDeriv4 = DERIV_BAD_RESULT
    Exit Function
  End If

  CDiffDeriv4 = (yArr(y0Idx - 2) - 4 * yArr(y0Idx - 1) + 6 * yArr(y0Idx) - 4 * yArr(y0Idx + 1) + yArr(y0Idx + 2)) / incr ^ 4
End Function

Function CUdfDeriv1 (x As Double, incr As Double) As Double

  If incr = 0 Then incr = DERIV_DEFAULT_INCR
  CUdfDeriv1 = (MyFx(x + incr) - MyFx(x - incr)) / 2 / incr
End Function

Function CUdfDeriv2 (x As Double, incr As Double) As Double

  If incr = 0 Then incr = DERIV_DEFAULT_INCR
  CUdfDeriv2 = (MyFx(x - incr) - 2 * MyFx(x) + MyFx(x + incr)) / incr ^ 2
End Function

Function CUdfDeriv3 (x As Double, incr As Double) As Double

  If incr = 0 Then incr = DERIV_DEFAULT_INCR
  CUdfDeriv3 = (-1 * MyFx(x - 2 * incr) + 2 * MyFx(x - incr) - 2 * MyFx(x + incr) + MyFx(x + 2 * incr)) / 2 / incr ^ 3
End Function

Function CUdfDeriv4 (x As Double, incr As Double) As Double

  If incr = 0 Then incr = DERIV_DEFAULT_INCR
  CUdfDeriv4 = (MyFx(x - 2 * incr) - 4 * MyFx(x - incr) + 6 * MyFx(x) - 4 * MyFx(x + incr) + MyFx(x + 2 * incr)) / incr ^ 4
End Function

Function FBDiffDeriv1 (yArr() As Double, incr As Double, n As Integer, y0Idx As Integer) As Double

  ' verify arguments for parameters n and y0Idx
  If n < 3 Or y0Idx < 0 Or y0Idx > (n - 3) Or incr = 0 Then
    FBDiffDeriv1 = DERIV_BAD_RESULT
    Exit Function
  End If

  FBDiffDeriv1 = (-3 * yArr(y0Idx) + 4 * yArr(y0Idx + 1) - yArr(y0Idx + 2)) / 2 / incr
End Function

Function FBDiffDeriv2 (yArr() As Double, incr As Double, n As Integer, y0Idx As Integer) As Double

  ' verify arguments for parameters n and y0Idx
  If n < 4 Or y0Idx < 0 Or y0Idx > (n - 4) Or incr = 0 Then
    FBDiffDeriv2 = DERIV_BAD_RESULT
    Exit Function
  End If

  FBDiffDeriv2 = (2 * yArr(y0Idx) - 5 * yArr(y0Idx + 1) + 4 * yArr(y0Idx + 2) - yArr(y0Idx + 3)) / incr ^ 2
End Function

Function FBDiffDeriv3 (yArr() As Double, incr As Double, n As Integer, y0Idx As Integer) As Double

  ' verify arguments for parameters n and y0Idx
  If n < 5 Or y0Idx < 0 Or y0Idx > (n - 5) Or incr = 0 Then
    FBDiffDeriv3 = DERIV_BAD_RESULT
    Exit Function
  End If

  FBDiffDeriv3 = (-5 * yArr(y0Idx) + 18 * yArr(y0Idx + 1) - 24 * yArr(y0Idx + 2) + 14 * yArr(y0Idx + 3) - 3 * yArr(y0Idx + 4)) / 2 / incr ^ 3
End Function

Function FBDiffDeriv4 (yArr() As Double, incr As Double, n As Integer, y0Idx As Integer) As Double

  ' verify arguments for parameters n and y0Idx
  If n < 6 Or y0Idx < 0 Or y0Idx > (n - 6) Or incr = 0 Then
    FBDiffDeriv4 = DERIV_BAD_RESULT
    Exit Function
  End If

  FBDiffDeriv4 = (3 * yArr(y0Idx) - 14 * yArr(y0Idx + 1) + 26 * yArr(y0Idx + 2) - 24 * yArr(y0Idx + 3) + 11 * yArr(y0Idx + 4) - 2 * yArr(y0Idx + 5)) / incr ^ 4
End Function

Function FBUdfDeriv1 (x As Double, incr As Double) As Double

  If incr = 0 Then incr = DERIV_DEFAULT_INCR
  FBUdfDeriv1 = (-3 * MyFx(x) + 4 * MyFx(x + incr) - MyFx(x + 2 * incr)) / 2 / incr
End Function

Function FBUdfDeriv2 (x As Double, incr As Double) As Double

  If incr = 0 Then incr = DERIV_DEFAULT_INCR
  FBUdfDeriv2 = (2 * MyFx(x) - 5 * MyFx(x + incr) + 4 * MyFx(x + 2 * incr) - MyFx(x + 3 * incr)) / incr ^ 2
End Function

Function FBUdfDeriv3 (x As Double, incr As Double) As Double

  If incr = 0 Then incr = DERIV_DEFAULT_INCR
  FBUdfDeriv3 = (-5 * MyFx(x) + 18 * MyFx(x + incr) - 24 * MyFx(x + 2 * incr) + 14 * MyFx(x + 3 * incr) - 3 * MyFx(x + 4 * incr)) / 2 / incr ^ 3
End Function

Function FBUdfDeriv4 (x As Double, incr As Double) As Double

  If incr = 0 Then incr = DERIV_DEFAULT_INCR
  FBUdfDeriv4 = (3 * MyFx(x) - 14 * MyFx(x + incr) + 26 * MyFx(x + 2 * incr) - 24 * MyFx(x + 3 * incr) + 11 * MyFx(x + 4 * incr) - 2 * MyFx(x + 5 * incr)) / incr ^ 4
End Function

Function XCDiffDeriv1 (yArr() As Double, incr As Double, n As Integer, y0Idx As Integer) As Double

  ' verify arguments for parameters n and y0Idx
  If n < 5 Or y0Idx < 2 Or y0Idx > (n - 3) Or incr = 0 Then
    XCDiffDeriv1 = DERIV_BAD_RESULT
    Exit Function
  End If

  XCDiffDeriv1 = (yArr(y0Idx - 2) - 8 * yArr(y0Idx - 1) + 8 * yArr(y0Idx + 1) - yArr(y0Idx + 2)) / 12 / incr
End Function

Function XCDiffDeriv2 (yArr() As Double, incr As Double, n As Integer, y0Idx As Integer) As Double

  ' verify arguments for parameters n and y0Idx
  If n < 5 Or y0Idx < 2 Or y0Idx > (n - 3) Or incr = 0 Then
    XCDiffDeriv2 = DERIV_BAD_RESULT
    Exit Function
  End If

  XCDiffDeriv2 = (-1 * yArr(y0Idx - 2) + 16 * yArr(y0Idx - 1) - 30 * yArr(y0Idx) + 16 * yArr(y0Idx + 1) - yArr(y0Idx + 2)) / 12 / incr ^ 2
End Function

Function XCDiffDeriv3 (yArr() As Double, incr As Double, n As Integer, y0Idx As Integer) As Double

  ' verify arguments for parameters n and y0Idx
  If n < 7 Or y0Idx < 3 Or y0Idx > (n - 4) Or incr = 0 Then
    XCDiffDeriv3 = DERIV_BAD_RESULT
    Exit Function
  End If

  XCDiffDeriv3 = (yArr(y0Idx - 3) - 8 * yArr(y0Idx - 2) + 13 * yArr(y0Idx - 1) - 13 * yArr(y0Idx + 1) + 8 * yArr(y0Idx + 2) - yArr(y0Idx + 3)) / 8 / incr ^ 3
End Function

Function XCDiffDeriv4 (yArr() As Double, incr As Double, n As Integer, y0Idx As Integer) As Double

  ' verify arguments for parameters n and y0Idx
  If n < 7 Or y0Idx < 3 Or y0Idx > (n - 4) Or incr = 0 Then
    XCDiffDeriv4 = DERIV_BAD_RESULT
    Exit Function
  End If

  XCDiffDeriv4 = (-1 * yArr(y0Idx - 3) + 12 * yArr(y0Idx - 2) - 39 * yArr(y0Idx - 1) + 56 * yArr(y0Idx) - 39 * yArr(y0Idx + 1) + 12 * yArr(y0Idx + 2) - yArr(y0Idx + 3)) / 6 / incr ^ 4
End Function

Function XCUdfDeriv1 (x As Double, incr As Double) As Double

  If incr = 0 Then incr = DERIV_DEFAULT_INCR
  XCUdfDeriv1 = (MyFx(x - 2 * incr) - 8 * MyFx(x - incr) + 8 * MyFx(x + incr) - MyFx(x + 2 * incr)) / 12 / incr
End Function

Function XCUdfDeriv2 (x As Double, incr As Double) As Double

  If incr = 0 Then incr = DERIV_DEFAULT_INCR
  XCUdfDeriv2 = (-MyFx(x - 2 * incr) + 16 * MyFx(x - incr) - 30 * MyFx(x) + 16 * MyFx(x + incr) - MyFx(x + 2 * incr)) / 12 / incr ^ 2
End Function

Function XCUdfDeriv3 (x As Double, incr As Double) As Double

  If incr = 0 Then incr = DERIV_DEFAULT_INCR
  XCUdfDeriv3 = (MyFx(x - 3 * incr) - 8 * MyFx(x - 2 * incr) + 13 * MyFx(x - incr) - 13 * MyFx(x + incr) + 8 * MyFx(x + 2 * incr) - MyFx(x + 3 * incr)) / 8 / incr ^ 3
End Function

Function XCUdfDeriv4 (x As Double, incr As Double) As Double

  If incr = 0 Then incr = DERIV_DEFAULT_INCR
  XCUdfDeriv4 = (-MyFx(x - 3 * incr) + 12 * MyFx(x - 2 * incr) - 39 * MyFx(x - incr) + 56 * MyFx(x) - 39 * MyFx(x + incr) + 12 * MyFx(x + 2 * incr) - MyFx(x + 3 * incr)) / 6 / incr ^ 4
End Function

