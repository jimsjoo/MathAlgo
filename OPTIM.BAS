Global Const OPTIM_BAD_RESULT# = -1E+31

Function BisectionMin (Xa As Double, Xb As Double, tolerance As Double) As Double

  Dim Xc As Double
  Dim Fa As Double, Fb As Double, Fc As Double

  Fa = MyFx(Xa)
  Fb = MyFx(Xb)
  Do
    Xc = (Xa + Xb) / 2
    Fc = MyFx(Xc)
    If slope(Xc) * slope(Xa) > 0 Then
      Xa = Xc
      Fc = Fa
    Else
      Xb = Xc
      Fc = Fb
    End If
  Loop While Abs(Xb - Xa) > tolerance
  BisectionMin = Xc
End Function

Sub CalcYSimplex (X() As Double, Y() As Double, numVars As Integer)

  Dim numPoints As Integer
  Dim j As Integer, i As Integer
  Dim Xarr() As Double

  numPoints = numVars + 1
  ReDim Xarr(numVars)
  
  ' calculate the values for y() using sumX()
  For i = 0 To numPoints - 1
    For j = 0 To numVars - 1
      Xarr(j) = X(i, j)
    Next j
    Y(i) = MySimplexFx(Xarr())
  Next i
End Sub

Function CubicIntMin (Xa As Double, Xb As Double, Xtol As Double, Gtol As Double, maxIter As Integer) As Double
  Dim Fa As Double, Fb As Double, Fx As Double
  Dim Ga As Double, Gb As Double, Gx As Double
  Dim w As Double, X As Double, v As Double
  Dim h As Double, Gmin As Double
  Dim iter As Integer

  iter = 0
  Fa = MyFx(Xa)
  Fb = MyFx(Xb)
  ' calculate slope at Xa
  Ga = slope(Xa)
  ' calculate slope at Xb
  Gb = slope(Xb)
  Do
    iter = iter + 1
    w = 3 / (Xb - Xa) * (Fa - Fb) + Ga + Gb
    v = Sqr(w ^ 2 - Ga * Gb)
    X = Xa + (Xb - Xa) * (1 - (Gb + v - w) / (Gb - Ga + 2 * v))
    Fx = MyFx(X)
    ' calculate slope at x
    Gx = slope(X)
    If Ga < 0 And (Gx > 0 Or Fx > Fa) Then
      Xb = X
      Fb = Fx
      Gb = Gx
    Else
      Xa = X
      Fa = Fx
      Ga = Gx
    End If
    If Abs(Ga) > Abs(Gb) Then
      Gmin = Abs(Gb)
    Else
      Gmin = Abs(Ga)
    End If
  Loop While Not (Abs(Xb - Xa) < Xtol Or Gmin < Gtol Or iter > maxIter)
  If iter <= maxIter Then
    CubicIntMin = X
  Else
    CubicIntMin = OPTIM_BAD_RESULT
  End If
End Function

Function ExNewtonMin (X() As Double, tolerance As Double, diff As Double, maxIter As Integer, index As Integer) As Integer

  Dim h As Double
  Dim Fm As Double, F0 As Double
  Dim Fp As Double
  Dim firstDeriv As Double, secondDeriv As Double
  Dim iter As Integer
  Dim xOld As Double
  Dim xx As Double

  iter = 0
  xOld = X(index)
  xx = xOld
  Do
    ' calculate increment
    If Abs(xx) > 1 Then h = .01 * xx Else h = .01
    ' calculate function values at x-h, X, and x+h
    X(index) = xx - h
    Fm = MySimplexFx(X())
    X(index) = xx + h
    Fp = MySimplexFx(X())
    X(index) = xx
    F0 = MySimplexFx(X())
    ' calculate the first derivative
    firstDeriv = (Fp - Fm) / 2 / h
    ' calculate the second derivative
    secondDeriv = (Fp - 2 * F0 + Fm) / h ^ 2
    ' calculate the guess refinement
    diff = firstDeriv / secondDeriv
    ' refine the guess
    xx = xx - diff
    X(index) = xx
    iter = iter + 1 ' the increment iteration counter
  Loop While Abs(diff) > tolerance And iter < maxIter

  If iter <= maxIter Then
    ExNewtonMin = True
  Else
    X(index) = xOld
    ExNewtonMin = False
  End If
End Function

Function GoldenSearchMin (Xa As Double, Xb As Double, tolerance As Double, maxIter As Integer) As Double

  Dim Xc As Double, Xd As Double
  Dim Fc As Double, Fd As Double
  Dim oneMinusTau As Double
  Dim iter As Integer

  iter = 0
  oneMinusTau = 1 - (Sqr(5) - 1) / 2
  Xc = Xa + oneMinusTau * (Xb - Xa)
  Fc = MyFx(Xc)
  Xd = Xb - oneMinusTau * (Xb - Xa)
  Fd = MyFx(Xd)
  Do
    iter = iter + 1
    If Fc < Fd Then
       Xb = Xd
       Xd = Xc
       Xc = Xa + oneMinusTau * (Xb - Xa)
       Fd = Fc
       Fc = MyFx(Xc)
    Else
       Xa = Xc
       Xc = Xd
       Xd = Xb - oneMinusTau * (Xb - Xa)
       Fc = Fd
       Fd = MyFx(Xd)
    End If
  Loop While Abs(Xb - Xa) > tolerance And iter < maxIter
  If iter <= maxIter Then
    GoldenSearchMin = Xc
  Else
    GoldenSearchMin = OPTIM_BAD_RESULT
  End If
End Function

Function NewtonMin (X As Double, tolerance As Double, maxIter As Integer) As Double

  Dim h As Double, diff As Double
  Dim Fm As Double, F0 As Double, Fp As Double
  Dim firstDeriv As Double, secondDeriv As Double
  Dim iter As Integer

  iter = 0
  Do
    ' caluclate increment
    If Abs(X) > 1 Then h = .01 * X Else h = .01
    ' calculate function values at X-h, X, and X+h
    Fm = MyFx(X - h)
    F0 = MyFx(X)
    Fp = MyFx(X + h)
    ' calculate the first derivative
    firstDeriv = (Fp - Fm) / 2 / h
    ' calculate the second derivative
    secondDeriv = (Fp - 2 * F0 + Fm) / h ^ 2
    ' calculate the guess refinement
    diff = firstDeriv / secondDeriv
    ' refine the guess
    X = X - diff
    iter = iter + 1 ' the increment iteration counter
  Loop While Abs(diff) > tolerance And iter < maxIter

  If iter <= maxIter Then
    NewtonMin = X
  Else
    NewtonMin = OPTIM_BAD_RESULT
  End If
End Function

Function NewtonMultiMin (X() As Double, numVars As Integer, tolerance As Double, maxIter As Integer) As Integer

  Dim i As Integer, j As Integer
  Dim k As Integer
  Dim diff As Double
  Dim ok As Integer

  Do
    ok = False
    For i = 0 To numVars - 1
      If ExNewtonMin(X(), tolerance, diff, maxIter, i) Then
        If diff > tolerance Then
          ok = True
        End If
      Else
        NewtonMultiMin = False
        Exit Function
      End If
    Next i
  Loop While ok
  NewtonMultiMin = True
End Function

Function QuadIntMin (Xa As Double, Xb As Double, Xc As Double, Xtol As Double, Ftol As Double, maxIter As Integer) As Double

  Dim Fa As Double, Fb As Double, Fc As Double
  Dim lastFc As Double, X As Double, Fx As Double
  Dim iter As Integer
  Dim ok As Integer

  iter = 0
  ok = False
  Fa = MyFx(Xa)
  Fb = MyFx(Xb)
  Fc = MyFx(Xc)

  Do
    iter = iter + 1
    lastFc = Fc
    X = .5 * ((Xb ^ 2 - Xc ^ 2) * Fa + (Xc ^ 2 - Xa ^ 2) * Fb + (Xa ^ 2 - Xb ^ 2) * Fc) / ((Xb - Xc) * Fa + (Xc - Xa) * Fb + (Xa - Xb) * Fc)
    Fx = MyFx(X)
    If X < Xc And Fx < Fc Then
      Xb = Xc
      Xc = X
      Fb = Fc
      Fc = Fx
    ElseIf X > Xc And Fx > Fc Then
      Xb = X
      Fb = Fx
    ElseIf X < Xc And Fx > Fc Then
      Xa = X
      Fa = Fx
    Else
      Xa = Xc
      Xc = X
      Fa = Fc
      Fc = Fx
    End If
    If Abs(Xa - Xc) < Xtol Or Abs(Xc - Xb) < Xtol Then
      ok = True
    End If
  Loop While Not (ok Or Abs((lastFc - Fc) / lastFc) < Ftol And iter > maxIter)
  If iter <= maxIter Then
    QuadIntMin = X
  Else
    QuadIntMin = OPTIM_BAD_RESULT
  End If
End Function

Sub Simplex (X() As Double, Y() As Double, numVars As Integer, tolerance As Double, reflectFact As Double, expandFact As Double, contractFact As Double, maxIter As Integer)

  Const EPS# = .00000001
  Const defReflectFact# = 1#
  Const defExpandFact# = 2#
  Const defContractFact# = .5

  Dim numIter As Integer
  Dim numPoints As Integer
  Dim j As Integer, i As Integer
  Dim worstI As Integer, bestI As Integer
  Dim goOn As Integer
  Dim flag As Integer
  Dim y1 As Double, y2 As Double, x0 As Double
  Dim yMean As Double, sum As Double
  Dim X1() As Double
  Dim X2() As Double
  Dim Centroid() As Double

  numIter = 0
  numPoints = numVars + 1
  goOn = True
  
  ' allocate dynamic arrays
  ReDim X1(numVars)
  ReDim X2(numVars)
  ReDim Centroid(numVars)

  ' check Simplex modification factors
  If reflectFact < EPS Then
    reflectFact = defReflectFact
  End If
  If expandFact < EPS Then
    expandFact = defExpandFact
  End If
  If contractFact < EPS Then
    contractFact = defContractFact
  End If
  ' calculate the values for y() using X1()
  For i = 0 To numPoints - 1
    For j = 0 To numVars - 1
      X1(j) = X(i, j)
    Next j
    Y(i) = MySimplexFx(X1())
  Next i

  Do While numIter < maxIter And goOn
    numIter = numIter + 1
    ' find worst and best point
    worstI = 0
    bestI = 0
    For i = 1 To numPoints - 1
      If Y(i) < Y(bestI) Then
        bestI = i
      ElseIf Y(i) > Y(worstI) Then
        worstI = i
      End If
    Next i
    ' calculate centroid (exclude worst point)
    For i = 0 To numVars - 1
      Centroid(i) = 0
      For j = 0 To numPoints - 1
        If j <> worstI Then
          Centroid(i) = Centroid(i) + X(j, i)
        End If
      Next j
      Centroid(i) = Centroid(i) / numVars
    Next i
    ' calculate reflected point
    For i = 0 To numVars - 1
      X1(i) = (1 + reflectFact) * Centroid(i) - reflectFact * X(worstI, i)
    Next i
    y1 = MySimplexFx(X1())

    If y1 < Y(bestI) Then
      ' calculate expanded point
      For i = 0 To numVars - 1
        X2(i) = (1 + expandFact) * X1(i) - expandFact * Centroid(i)
      Next i
      y2 = MySimplexFx(X2())
      If y2 < Y(bestI) Then
        ' replace worst point by X2
        For i = 0 To numVars - 1
          X(worstI, i) = X2(i)
        Next i
      Else
        ' replace worst point by X1
        For i = 0 To numVars - 1
          X(worstI, i) = X1(i)
        Next i
      End If
    Else
      flag = True
      For i = 0 To numPoints - 1
        If i <> worstI And y1 <= Y(i) Then
          flag = False
          Exit For
        End If
      Next i
      If flag Then
        If y1 < Y(worstI) Then
          ' replace worst point by X1
          For i = 0 To numVars - 1
            X(worstI, i) = X1(i)
          Next i
          Y(worstI) = y1
        End If
        ' calculate contracted point
        For i = 0 To numVars - 1
          X2(i) = contractFact * X(worstI, i) + (1 - contractFact) * Centroid(i)
        Next i
        y2 = MySimplexFx(X2())
        If y2 > Y(worstI) Then
          ' store best x in X1
          For i = 0 To numVars - 1
            X2(i) = X(bestI, i)
          Next i
          For j = 0 To numPoints - 1
            For i = 0 To numVars - 1
              X(j, i) = .5 * (X2(i) + X(j, i))
            Next i
          Next j
        Else
          ' replace worst point by X2
          For i = 0 To numVars - 1
            X(worstI, i) = X2(i)
          Next i
        End If
      Else
        ' replace worst point by X1
        For i = 0 To numVars - 1
          X(worstI, i) = X1(i)
        Next i
      End If
    End If
    ' calculate the values for y() using X1
    For i = 0 To numPoints - 1
      For j = 0 To numVars - 1
        X1(j) = X(i, j)
      Next j
      Y(i) = MySimplexFx(X1())
    Next i
    ' calculate mean y
    yMean = 0
    For i = 0 To numPoints - 1
      yMean = yMean + Y(i)
    Next i
    yMean = yMean / numPoints
    ' calculate deviation fron mean y
    sum = 0
    For i = 0 To numPoints - 1
      sum = sum + (Y(i) - yMean) ^ 2
    Next i
    ' test convergence
    If Sqr(sum / numPoints) > tolerance Then
      goOn = True
    Else
      goOn = False
    End If
  Loop

  ' find the best point
  bestI = 0
  For i = 1 To numPoints - 1
    If Y(i) < Y(bestI) Then bestI = i
  Next i
  If bestI <> 0 Then
    ' move best point to index 0
    For i = 0 To numVars - 1
       x0 = X(0, i)
       X(0, i) = X(bestI, i)
       X(bestI, i) = x0
    Next i
    y1 = Y(0)
    Y(0) = Y(bestI)
    Y(bestI) = y1
  End If
End Sub

Private Function slope (X As Double) As Double
  Dim h As Double
  If Abs(X) > 1 Then
    h = .01 * X
  Else
    h = .01
  End If
  slope = (MyFx(X + h) - MyFx(X - h)) / 2 / h
End Function

