Global Const POLYREG_BIG# = 1E+30

Type PolyReg_ANOVA
  Reg_df As Double
  Reg_SS As Double
  Residual_df As Double
  Residual_SS As Double
  Total_df As Double
  Total_SS As Double
  S2 As Double
  F As Double
End Type

Type PolyRegrec
  NumTerms As Integer
  ' summation block
  sum As Double
  sumY As Double
  sumYY As Double
  hasMissingData As Integer
  missingCode As Double
  Inverted_S As Integer
  Intercept As Double
  R2 As Double
  R2adj As Double
  ANOVA As PolyReg_ANOVA
End Type

Sub CalcPolyReg (r As PolyRegrec, A() As Double, B() As Double, S() As Double, Mean() As Double, StdErSlope() As Double, TCalc() As Double)

  Dim Diag As Double, tempo As Double
  Dim DegF As Double, MS As Double
  Dim i As Integer, j As Integer
  Dim n As Integer, k As Integer

  n = r.NumTerms + 1
  r.sum = A(0, 0)
  r.sumY = B(0)
  ' Form the A and S Matrices
  For i = 1 To n - 1
    A(i, 0) = B(i) - r.sumY * A(0, i) / r.sum
    For j = 1 To i
      A(j, i) = A(j, i) - A(0, i) * A(0, j) / r.sum
      ' Make a copy of the matrix A
      S(j, i) = A(j, i)
    Next j
  Next i

  ' Clear the matrix S inversion flag
  r.Inverted_S = False

  r.sumYY = r.sumYY - r.sumY ^ 2 / r.sum
  r.sumY = r.sumY / r.sum
  For i = 1 To n - 1
    A(0, i) = A(0, i) / r.sum
    ' Make copies of the mean and std error vectors
    Mean(i) = A(0, i)
    StdErSlope(i) = A(i, i)
  Next i

  ' Form correlation matrix
  For i = 1 To n - 1
    A(i, 0) = A(i, 0) / Sqr(r.sumYY * StdErSlope(i))
    For j = 1 To i - 1
      A(j, i) = A(j, i) / Sqr(StdErSlope(i) * StdErSlope(j))
    Next j
    A(i, i) = 1#
  Next i

  ' Complete symmetric part of the matrix
  For i = 2 To n - 1
    For j = 1 To i - 1
      A(i, j) = A(j, i)
    Next j
  Next i
  
  ' Start solving the simultaneous equations
  For j = 1 To n - 1
    Diag = A(j, j)
    A(j, j) = 1
    For k = 1 To n - 1
      A(j, k) = A(j, k) / Diag
    Next k

    For k = 1 To n - 1
      If j <> k Then
        tempo = A(k, j)
        A(k, j) = 0
        For i = 1 To n - 1
          A(k, i) = A(k, i) - A(j, i) * tempo
        Next i
      End If
    Next k
  Next j
  r.R2 = 0
  r.Intercept = 0
  For i = 1 To n - 1
    B(i) = 0
    For j = 1 To n - 1
      B(i) = B(i) + A(j, 0) * A(j, i)
    Next j
    r.R2 = r.R2 + B(i) * A(i, 0)
    B(i) = B(i) * Sqr(r.sumYY / StdErSlope(i))
    r.Intercept = r.Intercept + B(i) * A(0, i)
  Next i
  r.Intercept = r.sumY - r.Intercept
  DegF = r.sum - r.NumTerms - 1
  MS = (1 - r.R2) * r.sumYY / DegF
  For i = 1 To n - 1
    StdErSlope(i) = Sqr(MS * A(i, i) / StdErSlope(i))
    TCalc(i) = B(i) / StdErSlope(i)
  Next i
  r.R2adj = 1 - (1 - r.R2) * r.sum / DegF
  r.ANOVA.Total_df = r.NumTerms + 1
  r.ANOVA.Total_SS = r.sumYY
  r.ANOVA.Reg_df = r.NumTerms
  r.ANOVA.S2 = MS
  r.ANOVA.Reg_SS = r.sumYY - r.ANOVA.S2
  r.ANOVA.Residual_df = DegF
  r.ANOVA.Residual_SS = r.ANOVA.S2 * DegF
  r.ANOVA.F = DegF / r.NumTerms * r.R2 / (1 - r.R2)
End Sub

Sub CalcSumPolyReg (DataMat() As Double, r As PolyRegrec, NData As Integer, A() As Double, B() As Double, Index() As Integer)

  Dim i As Integer, j As Integer
  Dim n As Integer, k As Integer
  Static T(2) As Double
  Dim Yt As Double
  Dim Xt As Double
  Dim OK As Integer

  n = r.NumTerms + 1

  For i = 0 To NData - 1
    T(0) = DataMat(i, Index(0))
    T(1) = DataMat(i, Index(1))

    ' set ok flag
    OK = True

    ' are there any possible missing data?
    If r.hasMissingData Then
      ' search for missing data
      j = 0
      Do While (j < 2) And OK
       If T(j) > r.missingCode Then
         j = j + 1
       Else
         OK = False
       End If
     Loop
    End If

    If OK Then
       ' now transform the data
       T(0) = MyPolyFx(T(0))
       T(1) = MyPolyFy(T(1))

       Yt = T(0)
       Xt = T(1)
       r.sumYY = r.sumYY + Yt ^ 2
       For k = 0 To n - 1
         B(k) = B(k) + Yt * Xt ^ k
         For j = 0 To n - 1
           A(j, k) = A(j, k) + Xt ^ (j + k)
         Next j
       Next k
    End If
  Next i
End Sub

Sub InitializePolyReg (r As PolyRegrec, NumTerms As Integer, hasMissingData As Integer, missingCode As Double, A() As Double, B() As Double)

  Dim i As Integer
  Dim j As Integer
  Dim n As Integer

  r.NumTerms = NumTerms
  n = r.NumTerms + 1

  For i = 0 To n - 1
    ' initialize summations
    For j = 0 To n - 1
      A(i, j) = 0
    Next j
    B(i) = 0
  Next i
  r.sumYY = 0
  r.hasMissingData = hasMissingData
  r.missingCode = missingCode
End Sub

Sub PolyReg_R2_T_Test (r As PolyRegrec, probability As Double, calcT As Double, tableT As Double, passTest As Integer)
' test hypothesis H0 : R^2 = 0

  Dim p As Double, df As Double

  If probability > 1 Then
    p = .5 - probability / 200
  Else
    p = .5 - probability / 2
  End If
  df = r.sum - (r.NumTerms + 1)
  tableT = TInv(p, df)
  calcT = Sqr(r.R2 * df / (1 - r.R2))
  passTest = calcT <= tableT
End Sub

Sub PolyReg_Slope_T_test (r As PolyRegrec, probability As Double, testValue As Double, termNum As Integer, calcT As Double, tableT As Double, passTest As Integer, B() As Double, StdErSlope() As Double)
  Dim p As Double, df As Double

  If probability > 1 Then
    p = .5 - probability / 200
  Else
    p = .5 - probability / 2
  End If
  df = r.sum - (r.NumTerms + 1)
  tableT = TInv(p, df)
  calcT = (B(termNum) - testValue) / StdErSlope(termNum)
  passTest = Abs(calcT) <= tableT
End Sub

Sub PolyRegCoefCI (r As PolyRegrec, probability As Double, slopeHi() As Double, slopeLow() As Double, B() As Double, StdErSlope() As Double)

  Dim diff As Double, p As Double
  Dim df As Double, tableT As Double
  Dim j As Integer
  Dim n As Integer

  n = r.NumTerms + 1
  If probability > 1 Then
    p = .5 - probability / 200
  Else
    p = .5 - probability / 2
  End If
  df = r.sum - (r.NumTerms + 1)
  tableT = TInv(p, df)
  For j = 1 To n - 1
   diff = tableT * StdErSlope(j)
   slopeHi(j) = B(j) + diff
   slopeLow(j) = B(j) - diff
  Next j
End Sub

Sub PolyYHatCI (r As PolyRegrec, probability As Double, X As Double, yHatLow As Double, yHat As Double, yHatHigh As Double, B() As Double, S() As Double, Mean() As Double)

  Dim diff As Double, p As Double
  Dim df As Double, tableT As Double
  Dim traceMat As Double, pivot As Double
  Dim tempo As Double
  Dim i As Integer, j As Integer
  Dim k As Integer, m As Integer
  Dim XX() As Double
  Dim prodMat() As Double
  Dim Xpow() As Double
  Dim n As Integer
  
  n = r.NumTerms + 1
  ReDim XX(n, n)
  ReDim prodMat(n, n)
  ReDim Xpow(n)
 
  If probability > 1 Then
    p = .5 - probability / 200
  Else
    p = .5 - probability / 2
  End If
  df = r.sum - (r.NumTerms + 1)
  tableT = TInv(p, df)

  If Not r.Inverted_S Then ' Invert matrix S
    r.Inverted_S = True
    For j = 1 To n - 1
      pivot = S(j, j)
      S(j, j) = 1
      For k = 1 To n - 1
        S(j, k) = S(j, k) / pivot
      Next k
      For k = 1 To n - 1
        If k <> j Then
          tempo = S(k, j)
          S(k, j) = 0
          For m = 1 To n - 1
            S(k, m) = S(k, m) - S(j, m) * tempo
          Next m
        End If
     Next k
   Next j
  End If
  
  ' Calculate yHat
  yHat = r.Intercept
  X = MyPolyFx(X)
  For i = 1 To n - 1
    Xpow(i) = X ^ i
    yHat = yHat + B(i) * Xpow(i)
  Next i

  ' Form standarized vector
  For i = 1 To n - 1
    Xpow(i) = Xpow(i) - Mean(i)
  Next i

  ' Form X X' matrix
  For k = 1 To n - 1
    For j = 1 To n - 1
      XX(j, k) = 0
    Next j
  Next k
    
  For k = 1 To n - 1
    For j = 1 To n - 1
      XX(j, k) = XX(j, k) + Xpow(j) * Xpow(k)
    Next j
  Next k

  ' Multiply S_Inverse and XX'
  For i = 1 To n - 1
    For j = 1 To n - 1
      prodMat(i, j) = 0
      For k = 1 To n - 1
        prodMat(i, j) = prodMat(i, j) + S(i, k) * XX(k, j)
      Next k
    Next j
  Next i

  ' Calculate trace of prodMat
  traceMat = 1#
  For i = 1 To n - 1
    traceMat = traceMat * prodMat(i, i)
  Next i

  diff = tableT * Sqr(r.ANOVA.S2 * (1 + 1 / r.sum + traceMat))
  yHatLow = yHat - diff
  yHatHigh = yHat + diff

  yHat = MyPolyInvFy(yHat)
  yHatLow = MyPolyInvFy(yHatLow)
  yHatHigh = MyPolyInvFy(yHatHigh)
End Sub

