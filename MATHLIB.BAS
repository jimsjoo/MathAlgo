Function BesselJ (z As Double, v As Integer) As Double

  Const tolerance = .00000001
  Dim sum As Double
  Dim term As Double
  Dim i As Integer

  sum = 0
  i = 0

  Do
     term = (-z * z / 4) ^ i / factorial(i) / gamma(v + i + 1)
     i = i + 1
  Loop While Abs(term) > tolerance

  BeselJ = (z / 2) ^ v * sum
End Function

Function BesselY (z As Double, v As Integer) As Double

  BesselY = (BesselJ(z, v) * Cos(v * pi) - BesselJ(z, -v)) / Sin(v * pi)
End Function

Function beta (z As Double, w As Double) As Double
  beta = gamma(z) * gamma(w) / gamma(z + w)
End Function

Function Chebyshev (x As Double, n As Integer) As Double

  Dim T0 As Double
  Dim T1 As Double
  Dim T2 As Double
  Dim i As Integer

  T0 = 1
  T1 = x
  i = 2

  If n < 0 Then
     Chebyshev = -1
  ElseIf n = 0 Then
     Chebyshev = T0
  ElseIf n = 1 Then
     Chebyshev = T1
  Else
    Do While i < n
      T2 = 2 * x * T1 - 2 * T0
      T0 = T1
      T1 = T2
      i = i + 1
    Loop
    Chebyshev = T2
  End If
End Function

Function Ci (x As Double) As Double

  Const tolerance = .00000001
  Dim sum As Double
  Dim term As Double
  Dim pow As Double
  Dim chs As Double
  Dim n As Integer

  sum = 0
  pow = x * x
  chs = -1
  n = 1

  Do
    term = chs * pow / (2 * n) / factorial(2 * n)
    sum = sum + term
    ' update components of the next term
    n = n + 1
    chs = -chs
    pow = pow * x * x

  Loop While Abs(term) > tolerance

  Ci = .5772156649 + Log(x) + sum
End Function

Function combination (m As Integer, n As Integer) As Double

  Dim factM As Double
  Dim factMN As Double
  Dim factN As Double
  Dim i As Integer

  factM = 1
  factMN = 1
  factN = 1

  If (m <= n) Or (m < 0) Or (n < 0) Then
    combination = -1 ' error code
    Exit Function
  End If

  ' calculate (m - n)!
  For i = 1 To m - n
    factMN = factMN * i
  Next i

  factM = factMN
  ' calculate m!
  For i = m - n + 1 To m
    factM = factM * i
  Next i

  ' calculate n!
  For i = 1 To n
     factN = factN * i
  Next i

  combination = factM / factMN / factN
End Function

Function Ei (x As Double) As Double

  Const tolerance = .00000001
  Dim sum As Double
  Dim term As Double
  Dim pow As Double
  Dim fact As Double
  Dim n As Integer

  sum = 0
  pow = x
  fact = 1
  n = 1

  Do
    term = pow / n / fact
    sum = sum + term
    ' update components of the next term
    n = n + 1
    fact = fact * n
    pow = pow * x

  Loop While Abs(term) > tolerance

  Ei = .5772156649 + Log(x) + sum
End Function

Function erf (x As Double) As Double
' error function

  Dim t As Double
  Dim sum As Double
  Dim i As Integer
  Static a(5)

  t = 1 / (1 + .32759 * x)
  a(0) = 0
  a(1) = .254829592
  a(2) = -.284496736
  a(3) = 1.424143741
  a(4) = -1.453152027
  a(5) = 1.061405429

  sum = 0
  For i = 5 To 0 Step -1
    sum = (sum + a(i)) * t
  Next i

  erf = 1 - sum / Exp(x * x)
End Function

Function factorial (n As Integer) As Double
  Dim result As Double
  Dim i As Integer
  result = 1

  If n > -1 Then
    For i = 1 To n
      result = result * i
    Next i
  End If
  factorial = result
End Function

Function gamma (x As Double) As Double
  
  Dim sum As Double
  Dim i As Integer
  Static c(26) As Double
  c(0) = 0
  c(1) = 1
  c(2) = .5772156649
  c(3) = -.6558780715
  c(4) = -.042002635034
  c(5) = .16665386113
  c(6) = -.0421977345555
  c(7) = -.0096219715278
  c(8) = .00721894324666
  c(9) = -1.1651675918591E-03
  c(10) = -2.152416741149E-04
  c(11) = 1.280502823882E-04
  c(12) = -2.01348547807E-05
  c(13) = -1.2504934821E-06
  c(14) = .000001133027232
  c(15) = -2.056338417E-07
  c(16) = .000000006116095
  c(17) = 5.0020075E-09
  c(18) = -1.1812746E-09
  c(19) = 1.043427E-10
  c(20) = 7.7823E-12
  c(21) = -3.6968E-12
  c(22) = .00000000000051
  c(23) = -2.06E-14
  c(24) = -5.4E-15
  c(25) = 1.4E-16
  c(26) = 1E-16

  sum = 0

  For i = 26 To 0 Step -1
   sum = (sum + c(i)) * x
  Next i

  gamma = 1 / sum
End Function

Function Hermite (x As Double, n As Integer) As Double

  Dim H0 As Double
  Dim H1 As Double
  Dim H2 As Double
  Dim i As Integer

  H0 = 1
  H1 = 2 * x
  i = 2

  If n < 0 Then
    Hermite = -1
  ElseIf n = 0 Then
    Hermite = H0
  ElseIf n = 1 Then
    Hermite = H1
  Else
    Do While i < n
      H2 = 2 * x * H1 - 2 * i * H0
      H0 = H1
      H1 = H2
      i = i + 1
    Loop
    Hermite = H2
  End If
End Function

Function Laguerre (x As Double, n As Integer) As Double

  Dim L0 As Double
  Dim L1 As Double
  Dim L2 As Double
  Dim i As Integer

  L0 = 1
  L1 = 1 - x
  i = 2
  
  If n < 0 Then
    Laguerre = -1
  ElseIf n = 0 Then
    Laguerre = L0
  ElseIf n = 1 Then
    Laguerre = L1
  Else
    Do While i < n
      L2 = (1 + 2 * i - x) * L1 - i * i * L0
      L0 = L1
      L1 = L2
      i = i + 1
    Loop
    Laguerre = L2
  End If
End Function

Function Legrendre (x As Double, n As Double) As Double

  Dim L0 As Double
  Dim L1 As Double
  Dim L2 As Double
  Dim i As Integer

  L0 = 1
  L1 = x
  i = 2

  If n < 0 Then
    Legrendre = -1
  ElseIf n = 0 Then
    Legrendre = L0
  ElseIf n = 1 Then
    Legrendre = L1
  Else
    Do While i < n
      L2 = ((2 * i + 1) * L1 - i * L0) / (i + 1)
      L0 = L1
      L1 = L2
      i = i + 1
    Loop
    Legrendre = L2
  End If
End Function

Function permutation (m As Integer, n As Integer) As Double

  Dim factM As Double
  Dim factMN As Double
  Dim i As Integer

  factM = 1
  factMN = 1

  If (m <= n) Or (m < 0) Or (n < 0) Then
     permutation = -1 ' error code
     Exit Function
  End If

  ' calculate (m - n)!
  For i = 1 To (m - n)
    factMN = factMN * i
  Next i

  factM = factMN
  ' calculate m!
  For i = m - n + 1 To m
    factM = factM * i
  Next i

  permutation = factM / factMN
End Function

Function Si (x As Double) As Double

  Const tolerance = .00000001
  Dim sum As Double
  Dim term As Double
  Dim pow As Double
  Dim chs As Double
  Dim n As Integer

  sum = 0
  pow = x
  chs = 1
  n = 0
   
  Do
    term = chs * pow / (2 * n + 1) / factorial(2 * n + 1)
    sum = sum + term
    ' update components of the next term
    n = n + 1
    chs = -chs
    pow = pow * x * x
  Loop While Abs(term) > tolerance

  Si = sum
End Function

