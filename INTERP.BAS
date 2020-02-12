Const BAD_RESULT# = -1E+30

Function Barycentric (xarr() As Double, yarr() As Double, wt() As Double, calcWtFlag As Integer, n As Integer, x As Double) As Double

  Dim prod As Double, diff As Double
  Dim sum1 As Double, sum2 As Double
  Dim i As Integer, j As Integer
  
  sum1 = 0
  sum2 = 0

  If calcWtFlag Then
    For i = 0 To n - 1
      prod = 1
      For j = 0 To n - 1
	If i <> j Then
	  prod = prod * (xarr(i) - xarr(j))
	End If
      Next j
	wt(i) = 1 / prod
    Next i
  End If

  '  loop for each term
  For i = 0 To n - 1
    diff = wt(i) / (x - xarr(i))
    sum1 = sum1 + yarr(i) * diff
    sum2 = sum2 + diff
  Next i
  Barycentric = sum1 / sum2
End Function

Function ED_Barycentric (xarr() As Double, yarr() As Double, n As Integer, x As Double) As Double

  Dim h As Double, x0 As Double, diff As Double
  Dim w As Double, s As Double, t As Double
  Dim i As Integer, k As Integer, m As Integer

  i = 0
  w = 1
  k = 0
  m = n - 1
  s = 0
  t = 0
  x0 = xarr(0)
  h = xarr(1) - x0

  Do While True
    diff = x - x0
    If Abs(diff) < INTERP_EPS Then
      diff = INTERP_EPS
    End If
    t = t + w / diff
    s = s + w * yarr(i) / diff
    If m <> 0 Then
      w = w * m
      m = m - 1
      k = k - 1
      w = w / k
      x0 = x0 + h
      i = i + 1
    Else
      Exit Do
    End If
  Loop

  ED_Barycentric = s / t
End Function

Function getDerivatives (xarr() As Double, yarr() As Double, n As Integer, deriv() As Double, h() As Double, tolerance As Double) As Integer

  Dim a() As Double
  Dim b() As Double
  Dim c() As Double
  Dim i As Integer, j As Integer
  Dim k As Integer, result As Integer

  ReDim a(n)
  ReDim b(n)
  ReDim c(n)

  For i = 0 To n - 1
    h(i) = 0
    deriv(0) = 0
  Next i

  For i = 1 To n - 1
    j = i - 1
    h(i) = xarr(i) - xarr(j)
    deriv(i) = (yarr(i) - yarr(j)) / h(i)
  Next i

  For i = 1 To n - 2
    j = i - 1
    k = i + 1
    b(j) = 2
    c(j) = h(k) / (h(i) + h(k))
    a(j) = 1 - c(j)
    deriv(j) = 6 * (deriv(k) - deriv(i)) / (h(i) + h(k))
  Next i

  ' solve for the second derivatives
  getDerivatives = tridiagonal(n - 2, a(), b(), c(), deriv(), tolerance)
End Function

Function Lagrange (xarr() As Double, yarr() As Double, n As Integer, x As Double) As Double

  Dim prod As Double, yint As Double
  Dim i As Integer, j As Integer
  
  yint = 0

  '  loop for each term
  For i = 0 To n - 1
    ' initialize term with yint[i]
    prod = yarr(i)
    ' build each term
    For j = 0 To n - 1
      If i <> j Then
	prod = prod * (x - xarr(j)) / (xarr(i) - xarr(j))
      End If
    Next j
    yint = yint + prod
  Next i
  Lagrange = yint
End Function

Function NewtonDiff (x0 As Double, h As Double, yarr() As Double, diffTable() As Double, buildMatFlag As Integer, n As Integer, x As Double) As Double

  Dim yint As Double
  Dim i As Integer, j As Integer

  If buildMatFlag Then
    For i = 0 To n - 1
      diffTable(i) = yarr(i)
    Next i

    For i = 0 To n - 2
      For j = n - 1 To i + 1 Step -1
	diffTable(j) = (diffTable(j - 1) - diffTable(j)) / (-h * (1 + i))
      Next j
    Next i
  End If

  ' initialize interpolated value
  yint = diffTable(n - 1)

  For i = n - 2 To 0 Step -1
    yint = diffTable(i) + (x - (x0 + i * h)) * yint
  Next i

  NewtonDiff = yint
End Function

Function NewtonDivDiff (xarr() As Double, yarr() As Double, diffTable() As Double, buildMatFlag As Integer, n As Integer, x As Double) As Double

  Dim yint As Double
  Dim i As Integer, j As Integer

  If buildMatFlag Then
    For i = 0 To n - 1
      diffTable(i) = yarr(i)
    Next i
    For i = 0 To n - 2
      For j = n - 1 To i + 1 Step -1
	diffTable(j) = (diffTable(j - 1) - diffTable(j)) / (xarr(j - 1 - i) - xarr(j))
      Next j
    Next i
  End If

  ' initialize interpolated value
  yint = diffTable(n - 1)

  For i = n - 2 To 0 Step -1
    yint = diffTable(i) + (x - xarr(i)) * yint
  Next i

  NewtonDivDiff = yint
End Function

Function Spline (xarr() As Double, yarr() As Double, deriv() As Double, h() As Double, n As Integer, x As Double) As Double

  Dim j As Integer, i As Integer
  Dim delta1 As Double, delta2 As Double
  Dim y As Double

  i = 0
  ' locate subinterval containing x
  Do While ((x < xarr(i)) Or (x > xarr(i + 1))) And (i < (n - 1))
    i = i + 1
  Loop

  If (i < 1) Or (i >= (n - 1)) Then
    Spline = BAD_RESULT
    Exit Function
  End If

  j = i + 1
  delta1 = x - xarr(i)
  delta2 = xarr(j) - x
  y = (deriv(i - 1) * delta2 ^ 3) / (6 * h(j)) + (deriv(i) * delta1 ^ 3) / (6 * h(j)) + (yarr(j) / h(j) - deriv(i) * h(j) / 6) * delta1 + (yarr(i) / h(j) - deriv(i - 1) * h(j) / 6) * delta2

  Spline = y
End Function

Private Function tridiagonal (n As Integer, a() As Double, b() As Double, c() As Double, d() As Double, epsilon As Double) As Integer

  Dim isSingular As Integer
  Dim i As Integer

  isSingular = b(0) < epsilon

  ' carry out LU factorization
  For i = 1 To n - 1
    If isSingular Then Exit For
    a(i) = a(i) / b(i - 1)
    b(i) = b(i) - a(i) * c(i - 1)
    ' determine if diagonal element is too small
    isSingular = b(i) < epsilon
    d(i) = d(i) - a(i) * d(i - 1)
  Next i

  If Not isSingular Then
    ' carry out backward substitution
    d(n - 1) = d(n - 1) / b(n - 1)
    For i = n - 2 To 0 Step -1
      d(i) = (d(i) - c(i) * d(i + 1)) / b(i)
    Next i
    tridiagonal = True
  Else
    tridiagonal = False
  End If
End Function

