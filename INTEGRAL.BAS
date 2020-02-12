Global Const INTEGRAL_BAD_RESULT# = -1E+30
Global Const MAX_ROMBERG_TABLE% = 52




Function AltExtSimpson (numElem As Integer, yArr() As Double, h As Double) As Double

  Dim sum As Double
  Dim i As Integer, isOdd As Integer

  isOdd = (numElem Mod 2) <> 0
  
  If numElem < 10 Then
    AltExtSimpson = INTEGRAL_BAD_RESULT
    Exit Function
  End If

  ' decrement number of element if value is even
  If Not isOdd Then numElem = numElem - 1

  ' add terms
  sum = 0
  For i = 0 To numElem - 1
    sum = sum + yArr(i)
  Next i

  ' calculate integral
  sum = h / 48 * (48 * sum - 31 * yArr(0) + 11 * yArr(1) - 5 * yArr(2) + yArr(3) + yArr(numElem - 4) - 5 * yArr(numElem - 3) + 11 * yArr(numElem - 2) - 31 * yArr(numElem - 1))

  ' adjust for even number of observations
  If Not isOdd Then
    numElem = numElem + 1
    sum = sum + h / 2 * (yArr(numElem - 2) + yArr(numElem - 1))
  End If
  AltExtSimpson = sum
End Function

Function fAltExtSimpson (xFirst As Double, xLast As Double, n As Integer) As Double

  Dim h As Double
  Dim x As Double
  Dim sum As Double
  Dim i As Integer, isOdd As Integer

  isOdd = (n Mod 2) <> 0

  If n < 10 Then
    fAltExtSimpson = INTEGRAL_BAD_RESULT
    Exit Function
  End If

  ' increment n if it is even
  If Not isOdd Then n = n + 1

  ' calculate increment
  h = (xLast - xFirst) / (n - 1)

  ' add terms
  x = xFirst
  sum = 0
  For i = 0 To n - 1
    sum = sum + MyFx(x)
    x = x + h
  Next i

  ' calculate integral
  sum = h / 48 * (48 * sum - 31 * MyFx(xFirst) + 11 * MyFx(xFirst + h) - 5 * MyFx(xFirst + 2 * h) + MyFx(xFirst + 3 * h) + MyFx(xLast - 3 * h) - 5 * MyFx(xLast - 2 * h) + 11 * MyFx(xLast - h) - 31 * MyFx(xLast))

  fAltExtSimpson = sum
End Function

Function fSimpson (xFirst As Double, xLast As Double, n As Integer) As Double

  Dim h As Double
  Dim x As Double
  Dim sumEven As Double
  Dim sumOdd As Double
  Dim sum As Double
  Dim i As Integer, isOdd As Integer

  isOdd = (n Mod 2) <> 0

  If n < 3 Then
    fSimpson = INTEGRAL_BAD_RESULT
    Exit Function
  End If

  ' increment n if it is even
  If Not isOdd Then n = n - 1

  ' calculate increment
  h = (xLast - xFirst) / (n - 1)

  ' add even terms
  x = xFirst + h
  sumEven = 0
  For i = 1 To n - 2 Step 2
    sumEven = sumEven + MyFx(x)
    x = x + 2 * h
  Next i

  ' add odd terms
  x = xFirst + 2 * h
  sumOdd = 0
  For i = 2 To n - 3 Step 2
    sumOdd = sumOdd + MyFx(x)
    x = x + 2 * h
  Next i

  ' calculate integral
  sum = h / 3 * (MyFx(xFirst) + 4 * sumEven + 2 * sumOdd + MyFx(xLast))

  fSimpson = sum
End Function

Function GaussChebyshevQuadrature (n As Integer) As Double

  Dim sum As Double
  Dim x As Double
  Dim pi As Double
  Dim i As Integer

  sum = 0
  pi = 4 * Atn(1)
  For i = 1 To n
    x = Cos((2 * i - 1) * pi / 2 / n)
    sum = sum + MyChebyFx(x)
  Next i
  GaussChebyshevQuadrature = pi / n * sum
End Function

Function GaussHermiteQuadrature () As Double

  Dim sum As Double
  Static Xk(5) As Double
  Static Ak(5) As Double
  Dim i As Integer

  Xk(0) = -2.35060497
  Xk(1) = -1.33584907
  Xk(2) = -.43607741
  Xk(3) = .43607741
  Xk(4) = 1.33584907
  Xk(5) = 2.35060497
  Ak(0) = .00453001
  Ak(1) = .15706732
  Ak(2) = .7246296
  Ak(3) = .7246296
  Ak(4) = .15706732
  Ak(5) = .00453001
  sum = 0
  For i = 0 To 5
    sum = sum + Ak(i) * MyHermiteFx(Xk(i))
  Next i
  GaussHermiteQuadrature = sum
End Function

Function GaussLaguerreQuadrature () As Double

  Dim sum As Double
  Static Xk(5) As Double
  Static Ak(5) As Double
  Dim i As Integer

  Xk(0) = .2228466
  Xk(1) = 1.1889321
  Xk(2) = 2.99273663
  Xk(3) = 5.77514357
  Xk(4) = 9.83746742
  Xk(5) = 15.98287398
  Ak(0) = .45896467
  Ak(1) = .41700083
  Ak(2) = .11337333
  Ak(3) = .0103992
  Ak(4) = .00026102
  Ak(5) = .0000009
  sum = 0
  For i = 0 To 5
    sum = sum + Ak(i) * MyLaguerreFx(Xk(i))
  Next i
  GaussLaguerreQuadrature = sum
End Function

Function GaussLegendreQuadrature (xFirst As Double, xLast As Double, nSubIntervals As Integer) As Double

  Dim xA As Double, xB As Double, xJ As Double
  Dim h As Double, hDiv2 As Double
  Dim sum As Double, area As Double
  Static Xk(5) As Double
  Static Ak(5) As Double
  Dim n As Integer, i As Integer, j As Integer

  Xk(0) = -.9324695142
  Xk(1) = -.6612093865
  Xk(2) = -.2386191861
  Xk(3) = .2386191861
  Xk(4) = .6612093865
  Xk(5) = .9324695142
  Ak(0) = .1713244924
  Ak(1) = .360761573
  Ak(2) = .4679139346
  Ak(3) = .4679139346
  Ak(4) = .360761573
  Ak(5) = .1713244924
  area = 0
  n = nSubIntervals

  h = (xLast - xFirst) / n
  xA = xFirst
  For i = 1 To n
	 sum = 0
	 xB = xA + h
	 hDiv2 = h / 2
	 ' obtain area of sub-interval
	 For j = 0 To 5
		xJ = xA + hDiv2 * (Xk(j) + 1)
		sum = sum + Ak(j) * MyFx(xJ)
	 Next j
	 area = area + hDiv2 * sum
	 xA = xB
  Next i

  GaussLegendreQuadrature = area
End Function

Function Romberg (xA As Double, xB As Double, tolerance As Double, MaxIter As Integer, area As Double) As Integer
  Const goOn% = 1
  Const converge% = 2
  Const reachedMax% = 3

  Static T(MAX_ROMBERG_TABLE) As Double
  Static fnArr(MAX_ROMBERG_TABLE) As Double
  Dim fA As Double, h As Double, c1 As Double
  Dim i As Integer, j As Integer
  Dim nIter As Integer, nSteps As Integer
  Dim opFlag As Integer

  opFlag = goOn
  nIter = 0
  nSteps = 1
  h = xB - xA
  fA = (MyFx(xA) + MyFx(xB)) / 2
  T(0) = h * fA
  area = T(0)

  Do While opFlag = goOn
	 ' divide the step size by 2
	 nIter = nIter + 1
	 h = h / 2
	 nSteps = 2 * nSteps
	 c1 = T(0)

	 ' get new area estimate
	 j = 0
	 i = nSteps - 1
	 Do
		j = j + 1
		fnArr(i - 1) = MyFx(xA + i * h)
		If i > 1 Then
		  fnArr(i - 2) = fnArr(nSteps / 2 - j - 1)
	    End If
		i = i - 2
	 Loop While i >= 1

	 T(0) = fA
	 For i = 1 To nSteps - 1
		T(0) = T(0) + fnArr(i - 1)
	 Next i
	 T(0) = T(0) * h

	 ' perform Richardson's extrapolation
	 For i = 2 To (nIter + 1)
	   T(i - 1) = T(i - 2) + (T(i - 2) - c1) / (4 ^ (i - 1) - 1)
	 Next i

	 ' verify convergence
	 If Abs(T(nIter) - area) <= tolerance Then
		opFlag = converge
	 ElseIf nIter >= MaxIter Then
		opFlag = reachedMax
       End If
	 area = T(nIter)
  Loop
  Romberg = (opFlag = converge)
End Function

Function Simpson (numElem As Integer, yArr() As Double, h As Double) As Double

  Dim sumEven As Double
  Dim sumOdd As Double
  Dim sum As Double
  Dim i As Integer, isOdd As Integer
 
  isOdd = (numElem Mod 2) <> 0

  If numElem < 3 Then
    Simpson = INTEGRAL_BAD_RESULT
    Exit Function
  End If

  ' decrement number of element if value is even
  If Not isOdd Then numElem = numElem - 1

  ' add even terms
  sumEven = 0
  For i = 1 To numElem - 2 Step 2
    sumEven = sumEven + yArr(i)
  Next i

  ' add odd terms
  sumOdd = 0
  For i = 2 To numElem - 3 Step 2
    sumOdd = sumOdd + yArr(i)
  Next i

  ' calculate integral
  sum = h / 3 * (yArr(0) + 4 * sumEven + 2 * sumOdd + yArr(numElem - 1))

  ' adjust for even number of observations
  If Not isOdd Then
    numElem = numElem + 1
    sum = sum + h / 2 * (yArr(numElem - 2) + yArr(numElem - 1))
  End If
  Simpson = sum
End Function

