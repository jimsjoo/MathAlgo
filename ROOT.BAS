
Global Const ROOT_EPS# = 1E-30

Type polyRoot
  real As Double
  imag As Double
  isComplex As Integer
End Type

Function Bisection (low As Double, high As Double, tolerance As Double, root As Double, params() As Double) As Integer

  If MyFx(low, params()) * MyFx(high, params()) > 0 Then
     Bisection = False
     Exit Function
  End If

  Do While Abs(high - low) > tolerance
    ' update guess
    root = (low + high) / 2
    If (MyFx(root, params()) * MyFx(high, params())) > 0 Then
      high = root
    Else
      low = root
    End If
  Loop
  Bisection = True
End Function

Function Brent (low As Double, high As Double, tolerance As Double, root As Double, maxIter As Integer, params() As Double) As Integer

  Const SMALL# = .0000001' epsilon
  Const VERY_SMALL# = .0000000001' near zero
  Dim iter As Integer
  Dim a As Double, b As Double, c As Double
  Dim fa As Double, fb As Double, fc As Double
  Dim d As Double, e As Double, tol As Double
  Dim small1 As Double, small2 As Double, small3 As Double
  Dim p As Double, q As Double, r As Double
  Dim s As Double, xm As Double

  iter = 1
  a = low
  b = high
  c = high
  fa = MyFx(low, params())
  fb = MyFx(high, params())
  fc = fb

  ' check that the guesses contain the root
  If (fa * fb) > 0 Then
     Brent = False ' bad guesses
     Exit Function
  End If

  ' start loop to refine the guess for the root
  Do While iter <= maxIter
    iter = iter + 1
    If (fb * fc) > 0 Then
      c = a
      fc = fa
      e = b - a
      d = e
    End If
    If Abs(fc) < Abs(fb) Then
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
    End If
    tol = 2 * SMALL * Abs(b) + tolerance / 2
    xm = (c - b) / 2
    If (Abs(xm) <= tol) Or (Abs(fb) <= VERY_SMALL) Then
      root = b
      Brent = True
      Exit Function
    End If
    If (Abs(e) >= tol) And (Abs(fa) > Abs(fb)) Then
      ' perform the inverse quadratic interpolation
      s = fb / fa
      If Abs(a - c) <= VERY_SMALL Then
        p = 2 * xm * s
        q = 1 - s
      Else
        q = fa / fc
        r = fb / fc
        q = s * (2 * xm * q * (q - r) - (b - a) * (r - 1))
        q = (q - 1) * (r - 1) * (s - 1)
      End If
      ' determine if improved guess is inside
      ' the range
      If p > 0 Then q = -q
      p = Abs(p)
      small1 = 3 * xm * q * Abs(tol * q)
      small2 = Abs(e * q)
      If small1 < small2 Then
        small3 = small1
      Else
        small3 = small2
      End If
      If (2 * p) < small3 Then
        ' interpolation successfull
        e = d
        d = p / q
      Else
        ' use bisection because the interpolation
        ' did not succeed
        d = xm
        e = d
      End If
    Else
      ' use bisection because the range
      ' is slowly decreasing
      d = xm
      e = d
    End If
    ' copy most recent guess to variable a
    a = b
    fa = fb
    ' evaluate improved guess for the root
    If Abs(d) > tol Then
      b = b + d
    Else
      If xm > 0 Then
        b = b + Abs(tol)
      Else
        b = b - Abs(tol)
      End If
    End If
    fb = MyFx(b, params())
  Loop
  Brent = False
End Function

Function Combined (low As Double, high As Double, tolerance As Double, root As Double, maxIter As Integer, params() As Double) As Integer

  Dim iter As Integer
  Dim h As Double, diff As Double

  iter = 0
  Do
     iter = iter + 1
     h = .01 * root
     If Abs(root) < 1 Then h = .01
     ' calculate guess refinement
     diff = 2 * h * MyFx(root, params()) / (MyFx(root + h, params()) - MyFx(root - h, params()))
     root = root - diff
     ' check if Newton's method yields a refined guess
     ' outside the range (low, high)
     If (root < low) Or (root > high) Then
       ' apply Bisection method for this iteration
       root = (low + high) / 2
       If (MyFx(root, params()) * MyFx(high, params())) > 0 Then
         high = root
       Else
         low = root
       End If
     End If
  Loop While (iter <= maxIter) And (Abs(diff) > tolerance)

  If Abs(diff) <= tolerance Then
    Combined = True
  Else
    Combined = False
  End If
End Function

Function DeflatePolyRoots (coeff() As Double, initGuess As Double, roots() As Double, numRoots As Integer, polyOrder As Integer, maxIter As Integer, tolerance As Double) As Integer

  Dim a() As Double
  Dim b() As Double
  Dim c() As Double
  Dim diff As Double
  Dim z As Double, X As Double
  Dim i As Integer, iter As Integer, n As Integer

  iter = 1
  n = polyOrder + 1
  X = initGuess
  numRoots = 0

  ' allocate dynamic coefficients
  ReDim a(n)
  ReDim b(n)
  ReDim c(n)

  For i = 0 To n - 1
     a(i) = coeff(i)
  Next i

  Do While (iter <= maxIter) And (n > 1)
     iter = iter + 1
     z = X
     b(n - 1) = a(n - 1)
     c(n - 1) = a(n - 1)
     For i = n - 2 To 1 Step -1
       b(i) = a(i) + z * b(i + 1)
       c(i) = b(i) + z * c(i + 1)
     Next i
     b(0) = a(0) + z * b(1)
     diff = b(0) / c(1)
     X = X - diff
     If Abs(diff) <= tolerance Then
       iter = 1 ' reset iteration counter
       n = n - 1
       roots(numRoots) = X
       numRoots = numRoots + 1
       ' update deflated roots
       For i = 0 To n - 1
         a(i) = b(i + 1)
         ' get the last root
         If n = 2 Then
           n = n - 1
           roots(numRoots) = -a(0)
           numRoots = numRoots + 1
         End If
       Next i
     End If
  Loop
  DeflatePolyRoots = True
End Function

Function LBPolyRoots (coeff() As Double, roots() As polyRoot, polyOrder As Integer, tolerance As Double) As Integer
'
'  solves for the roots of the following polynomial
'
'     y = coeff(0) + coeff(1) X + coeff(2) X^2 +...+ coeff(n) X^n
'
' Parameters:
'
'   coeff      must be an array with at least polyOrder+1 elements.
'
'   roots      output array of roots
'
'   polyOrder  order of polynomial
'
'   tolerance  tolerance of solutions
'
'

  Const SMALL# = .00000001
  Dim a() As Double
  Dim b() As Double
  Dim c() As Double
  Dim d() As Double
  Dim alfa1 As Double, alfa2 As Double
  Dim beta1 As Double, beta2 As Double
  Dim delta1 As Double, delta2 As Double, delta3 As Double
  Dim i As Integer, j As Integer, k As Integer
  Dim count As Integer
  Dim n As Integer
  Dim n1 As Integer
  Dim n2 As Integer

  n = polyOrder
  n1 = n + 1
  n2 = n + 2
  ' is the coefficient of the highest term zero?
  If Abs(coeff(0)) < SMALL Then
    LBPolyRoots = False
    Exit Function
  End If

  ' allocate dynamic coefficients
  ReDim a(n1)
  ReDim b(n1)
  ReDim c(n1)
  ReDim d(n1)

  For i = 0 To n
     a(n1 - i) = coeff(i)
  Next i

  ' is highest coeff not close to 1?
  If Abs(a(1) - 1) > SMALL Then
     ' adjust coefficients because a(1) != 1
     For i = 2 To n1
       a(i) = a(i) / a(1)
     Next i
     a(1) = 1#
  End If

  ' initialize root counter
  count = 0

  Do
     '
     ' start the main Lin-Bairstow iteration loop
     ' initialize the counter and guesses for the
     ' coefficients of quadratic factor:
     '
     ' p(x) = x^2 + alfa1 * x + beta1
     '

     alfa1 = 1
     beta1 = 1

     Do
       b(0) = 0
       d(0) = 0
       b(1) = 1
       d(1) = 1

       j = 1
       k = 0
       For i = 2 To n1
         b(i) = a(i) - alfa1 * b(j) - beta1 * b(k)
         d(i) = b(i) - alfa1 * d(j) - beta1 * d(k)
         j = j + 1
         k = k + 1
       Next i
       j = n - 1
       k = n - 2
       delta1 = d(j) ^ 2 - (d(n) - b(n)) * d(k)
       alfa2 = (b(n) * d(j) - b(n1) * d(k)) / delta1
       beta2 = (b(n1) * d(j) - (d(n) - b(n)) * b(n)) / delta1
       alfa1 = alfa1 + alfa2
       beta1 = beta1 + beta2
     Loop While (Abs(alfa2) > tolerance) Or (Abs(beta2) > tolerance)

     delta1 = alfa1 ^ 2 - 4 * beta1

     If delta1 < 0 Then
       ' imaginary roots
       delta2 = Sqr(Abs(delta1)) / 2
       delta3 = -alfa1 / 2
       For i = 0 To 1
         roots(count + i).isComplex = True
         roots(count + i).real = delta3
         roots(count + i).imag = Sign(i) * delta2
       Next i
     Else
       delta2 = Sqr(delta1)
       ' roots are real
       For i = 0 To 1
         roots(count + i).isComplex = False
         roots(count + i).imag = 0
       Next i
       roots(count).real = (delta2 - alfa1) / 2
       roots(count + 1).real = (delta2 + alfa1) / (-2)
     End If
     ' update root counter
     count = count + 2

     ' reduce polynomial order
     n = n - 2
     n1 = n1 - 2
     n2 = n2 - 2

     ' for n >= 2 calculate coefficients of
     '  the new polynomial
     If n >= 2 Then
       For i = 1 To n1
         a(i) = b(i)
       Next i
     End If
  Loop While n >= 2

  If n = 1 Then ' obtain last single real root
    roots(count).isComplex = False
    roots(count).real = -b(2)
    roots(count).imag = 0
  End If

  LBPolyRoots = True
End Function

Function Newton2Functions (rootX As Double, rootY As Double, tolerance As Double, maxIter As Integer, paramsX() As Double, paramsY() As Double) As Integer

  Dim Jacob As Double
  Dim fx0 As Double, fy0 As Double
  Dim hX As Double, hY As Double
  Dim diffX As Double, diffY As Double
  Dim fxy As Double, fxx As Double
  Dim fyy As Double, fyx As Double
  Dim X As Double
  Dim Y As Double
  Dim iter As Integer

  X = rootX
  Y = rootY
  iter = 1

  Do
     hX = .01 * rootX
     If Abs(rootX) < 1 Then hX = .01
     hY = .01 * rootY
     If Abs(rootY) < 1 Then hY = .01
     fx0 = My2Fx(X, Y, paramsX())
     fy0 = My2Fy(X, Y, paramsY())
     fxx = (My2Fx(X + hX, Y, paramsX()) - My2Fx(X - hX, Y, paramsX())) / 2 / hX
     fyx = (My2Fy(X + hX, Y, paramsY()) - My2Fy(X - hX, Y, paramsY())) / 2 / hX
     fxy = (My2Fx(X, Y + hY, paramsX()) - My2Fx(X, Y - hY, paramsX())) / 2 / hY
     fyy = (My2Fy(X, Y + hY, paramsY()) - My2Fy(X, Y - hY, paramsY())) / 2 / hY
     Jacob = fxx * fyy - fxy * fyx
     diffX = (fx0 * fyy - fy0 * fxy) / Jacob
     diffY = (fy0 * fxx - fx0 * fyx) / Jacob
     X = X - diffX
     Y = Y - diffY
     iter = iter + 1
  Loop While (iter <= maxIter) And ((Abs(diffX) > tolerance) Or (Abs(diffY) > tolerance))

  rootX = X
  rootY = Y
  If (Abs(diffX) <= tolerance) And (Abs(diffY) <= tolerance) Then
    Newton2Functions = True
  Else
    Newton2Functions = False
  End If
End Function

Function NewtonApprox (root As Double, tolerance As Double, maxIter As Integer, params() As Double) As Integer

  Dim iter As Integer
  Dim h As Double, diff As Double

  iter = 0
  Do
     h = .01 * root
     If Abs(root) < 1 Then h = .01
     ' calculate guess refinement
     diff = 2 * h * MyFx(root, params()) / (MyFx(root + h, params()) - MyFx(root - h, params()))
     ' update guess
     root = root - diff
     iter = iter + 1
  Loop While (iter <= maxIter) And (Abs(diff) > tolerance)

  If Abs(diff) <= tolerance Then
    NewtonApprox = True
  Else
    NewtonApprox = False
  End If
End Function

Function NewtonExact (root As Double, tolerance As Double, maxIter As Integer, params() As Double) As Integer

  Dim iter As Integer
  Dim diff As Double

  Do
     ' calculate guess refinement
     diff = MyFx(root, params()) / MyDeriv1(root, params())
     ' update guess
     root = root - diff
     iter = iter + 1
  Loop While (iter <= maxIter) And (Abs(diff) > tolerance)

  If Abs(diff) <= tolerance Then
    NewtonExact = True
  Else
    NewtonExact = False
  End If
End Function

Function NewtonMultiRoots (roots() As Double, numRoots As Integer, maxRoots As Integer, tolerance As Double, maxIter As Integer, params() As Double) As Integer

  Dim iter As Integer, i As Integer
  Dim h As Double, diff As Double
  Dim f1 As Double, f2 As Double, f3 As Double
  Dim root As Double

  numRoots = 0
  root = roots(0)

  Do
    iter = 0
    Do
      h = .01 * root
      If Abs(root) < 1 Then h = .01
      f1 = MyFx(root - h, params())
      f2 = MyFx(root, params())
      f3 = MyFx(root + h, params())
      If numRoots > 0 Then
        For i = 0 To numRoots - 1
          f1 = f1 / (root - h - roots(i))
          f2 = f2 / (root - roots(i))
          f3 = f3 / (root + h - roots(i))
        Next i
      End If
      ' calculate guess refinement
      diff = 2 * h * f2 / (f3 - f1)
      ' update guess
      root = root - diff
      iter = iter + 1
    Loop While (iter <= maxIter) And (Abs(diff) > tolerance)

    If Abs(diff) <= tolerance Then
      roots(numRoots) = root
      numRoots = numRoots + 1
      If root < 0 Then
        root = .95 * root
      ElseIf root > 0 Then
        root = 1.05 * root
      Else
        root = .05
      End If
    End If
  Loop While (iter <= maxIter) And (numRoots < maxRoots)

  If numRoots > 0 Then
    NewtonMultiRoots = True
  Else
    NewtonMultiRoots = False
  End If
End Function

Function NewtonSimNLE (X() As Double, numEqns As Integer, tolerance As Double, maxIter As Integer) As Integer

  Dim Xdash() As Double
  Dim Fvector() As Double
  Dim index() As Integer
  Dim Jmat() As Double
  Dim i As Integer, j As Integer
  Dim moreIter As Integer, iter As Integer
  Dim rowSwapFlag As Integer
  Dim h As Double
  Dim dummy As Integer

  iter = 0
  ReDim Jmat(numEqns, numEqns)
  ReDim Xdash(numEqns)
  ReDim Fvector(numEqns)
  ReDim index(numEqns)

  Do
    iter = iter + 1
    ' copy the values of array X into array Xdash
    For i = 0 To numEqns - 1
      Xdash(i) = X(i)
    Next i
    ' calculate the array of function values
    For i = 0 To numEqns - 1
      Fvector(i) = MySNE(X(), i)
    Next i
    ' calculate the Jmat matrix
    For i = 0 To numEqns - 1
      For j = 0 To numEqns - 1
        ' calculate increment in variable number j
        If Abs(X(j)) > 1 Then
          h = .01 * X(j)
        Else
          h = .01
        End If
        Xdash(j) = Xdash(j) + h
        Jmat(i, j) = (MySNE(Xdash(), i) - Fvector(i)) / h
        ' restore incremented value
        Xdash(j) = X(j)
      Next j
    Next i
    ' solve for the guess refinement vector
    dummy = LUDecomp(Jmat(), index(), numEqns, rowSwapFlag)
    LUBackSubst Jmat(), index(), numEqns, Fvector()
    ' clear the more-iteration flag
    moreIter = False
    ' update guess and test for convergence
    For i = 0 To numEqns - 1
      X(i) = X(i) - Fvector(i)
      If Abs(Fvector(i)) > tolerance Then moreIter = True
    Next i
    ' check iteration limit
    If moreIter Then
      If iter > maxIter Then
        moreIter = False
      Else
        moreIter = True
      End If
    End If
  Loop While moreIter

  NewtonSimNLE = Not moreIter
End Function

Function RichmondApprox (root As Double, tolerance As Double, maxIter As Integer, params() As Double) As Integer

  Dim iter As Integer
  Dim h As Double, diff As Double
  Dim f1 As Double, f2 As Double, f3 As Double
  Dim fd1 As Double, fd2 As Double

  Do
     h = .01 * root
     If Abs(root) < 1 Then h = .01
     f1 = MyFx(root - h, params())      ' f(x-h)
     f2 = MyFx(root, params())          ' f(x)
     f3 = MyFx(root + h, params())      ' f(x+h)
     fd1 = (f3 - f1) / (2 * h)          ' f'(x)
     fd2 = (f3 - 2 * f2 + f1) / Sqr(h)  ' f''(x)
     ' calculate guess refinement
     diff = f1 * fd1 / (fd1 ^ 2 - .5 * f1 * fd2)
     ' update guess
     root = root - diff
     iter = iter + 1
  Loop While (iter <= maxIter) And (Abs(diff) > tolerance)

  If Abs(diff) <= tolerance Then
    RichmondApprox = True
  Else
    RichmondApprox = False
  End If
End Function

Function RichmondExact (root As Double, tolerance As Double, maxIter As Integer, params() As Double) As Integer

  Dim iter As Integer
  Dim diff As Double, f1 As Double
  Dim fd1 As Double, fd2 As Double

  iter = 0
  Do
     f1 = MyFx(root, params())
     fd1 = MyDeriv1(root, params())
     fd2 = MyDeriv2(root, params())
     ' calculate guess refinement
     diff = f1 * fd1 / (fd1 ^ 2 - .5 * f1 * fd2)
     ' update guess
     root = root - diff
     iter = iter + 1
  Loop While (iter <= maxIter) And (Abs(diff) > tolerance)

  If Abs(diff) <= tolerance Then
    RichmondExact = True
  Else
    RichmondExact = False
  End If
End Function

Function Sign (X As Integer) As Double
  If X = 0 Then
    Sign = 1
  Else
    Sign = -1
  End If
End Function

