Global Const MATVECT_EPSILON# = .000000000000001
Global Const MATVECT_BAD_RESULT# = -1E+30
Global Const matErr_None% = 0
Global Const matErr_Size% = 1
Global Const matErr_Singular% = 2
Global Const matErr_IllConditioned% = 3
Global Const matErr_IterLimit% = 4

Function AddMat (MatC() As Double, MatA() As Double, MatB() As Double, numRows As Integer, numCols As Integer) As Integer

  Dim row As Integer, col As Integer

  If Not (checkRowCol(MatA(), numRows, numCols) And checkRowCol(MatB(), numRows, numCols) And checkRowCol(MatC(), numRows, numCols) And (numCols = numRows)) Then
    AddMat = matErr_Size
    Exit Function
  End If

  For row = 0 To numRows - 1
    For col = 0 To numCols - 1
      MatC(row, col) = MatA(row, col) + MatB(row, col)
    Next col
  Next row
  AddMat = matErr_None
End Function

Function checkRowCol (Mat() As Double, row As Integer, col As Integer) As Integer
  If (row >= 0) And (col >= 0) And (row < UBound(Mat, 1)) And (col < UBound(Mat, 2)) Then
    checkRowCol = True
  Else
    checkRowCol = False
  End If
End Function

Function CopyMat (MatB() As Double, MatA() As Double, numRows As Integer, numCols As Integer) As Integer

  Dim row As Integer, col As Integer

  If Not (checkRowCol(MatA(), numRows, numCols) And checkRowCol(MatB(), numRows, numCols)) Then
      CopyMat = matErr_Size
      Exit Function
  End If

  For row = 0 To numRows - 1
    For col = 0 To numCols
      MatB(row, col) = MatA(row, col)
    Next col
  Next row
  CopyMat = matErr_None
End Function

Function GaussJordan (A() As Double, B() As Double, numRows As Integer, numCols As Integer) As Integer

  Dim rowIndex() As Integer
  Dim colIndex() As Integer
  Dim pivotIndex() As Integer
  Dim i As Integer, j As Integer, k As Integer
  Dim n As Integer, m As Integer
  Dim row As Integer, col As Integer
  Dim large As Double, z As Double, oneOverPiv As Double

  If Not checkRowCol(A(), numRows, numCols) Then
    GaussJordan = matErr_Size
    Exit Function
  End If

  ' resize local dynamic arrays
  ReDim rowIndex(numRows)
  ReDim colIndex(numRows)
  ReDim pivotIndex(numRows)

  ' initialize the row and column indices
  For i = 0 To numRows - 1
    rowIndex(i) = i
    colIndex(i) = i
  Next i

  ' initialize the pivot index array
  For i = 0 To numRows - 1
    pivotIndex(i) = -1
  Next i
  
  For i = 0 To numRows - 1
    large = 0
    For j = 0 To numRows - 1
      If pivotIndex(j) <> 0 Then
	For k = 0 To numRows - 1
	  If pivotIndex(k) = -1 Then
	    If Abs(A(j, k)) >= large Then
	      large = Abs(A(j, k))
	      row = j
	      col = k
	    End If
	  ElseIf pivotIndex(k) > 0 Then
	    GaussJordan = matErr_Singular
	    Exit Function
	  End If
	Next k
      End If
    Next j
    pivotIndex(col) = pivotIndex(col) + 1
    If row <> col Then
      For n = 0 To numRows - 1
	swapDouble A(row, n), A(col, n)
      Next n
      For n = 0 To numCols - 1
	swapDouble B(row, n), B(col, n)
      Next n
    End If
    rowIndex(i) = row
    colIndex(i) = col
    If Abs(A(col, col)) < .0000000001 Then
      GaussJordan = matErr_Singular
      Exit Function
    End If
    oneOverPiv = 1 / A(col, col)
    A(col, col) = 1
    For n = 0 To numRows - 1
      A(col, n) = A(col, n) * oneOverPiv
    Next n
    For n = 0 To numCols - 1
      B(col, n) = B(col, n) * oneOverPiv
    Next n
    For m = 0 To numRows - 1
      If m <> col Then
	z = A(m, col)
	A(m, col) = 1
	For n = 0 To numRows - 1
	  A(m, n) = A(m, n) - A(col, n) * z
	Next n
	For n = 0 To numCols - 1
	  B(m, n) = B(m, n) - B(col, n) * z
	Next n
      End If
    Next m
  Next i
  For n = numRows - 1 To 0 Step -1
    If rowIndex(n) <> colIndex(n) Then
      For k = 0 To numRows - 1
	swapDouble A(k, rowIndex(n)), A(k, colIndex(n))
      Next k
    End If
  Next n
  GaussJordan = matErr_None
End Function

Function GaussSeidel (A() As Double, B() As Double, X() As Double, numRows As Integer, maxIter As Integer, eps1 As Double, eps2 As Double)

  Const opContinue% = 1
  Const opConverge% = 2
  Const opSingular% = 3
  Const opError% = 4
  
  Dim Xold() As Double
  Dim denom As Double, sum As Double
  Dim dev As Double, devMax As Double
  Dim i As Integer, j As Integer
  Dim iter As Integer
  Dim operType As Integer

  iter = 0
  operType = opContinue
  ReDim Xold(numRows)
  ' normalize matrix A and vector B
  For i = 0 To numRows - 1
    denom = A(i, i)
    If denom < eps1 Then
      GaussSeidel = matErr_Singular
      Exit Function
    End If
    B(i) = B(i) / denom
    For j = 0 To numRows - 1
      A(i, j) = A(i, j) / denom
    Next j
  Next i

  ' perform Gauss-Seidel iteration
  Do While operType = opContinue
    For i = 0 To numRows - 1
      Xold(i) = X(i)
      X(i) = 0
      For j = 0 To numRows - 1
	If j <> i Then
	  X(i) = X(i) - A(i, j) * X(j)
	End If
      Next j
      X(i) = X(i) + B(i)
    Next i

    ' check for the convergence
    devMax = Abs(Xold(0) - X(0)) / X(0)
    For i = 1 To numRows - 1
      dev = Abs(Xold(i) - X(i)) / X(i)
      If dev > devMax Then devMax = dev
    Next i
    If devMax <= eps2 Then
      operType = opConverge
    Else
      iter = iter + 1
      If iter > maxIter Then operType = opError
    End If
  Loop

  Select Case operType
    Case Is = opConverge
      GaussSiedel = matErr_None

    Case Is = opSingular
      GaussSiedel = matErr_Singular

    Case Is = opError
      GaussSiedel = matErr_IterLimit

    Case Else
      GaussSiedel = matErr_None
  End Select
End Function

Sub LUBackSubst (A() As Double, Index() As Integer, numRows As Integer, B() As Double)

  Dim i As Integer, j As Integer
  Dim idx As Integer, k As Integer
  Dim sum As Double

  k = -1
  For i = 0 To numRows - 1
    idx = Index(i)
    sum = B(idx)
    B(idx) = B(i)
    If k > -1 Then
      For j = k To i - 1
	sum = sum - A(i, j) * B(j)
      Next j
    ElseIf sum <> 0 Then
      k = i
    End If
    B(i) = sum
  Next i
  For i = numRows - 1 To 0 Step -1
    sum = B(i)
    For j = i + 1 To numRows - 1
      sum = sum - A(i, j) * B(j)
    Next j
    B(i) = sum / A(i, i)
  Next i
End Sub

Function LUDecomp (A() As Double, Index() As Integer, numRows As Integer, rowSwapFlag As Integer) As Integer

  Dim i As Integer, j As Integer
  Dim k As Integer, iMax As Integer
  Dim large As Double, sum As Double
  Dim z As Double, z2 As Double
  Dim scaleVect() As Double

  ReDim scaleVect(UBound(A, 2))
  ' initialize row interchange flag
  rowSwapFlag = 1
  ' loop to obtain the scaling element
  For i = 0 To numRows - 1
    large = 0
    For j = 0 To numRows - 1
      z2 = Abs(A(i, j))
      If z2 > large Then large = z2
    Next j
    ' no non-zero large value? then exit with an error code
    If large = 0 Then
      LUDecomp = matErr_Singular
      Exit Function
    End If
    scaleVect(i) = 1 / large
  Next i
  For j = 0 To numRows - 1
    For i = 0 To j - 1
      sum = A(i, j)
      For k = 0 To i - 1
	sum = sum - A(i, k) * A(k, j)
      Next k
      A(i, j) = sum
    Next i
    large = 0
    For i = j To numRows - 1
      sum = A(i, j)
      For k = 0 To j - 1
	sum = sum - A(i, k) * A(k, j)
      Next k
      A(i, j) = sum
      z = scaleVect(i) * Abs(sum)
      If z >= large Then
	large = z
	iMax = i
      End If
    Next i
    If j <> iMax Then
      For k = 0 To numRows - 1
	z = A(iMax, k)
	A(iMax, k) = A(j, k)
	A(j, k) = z
      Next k
      rowSwapFlag = -1 * rowSwapFlag
      scaleVect(iMax) = scaleVect(j)
    End If
    Index(j) = iMax
    If A(j, j) = 0 Then
      A(j, j) = MATVECT_EPSILON
    End If
    If j <> numRows Then
      z = 1 / A(j, j)
      For i = j + 1 To numRows - 1
	A(i, j) = A(i, j) * z
      Next i
    End If
  Next j
  LUDecomp = matErr_None
End Function

Function LUDeterminant (A() As Double, numRows As Integer, rowSwapFlag As Integer) As Double

  Dim result As Double
  Dim i As Integer
  
  result = rowSwapFlag
  For i = 0 To numRows - 1
    result = result * A(i, i)
  Next i

  LUDeterminant = result
End Function

Sub LUInverse (A() As Double, InvA() As Double, Index() As Integer, numRows As Integer)

  Dim colVect() As Double
  Dim i As Integer, j As Integer

  ReDim colVect(numRows)
  For j = 0 To numRows - 1
    For i = 0 To numRows - 1
      colVect(i) = 0
    Next i
    colVect(j) = 1
    LUBackSubst A(), Index(), numRows, colVect()
    For i = 0 To numRows - 1
      InvA(i, j) = colVect(i)
    Next i
  Next j
End Sub

Function MatDeterminant (A() As Double, numRows As Integer) As Double

  Dim Index() As Integer
  Dim i As Integer, j As Integer
  Dim rowSwapFlag As Integer
  Dim fxErr As Integer
  Dim result As Double

  ReDim Index(numRows)
  fxErr = LUDecomp(A(), Index(), numRows, rowSwapFlag)
  If Err <> matErr_None Then
    MatDeterminant = MATVECT_BAD_RESULT
    Exit Function
  End If

  result = rowSwapFlag

  For i = 0 To numRows - 1
    result = result * A(i, i)
  Next i
  MatDeterminant = result
End Function

Function MatInverse (A() As Double, numRows As Integer) As Integer

  Dim colVect() As Double
  Dim Index() As Integer
  Dim i As Integer, j As Integer
  Dim rowSwapFlag As Integer
  Dim fxErr As Integer

  ReDim colVect(numRows)
  ReDim Index(numRows)
  Err = LUDecomp(A(), Index(), numRows, rowSwapFlag)
  If Err <> matErr_None Then
    MatInverse = Err
    Exit Function
  End If
    
  For j = 0 To numRows - 1
    For i = 0 To numRows - 1
      colVect(i) = 0
    Next i
    colVect(j) = 1
    LUBackSubst A(), Index(), numRows, colVect()
    For i = 0 To numRows - 1
      A(i, j) = colVect(i)
    Next i
  Next j

  MatInverse = matErr_None
End Function

Function MulMat (MatC() As Double, MatA() As Double, MatB() As Double, numRows As Integer, numCols As Integer) As Integer

  Dim row As Integer, col As Integer, k As Integer

  If Not checkRowCol(MatA(), numRows, numCols) Then
    MulMat = matErr_Size
    Exit Function
  End If

  If Not checkRowCol(MatB(), numRows, numCols) Then
    MulMat = matErr_Size
    Exit Function
  End If

  If numColsA <> numRowsB Then
    MulMat = matErr_Size
    Exit Function
  End If

  For row = 0 To numRowsB - 1
    For col = 0 To numColsA - 1
      MatC(row, col) = 0
      For k = 0 To numColsA - 1
	MatC(row, col) = MatC(row, col) + MatA(row, k) * MatB(row, k)
      Next k
    Next col
  Next row
  MulMat = matErr_None
End Function

Function SubMat (MatC() As Double, MatA() As Double, MatB() As Double, numRows As Integer, numCols As Integer) As Integer

  Dim row As Integer, col As Integer

  If Not (checkRowCol(MatA(), numRows, numCols) And checkRowCol(MatB(), numRows, numCols) And checkRowCol(MatC(), numRows, numCols) And (numCols = numRows)) Then
    SubMat = matErr_Size
    Exit Function
  End If

  For row = 0 To numRows - 1
    For col = 0 To numCols - 1
      MatC(row, col) = MatA(row, col) - MatB(row, col)
    Next col
  Next row
  SubMat = matErr_None
End Function

Sub swapDouble (d1 As Double, d2 As Double)

  Dim t As Double
  
  t = d1
  d1 = d2
  d2 = t
End Sub

Sub swapInt (i1 As Integer, i2 As Integer)
' swap two integers
  Dim t As Integer
  t = i1
  i1 = i2
  i2 = t
End Sub

