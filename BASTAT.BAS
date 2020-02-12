Global Const BASTAT_EPS# = 1E-30

Type basicStat
  hasMissingData As Integer
  countData As Integer
  sum As Double
  sumX As Double
  sumX2 As Double
  sumX3 As Double
  sumX4 As Double
  missingCode As Double
  mean As Double
  sdev As Double
End Type

Function getMeanSdev (x() As Double, numData As Integer, B As basicStat) As Integer

  Dim i As Integer
  Dim xx As Double, x2 As Double

  If numData > 1 Then
     B.countData = B.countData + numData
     If B.hasMissingData Then
        For i = 0 To numData - 1
          xx = x(i)
          If xx > B.missingCode Then
             x2 = xx * xx
             B.sum = B.sum + 1
             B.sumX = B.sumX + xx
             B.sumX2 = B.sumX2 + x2
             B.sumX3 = B.sumX3 + xx * x2
             B.sumX4 = B.sumX4 + x2 * x2
          End If
       Next i
     Else
        For i = 0 To numData - 1
          xx = x(i)
          x2 = xx * xx
          B.sum = B.sum + 1
          B.sumX = B.sumX + xx
          B.sumX2 = B.sumX2 + x2
          B.sumX3 = B.sumX3 + xx * x2
          B.sumX4 = B.sumX4 + x2 * x2
        Next i
     End If
     If B.sum > 1 Then
        B.mean = B.sumX / B.sum
        B.sdev = Sqr((B.sumX2 - B.sumX ^ 2 / B.sum) / (B.sum - 1#))
     Else
        getMeanSdev = False
        Exit Function
     End If
  End If
  getMeanSdev = False
End Function

Sub initializeBasicStat (B As basicStat, hasMissingData As Integer, missingCode As Double)

  B.hasMissingData = hasMissingData
  B.missingCode = missingCode
  B.countData = 0
  B.sum = 0
  B.sumX = 0
  B.sumX2 = 0
  B.sumX3 = 0
  B.sumX4 = 0
  B.mean = 0
  B.sdev = 0
End Sub

Function meanCI (B As basicStat, probability As Double, meanXLow As Double, meanXHi As Double) As Integer
' calculate the confidence interval of the mean value

  Dim tableT As Double, df As Double
  Dim delta As Double, p As Double

  If B.sum < 2 Then
     meanCI = False
     Exit Function
  End If
  df = B.sum - 1#
  If probability > 1# Then
     p = .5 - probability / 200#
  Else
     p = .5 - probability / 2#
  End If
  tableT = TInv(p, df)
  delta = tableT * B.sdev / Sqr(B.sum)
  meanXHi = B.mean + delta
  meanXLow = B.mean - delta
  meanCI = True
End Function

Function meanTest (B1 As basicStat, B2 As basicStat, testedDifference As Double, probability As Double, calcT As Double, tableT As Double, passTest As Integer) As Integer
' test means of two samples

  Dim df As Double, p As Double
  Dim factor1 As Double, factor2 As Double
  Dim n1 As Double, n2 As Double

  If B1.sum < 2 Or B2.sum < 2 Then
     meanTest = False
     Exit Function
  End If

  n1 = B1.sum
  n2 = B2.sum

  If probability > 1# Then
     p = .5 - probability / 200#
  Else
     p = .5 - probability / 2#
  End If
  df = n1 + n2 - 2#
  factor1 = Sqr(1# / n1 + 1# / n2)
  factor2 = Sqr((B1.sumX2 - B1.sumX ^ 2 / n1 + B2.sumX2 - B2.sumX ^ 2 / n2) / df)
  calcT = (B1.sumX / n1 - B2.sumX / n2 - testedDifference) / (factor1 * factor2)
  tableT = TInv(p, df)
  passTest = Abs(calcT) <= tableT
  meanTest = True
End Function

Function moments (B As basicStat, meanX As Double, varianceX As Double, moment3 As Double, moment4 As Double, skewnessCoeff As Double, kurtosisCoeff As Double) As Integer
' Function to calculate the first four moments and
'    the coefficients of skewness and kurtosis.

  Dim meanSqrd As Double

  If B.sum < 2 Then
     moments = False
     Exit Function
  End If

  meanX = B.sumX / B.sum
  meanSqrd = meanX ^ 2
  varianceX = B.sumX2 / B.sum - meanSqrd
  moment3 = B.sumX3 / B.sum - 3# / B.sum * meanX * B.sumX2 + 2# * meanX * meanSqrd
  moment4 = B.sumX4 / B.sum - 4# / B.sum * meanX * B.sumX3 + 6# / B.sum * meanSqrd * B.sumX2 - 3# * meanSqrd ^ 2
  skewnessCoeff = moment3 / varianceX ^ 1.5
  kurtosisCoeff = moment4 / varianceX ^ 2
  moments = True
End Function

Function sdevCI (B As basicStat, probability As Double, sdevXLow As Double, sdevXHi As Double) As Integer
'  calculate the confidence interval of the
'     standard deviation value.

  Dim df As Double, p As Double

  If B.sum < 2 Then
     sdevCI = False
     Exit Function
  End If
  df = B.sum - 1#
  If probability > 1# Then
     p = .5 - probability / 200#
  Else
     p = .5 - probability / 2#
  End If
  sdevXLow = B.sdev * Sqr(df / ChiInv(p, df))
  sdevXHi = B.sdev * Sqr(df / ChiInv((1 - p), df))
  sdevCI = True
End Function

