
' ANOVA table record Typeure
Type LinRegANOVA
  Reg_df As Double
  Reg_SS As Double
  Residual_df As Double
  Residual_SS As Double
  Total_df As Double
  Total_SS As Double
  S2 As Double
  LINREG_F As Double
End Type

' Linear regression record Typeure
Type LinRegrec
  Xindex As Integer
  Yindex As Integer
  hasMissingData As Integer
  missingCode As Double
  sumX As Double
  sumXX As Double
  sum As Double
  sumY As Double
  sumYY As Double
  sumXY As Double
  ' basic statsitics
  MeanX As Double
  MeanY As Double
  SdevX As Double
  SdevY As Double
  Slope As Double
  Intercept As Double
  R2 As Double
  ' Regression results
  StdErSlope As Double
  StdErIntercept As Double
  ANOVA As LinRegANOVA
End Type

Type BestLRrec
  fxIndex As Integer
  R2 As Double
  Slope As Double
  Intercept As Double
End Type

Function BestFit (DataMat() As Double, NData As Integer, numTrnsFx As Integer, Xindex As Integer, Yindex As Integer, r() As BestLRrec) As Integer
' find the best model to fit the data

  Dim sum As Double, sumX As Double
  Dim sumXX As Double, sumY As Double
  Dim sumYY As Double, sumXY As Double
  Dim MeanX As Double, MeanY As Double
  Dim SdevX As Double, SdevY As Double
  Dim x As Double, y As Double
  Dim i As Integer, j As Integer
  Dim offset As Integer, inOrder As Integer
  Dim temp As BestLRrec

  If (Xindex = Yindex) Or (NData < 3) Or (numTrnsFx < 2) Then
     BestFit = False
     Exit Function
  End If

  ' iterate over the functions
  For i = 0 To numTrnsFx - 1
     ' initialize summations
     sum = 0
     sumX = 0
     sumXX = 0
     sumY = 0
     sumYY = 0
     sumXY = 0
     ' process the observations
     For j = 0 To NData - 1
        x = DataMat(j, Xindex)
        y = DataMat(j, Yindex)
        x = MyBestLRFx(x, i)
        y = MyBestLRFy(y, i)
        sum = sum + 1
        sumX = sumX + x
        sumY = sumY + y
        sumXX = sumXX + x * x
        sumYY = sumYY + y * y
        sumXY = sumXY + x * y
     Next j
     ' calculate the results
     MeanX = sumX / sum
     MeanY = sumY / sum
     SdevX = Sqr((sumXX - sumX ^ 2 / sum) / (sum - 1))
     SdevY = Sqr((sumYY - sumY ^ 2 / sum) / (sum - 1))
     r(i).Slope = (sumXY - MeanX * MeanY * sum) / SdevX ^ 2 / (sum - 1)
     r(i).Intercept = MeanY - r(i).Slope * MeanX
     r(i).R2 = (r(i).Slope * SdevX / SdevY) ^ 2
     r(i).fxIndex = i
  Next i
  ' sort the results in descending order
  offset = numTrnsFx
  Do
     offset = (offset * 8) \ 11
     If offset = 0 Then offset = 1
     inOrder = True
     For i = 0 To numTrnsFx - offset - 1
        j = i + offset
        ' swap elements?
        If r(i).R2 < r(j).R2 Then
          inOrder = False
          temp = r(i)
          r(i) = r(j)
          r(j) = temp
        End If
     Next i
  Loop Until (offset = 1) And inOrder
  BestFit = True
End Function

Sub InitializeLinReg (r As LinRegrec, Xindex As Integer, Yindex As Integer, hasMissingData As Integer, missingCode As Double)
' initialize statistical summations and data range

  r.sumX = 0
  r.sumXX = 0
  r.sum = 0
  r.sumY = 0
  r.sumYY = 0
  r.sumXY = 0
  r.hasMissingData = hasMissingData
  r.missingCode = missingCode
  r.Xindex = Xindex
  r.Yindex = Yindex
End Sub

Sub LinReg (r As LinRegrec)
' calculate regression coefficients and related results

  Dim S As Double

  r.MeanX = r.sumX / r.sum
  r.MeanY = r.sumY / r.sum
  r.SdevX = Sqr((r.sumXX - r.sumX ^ 2 / r.sum) / (r.sum - 1))
  r.SdevY = Sqr((r.sumYY - r.sumY ^ 2 / r.sum) / (r.sum - 1))
  r.Slope = (r.sumXY - r.MeanX * r.MeanY * r.sum) / r.SdevX ^ 2 / (r.sum - 1)
  r.Intercept = r.MeanY - r.Slope * r.MeanX
  r.R2 = (r.Slope * r.SdevX / r.SdevY) ^ 2
  r.ANOVA.Reg_SS = (r.sumXY - r.sumY * r.MeanX) ^ 2 / ((r.SdevX) ^ 2 * (r.sum - 1))
  r.ANOVA.Total_SS = r.sumYY - r.sumY ^ 2 / r.sum
  r.ANOVA.Residual_SS = r.ANOVA.Total_SS - r.ANOVA.Reg_SS
  r.ANOVA.Residual_df = r.sum - 2
  S = Sqr(r.ANOVA.Residual_SS / r.ANOVA.Residual_df)
  r.StdErSlope = S / Sqr(r.sumXX - r.sumX ^ 2 / r.sum)
  r.StdErIntercept = S * Sqr(r.sumXX / r.sum / r.SdevX ^ 2 / (r.sum - 1))
  r.ANOVA.Reg_df = 1
  r.ANOVA.Total_df = r.sum - 1
  r.ANOVA.S2 = S * S
  r.ANOVA.LINREG_F = r.ANOVA.Reg_SS / r.ANOVA.S2
End Sub

Sub LinRegCoefCI (r As LinRegrec, probability As Double, slopeHi As Double, slopeLow As Double, intHi As Double, intLow As Double)
' calculate confidence interval for slope and intercept

  Dim Df As Double, tableT As Double
  Dim diff As Double, p As Double

  If probability > 1 Then
     p = .5 - probability / 200
  Else
     p = .5 - probability / 2
  End If
  Df = r.sum - 2
  tableT = TInv(p, Df)
  diff = tableT * r.StdErSlope
  slopeHi = r.Slope + diff
  slopeLow = r.Slope - diff

  diff = tableT * r.StdErIntercept
  intHi = r.Intercept + diff
  intLow = r.Intercept - diff
End Sub

Sub LR_Int_T_test (r As LinRegrec, probability As Double, testValue As Double, calcT As Double, tableT As Double, passTest As Integer)
' compare intercept value with a tested value
'   Hypothesis tested is H0: Intercept = testValue

  Dim Df As Double, p As Double

  If probability > 1 Then
     p = .5 - probability / 200
  Else
     p = .5 - probability / 2
  End If
  Df = r.sum - 2
  tableT = TInv(Df, p)
  calcT = (r.Intercept - testValue) / r.StdErIntercept
  passTest = Abs(calcT) <= tableT
End Sub

Sub LR_R2_T_Test (r As LinRegrec, probability As Double, calcT As Double, tableT As Double, passTest As Integer)
' Procedure to test hypothesis H0 : R^2 = 0

  Dim Df As Double, p As Double

  If probability > 1 Then
     p = .5 - probability / 200
  Else
     p = .5 - probability / 2
  End If
  Df = r.sum - 2
  tableT = TInv(p, Df)
  calcT = Sqr(r.R2 * Df / (1 - r.R2))
  passTest = calcT <= tableT
End Sub

Sub LR_Slope_T_test (r As LinRegrec, probability As Double, testValue As Double, calcT As Double, tableT As Double, passTest As Integer)
' compare slope value with a tested value
'   Hypothesis tested is H0: Slope = testValue

  Dim Df As Double, p As Double

  If probability > 1 Then
     p = .5 - probability / 200
  Else
     p = .5 - probability / 2
  End If
  Df = r.sum - 2
  tableT = TInv(p, Df)
  calcT = (r.Slope - testValue) / r.StdErSlope
  passTest = Abs(calcT) <= tableT
End Sub

Sub sumLinReg (DataMat() As Double, r As LinRegrec, NData As Integer)
' update statistical summations

  Dim i As Integer
  Dim xr As Double, yr As Double

  For i = 0 To NData - 1
     xr = DataMat(i, r.Xindex)
     yr = DataMat(i, r.Yindex)
     If Not (r.hasMissingData And ((xr <= r.missingCode) Or (yr <= r.missingCode))) Then
        ' transform x and y data
        xr = MyLRFx(xr)
        yr = MyLRFy(yr)
        ' Update summations
        r.sum = r.sum + 1
        r.sumXY = r.sumXY + xr * yr
        r.sumX = r.sumX + xr
        r.sumXX = r.sumXX + xr ^ 2
        r.sumY = r.sumY + yr
        r.sumYY = r.sumYY + yr ^ 2
     End If
  Next i
End Sub

Sub YhatCI (r As LinRegrec, xr As Double, probability As Double, yHat As Double, yHi As Double, yLow As Double)
' calculate projections and their confidence interval

  Dim Df As Double, deltaY As Double
  Dim p As Double, tableT As Double

  If probability > 1 Then
     p = .5 - probability / 200
  Else
     p = .5 - probability / 2
  End If
  Df = r.sum - 2
  tableT = TInv(p, Df)
  xr = MyLRFx(xr) ' transform xr
  deltaY = Sqr((xr - r.MeanX) ^ 2 / (r.sumX ^ 2 * (r.sum - 1)) + 1 / (r.sum + 1)) * Sqr(r.ANOVA.S2) * tableT
  yHat = r.Intercept + r.Slope * xr
  yHi = MyLRInvFy(yHat + deltaY)
  yLow = MyLRInvFy(yHat - deltaY)
  yHat = MyLRInvFy(yHat)
End Sub

