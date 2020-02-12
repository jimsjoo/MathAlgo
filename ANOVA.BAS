Global Const ANOVA_EPS# = 1E-30

Type ANOVA1rec
  NVar As Integer
  Treatment_SS As Double
  Treatment_df As Double
  Treatment_MS As Double
  Error_SS As Double
  Error_df As Double
  Error_MS As Double
  Total_SS As Double
  Total_df As Double
  ANOVA1_F As Double
  hasMissingData As Integer
  missingCode As Double
End Type

Type ANOVA2rec
  Num_Row As Integer
  Num_Col As Integer
  GrandSumX As Double
  GrandSumX2 As Double
  Row_SS As Double
  Row_df As Double
  Row_F As Double
  Col_SS As Double
  Col_df As Double
  Col_F As Double
  Error_SS As Double
  Error_df As Double
  Total_SS As Double
  Total_df As Double
End Type

Type ANOVA2Rrec
  NReplicate As Integer
  Num_Row As Integer
  Num_Col As Integer
  GrandSumX As Double
  GrandSumX2 As Double
  Row_SS As Double
  Row_df As Double
  Row_F As Double
  Col_SS As Double
  Col_df As Double
  Col_F As Double
  Interact_SS As Double
  Interact_df As Double
  Interact_F As Double
  Error_SS As Double
  Error_df As Double
  Total_SS As Double
  Total_df As Double
End Type

Type LatinSqrrec
  N As Integer
  Num As Double
  GrandSumX As Double
  GrandSumX2 As Double
  Row_SS As Double
  Row_df As Double
  Row_F As Double
  Col_SS As Double
  Col_df As Double
  Col_F As Double
  Trt_SS As Double
  Trt_df As Double
  Trt_F As Double
  Error_SS As Double
  Error_df As Double
  Total_SS As Double
  Total_df As Double
End Type

Type ANOCOVrec
  Num_Set As Integer
  hasMissingData As Integer
  missingCode As Double
  GrandSum As Double
  Df1 As Double
  Df2 As Double
  Df3 As Double
  Df4 As Double
  ASSx As Double
  ASP As Double
  ASSy As Double
  ASSyhat As Double
  AMSyhat As Double
  WSSx As Double
  WSP As Double
  WSSy As Double
  WSSyhat As Double
  WMSyhat As Double
  TSSx As Double
  TSP As Double
  TSSy As Double
  TSSyhat As Double
  ANOCOV_F As Double
  AMSx As Double
  WMSx As Double
  AMSy As Double
  WMSy As Double
  Fx As Double
  Fy As Double
End Type





Sub ANOCOV (r As ANOCOVrec, DataMat() As Double, NData As Integer, NumSets As Integer, hasMissingData As Integer, missingCode As Double, sum() As Double, sumX() As Double, sumY() As Double, sumX2() As Double, sumY2() As Double, sumXY() As Double)
  Dim i As Integer, j As Integer
  Dim jX As Integer, jY As Integer
  Dim x As Double, y As Double
  Dim dSumX As Double, dSumY As Double
  Dim dSumX2 As Double, dSumY2 As Double
  Dim dSumXY As Double
  Dim N As Integer

  N = 2 * NumSets

  ' initilaize statistical summations
  For i = 0 To N - 1
     sum(i) = 0
     sumX(i) = 0
     sumX2(i) = 0
     sumY(i) = 0
     sumY2(i) = 0
     sumXY(i) = 0
  Next i
  r.GrandSum = 0
  r.ASSx = 0
  r.ASSy = 0
  r.ASP = 0
  r.hasMissingData = hasMissingData
  r.missingCode = missingCode

  ' update statistical summations
  r.Num_Set = NumSets
  For i = 0 To NData - 1
     For j = 0 To r.Num_Set - 1
        jY = 2 * j + 1
        jX = jY - 1
        x = DataMat(i, jX)
        If Not (r.hasMissingData And x <= r.missingCode) Then
          y = DataMat(i, jY)
          sum(j) = sum(j) + 1
          sumX(j) = sumX(j) + x
          sumX2(j) = sumX2(j) + x ^ 2
          sumY(j) = sumY(j) + y
          sumY2(j) = sumY2(j) + y ^ 2
          sumXY(j) = sumXY(j) + x * y
        End If
     Next j
  Next i
  ' carry out ANOCOV calculations
  dSumX = 0
  dSumY = 0
  dSumX2 = 0
  dSumY2 = 0
  dSumXY = 0
  For j = 0 To r.Num_Set - 1
     r.ASSx = r.ASSx + sumX(j) ^ 2 / sum(j)
     r.ASSy = r.ASSy + sumY(j) ^ 2 / sum(j)
     r.ASP = r.ASP + sumX(j) * sumY(j) / sum(j)
     r.GrandSum = r.GrandSum + sum(j)
     dSumX = dSumX + sumX(j)
     dSumX2 = dSumX2 + sumX2(j)
     dSumY = dSumY + sumY(j)
     dSumY2 = dSumY2 + sumY2(j)
     dSumXY = dSumXY + sumXY(j)
  Next j
  r.TSSx = dSumX2 - dSumX ^ 2 / r.GrandSum
  r.TSSy = dSumY2 - dSumY ^ 2 / r.GrandSum
  r.ASSx = r.ASSx - dSumX ^ 2 / r.GrandSum
  r.ASSy = r.ASSy - dSumY ^ 2 / r.GrandSum
  r.WSSx = r.TSSx - r.ASSx
  r.WSSy = r.TSSy - r.ASSy
  r.Df1 = r.Num_Set - 1
  r.Df3 = r.Df1
  r.Df2 = r.GrandSum - r.Num_Set
  r.Df4 = r.Df2 - 1
  r.TSP = dSumXY - dSumX * dSumY / r.GrandSum
  r.ASP = r.ASP - dSumX * dSumY / r.GrandSum
  r.WSP = r.TSP - r.ASP
  r.AMSx = r.ASSx / r.Df1
  r.AMSy = r.ASSy / r.Df1
  r.WMSx = r.WSSx / r.Df2
  r.WMSy = r.WSSy / r.Df2
  r.Fx = r.AMSx / r.WMSx
  r.Fy = r.AMSy / r.WMSy
  r.TSSyhat = r.TSSy - r.TSP ^ 2 / r.TSSx
  r.WSSyhat = r.WSSy - r.WSP ^ 2 / r.WSSx
  r.ASSyhat = r.TSSyhat - r.WSSyhat
  r.AMSyhat = r.ASSyhat / r.Df3
  r.WMSyhat = r.WSSyhat / r.Df4
  r.ANOCOV_F = r.AMSyhat / r.WMSyhat
End Sub

Sub ANOVA1 (r As ANOVA1rec, DataMat() As Double, NData As Integer, NVar As Integer, hasMissingData As Integer, missingCode As Double, sum() As Double, sumX() As Double, sumX2() As Double)
  Dim i As Integer, j As Integer
  Dim x As Double
  Dim TS As Double, TSS As Double
  Dim GrandTotal As Double
  Dim Mean() As Double
  Dim Sdev() As Double

  ReDim Mean(NVar)
  ReDim Sdev(NVar)

  For j = 0 To NVar - 1
     sum(j) = 0
     sumX(j) = 0
     sumX2(j) = 0
  Next j
  r.hasMissingData = hasMissingData
  r.missingCode = missingCode

  ' update statistical summations
  r.NVar = NVar
  For i = 0 To NData - 1
     For j = 0 To NVar - 1
        x = DataMat(i, j)
        If Not (r.hasMissingData And x <= r.missingCode) Then
          sum(j) = sum(j) + 1
          sumX(j) = sumX(j) + x
          sumX2(j) = sumX2(j) + x * x
        End If
     Next j
  Next i

  ' perform one way ANOVA calculations
  TS = 0
  TSS = 0
  GrandTotal = 0
  For i = 0 To r.NVar - 1
     Mean(i) = sumX(i) / sum(i)
     Sdev(i) = Sqr((sumX2(i) - sumX(i) ^ 2 / sum(i)) / (sum(i) - 1))
     TS = TS + sumX(i)
     TSS = TSS + sumX2(i)
     GrandTotal = GrandTotal + sum(i)
  Next i
  r.Total_SS = TSS - TS ^ 2 / GrandTotal
  r.Treatment_SS = 0
  For i = 0 To r.NVar - 1
     r.Treatment_SS = r.Treatment_SS + sumX(i) ^ 2 / sum(i)
  Next i
  r.Treatment_SS = r.Treatment_SS - TS ^ 2 / GrandTotal
  r.Error_SS = r.Total_SS - r.Treatment_SS
  r.Treatment_df = r.NVar - 1
  r.Error_df = GrandTotal - r.Treatment_df - 1
  r.Total_df = r.Error_df + r.Treatment_df
  r.Treatment_MS = r.Treatment_SS / r.Treatment_df
  r.Error_MS = r.Error_SS / r.Error_df
  r.ANOVA1_F = r.Treatment_MS / r.Error_MS
End Sub

Sub ANOVA2 (r As ANOVA2rec, DataMat() As Double, NData As Integer, NVar As Integer, ColSumX() As Double, ColSumX2() As Double, RowSumX() As Double, RowSumX2() As Double)
  Dim i As Integer, j As Integer
  Dim nCol As Integer, nRow As Integer
  Dim RowCol As Double
  Dim x As Double, sx As Double

  ' initilaize statistical summations
  For j = 0 To NVar - 1
     ColSumX(j) = 0
     ColSumX2(j) = 0
  Next j

  For j = 0 To NData - 1
     RowSumX(j) = 0
     RowSumX2(j) = 0
  Next j
  r.Num_Row = NData
  r.Num_Col = NVar

  ' update statistical summations
  r.GrandSumX = 0
  r.GrandSumX2 = 0
  r.Num_Col = NVar
  r.Num_Row = NData
  For i = 0 To NData - 1
     For j = 0 To NVar - 1
        x = DataMat(i, j)
        sx = x * x
        RowSumX(i) = RowSumX(i) + x
        ColSumX(j) = ColSumX(j) + x
        RowSumX2(i) = RowSumX2(i) + sx
        ColSumX2(j) = ColSumX2(j) + sx
     Next j
  Next i

  For i = 0 To NData - 1
     r.GrandSumX = r.GrandSumX + RowSumX(i)
     r.GrandSumX2 = r.GrandSumX2 + RowSumX2(i)
  Next i
  ' perform two-way ANOVA calculations
  nCol = r.Num_Col
  nRow = r.Num_Row
  RowCol = r.Num_Row * r.Num_Col
  r.Total_SS = r.GrandSumX2 - r.GrandSumX ^ 2 / RowCol
  r.Row_SS = 0
  For i = 0 To nRow - 1
     r.Row_SS = r.Row_SS + RowSumX(i) ^ 2 / nCol
  Next i
  r.Row_SS = r.Row_SS - r.GrandSumX ^ 2 / RowCol
  r.Col_SS = 0
  For i = 0 To nCol - 1
     r.Col_SS = r.Col_SS + ColSumX(i) ^ 2 / nRow
  Next i
  r.Col_SS = r.Col_SS - r.GrandSumX ^ 2 / RowCol
  r.Error_SS = r.Total_SS - r.Row_SS - r.Col_SS
  r.Row_df = r.Num_Row - 1
  r.Col_df = r.Num_Col - 1
  r.Error_df = r.Row_df * r.Col_df
  r.Total_df = RowCol - 1
  r.Row_F = (r.Row_SS / r.Row_df) / (r.Error_SS / r.Error_df)
  r.Col_F = (r.Col_SS / r.Col_df) / (r.Error_SS / r.Error_df)
End Sub

Sub ANOVA2R (r As ANOVA2Rrec, DataMat() As Double, NData As Integer, NVar As Integer, NumReplicates As Integer, ColSumX() As Double, RowSumX() As Double)
  Dim x As Double, sx As Double
  Dim nCol As Integer, nRow As Integer
  Dim RowCol As Double, sum As Double
  Dim i As Integer, j As Integer
  Dim k As Integer, iRow As Integer
  Dim CellSum() As Double

  ' initiliaze statistical summations
  ReDim CellSum(0 To NData, 0 To NVar)
  For j = 0 To NVar - 1
     ColSumX(j) = 0
  Next j
  For j = 0 To NData - 1
     RowSumX(j) = 0
  Next j
  For i = 0 To NData - 1
     For j = 0 To NVar - 1
        CellSum(i, j) = 0
     Next j
  Next i
  r.Num_Row = NData
  r.Num_Col = NVar
  r.GrandSumX = 0
  r.GrandSumX2 = 0

  ' update statistical summations
  r.NReplicate = NumReplicates
  For i = 0 To NData - 1
     iRow = i \ r.NReplicate
     For j = 0 To NVar - 1
        x = DataMat(i, j)
        sx = x * x
        CellSum(iRow, j) = CellSum(iRow, j) + x
        RowSumX(iRow) = RowSumX(iRow) + x
        ColSumX(j) = ColSumX(j) + x
        r.GrandSumX2 = r.GrandSumX2 + sx
     Next j
  Next i
  k = NData / r.NReplicate
  r.Num_Row = k
  r.Num_Col = NVar
  For j = 0 To NVar - 1
     r.GrandSumX = r.GrandSumX + ColSumX(j)
  Next j
  ' perform two-way ANOVA calculations
  nCol = r.Num_Col
  nRow = r.Num_Row
  RowCol = r.Num_Row * r.Num_Col * r.NReplicate
  r.Total_SS = r.GrandSumX2 - r.GrandSumX ^ 2 / RowCol
  r.Row_SS = 0
  For i = 0 To nRow - 1
     r.Row_SS = r.Row_SS + RowSumX(i) ^ 2 / (r.Num_Col * r.NReplicate)
  Next i
  r.Row_SS = r.Row_SS - r.GrandSumX ^ 2 / RowCol
  r.Col_SS = 0
  For i = 0 To nCol - 1
     r.Col_SS = r.Col_SS + ColSumX(i) ^ 2 / (r.Num_Row * r.NReplicate)
  Next i
  r.Col_SS = r.Col_SS - r.GrandSumX ^ 2 / RowCol
  sum = 0
  For i = 0 To nRow - 1
     For j = 0 To nCol - 1
        sum = sum + CellSum(i, j) ^ 2
     Next j
  Next i
  r.Interact_SS = sum / r.NReplicate - r.Row_SS - r.Col_SS - r.GrandSumX ^ 2 / RowCol
  r.Error_SS = r.GrandSumX2 - sum / r.NReplicate
  r.Total_SS = r.GrandSumX2 - r.GrandSumX ^ 2 / RowCol
  r.Row_df = r.Num_Row - 1
  r.Col_df = r.Num_Col - 1
  r.Interact_df = r.Row_df * r.Col_df
  r.Total_df = RowCol - 1
  r.Error_df = r.Num_Row * r.Num_Col * (r.NReplicate - 1)
  r.Row_F = (r.Row_SS / r.Error_SS) * (r.Error_df / r.Row_df)
  r.Col_F = (r.Col_SS / r.Error_SS) * (r.Error_df / r.Col_df)
  r.Interact_F = (r.Interact_SS / r.Error_SS) * (r.Error_df / r.Interact_df)
End Sub

Sub Latin (r As LatinSqrrec, DataMat() As Double, Map() As Integer, NData As Integer, ColSumX() As Double, RowSumX() As Double, TrtSumX() As Double)
  Dim i As Integer, j As Integer, k As Integer
  Dim x As Double, sx As Double, g As Double

  ' initilaize statistical summations
  For i = 0 To NData - 1
     ColSumX(i) = 0
     RowSumX(i) = 0
     TrtSumX(i) = 0
  Next i
  r.GrandSumX = 0
  r.GrandSumX2 = 0
  r.Num = 0

  ' update statistical summations
  r.N = NData
  ' Obtain sums of colums, rows and treatments
  For i = 0 To r.N - 1
     r.Num = r.Num + 1
     For j = 0 To r.N - 1
        x = DataMat(i, j)
        sx = x * x
        ColSumX(j) = ColSumX(j) + x
        RowSumX(i) = RowSumX(i) + x
        k = Map(i + NData * j)
        TrtSumX(k) = TrtSumX(k) + x
        r.GrandSumX2 = r.GrandSumX2 + sx
     Next j
     r.GrandSumX = r.GrandSumX + RowSumX(i)
  Next i

  ' perform Latin Squares ANOVA
  r.Row_SS = 0
  r.Col_SS = 0
  r.Trt_SS = 0
  For i = 0 To r.N - 1
     r.Row_SS = r.Row_SS + RowSumX(i) ^ 2
     r.Col_SS = r.Col_SS + ColSumX(i) ^ 2
     r.Trt_SS = r.Trt_SS + TrtSumX(i) ^ 2
  Next i
  g = (r.GrandSumX / r.Num) ^ 2
  r.Row_SS = r.Row_SS / r.Num - g
  r.Col_SS = r.Col_SS / r.Num - g
  r.Trt_SS = r.Trt_SS / r.Num - g
  r.Error_SS = r.GrandSumX2 - r.Row_SS - r.Col_SS - r.Trt_SS - g
  r.Total_SS = r.GrandSumX2 - g
  r.Row_df = r.Num - 1
  r.Col_df = r.Num - 1
  r.Trt_df = r.Num - 1
  r.Error_df = (r.Num - 1) * (r.Num - 2)
  r.Total_df = r.Num ^ 2 - 1
  r.Row_F = (r.Num - 2) * r.Row_SS / r.Error_SS
  r.Col_F = (r.Num - 2) * r.Col_SS / r.Error_SS
  r.Trt_F = (r.Num - 2) * r.Trt_SS / r.Error_SS
End Sub

