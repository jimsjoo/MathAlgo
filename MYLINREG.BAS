

Function MyBestLRFx (X As Double, Index As Integer) As Double

   Select Case Index
   
     Case 0, 1
       MyBestLRFx = X
       
     Case 2, 3
       MyBestLRFx = Log(X)

     Case Else
       MyBestLRFx = X
    End Select
End Function

Function MyBestLRFy (Y As Double, Index As Integer) As Double

   Select Case Index
   
     Case 0, 2
       MyBestLRFy = Y
       
     Case 1, 3
       MyBestLRFy = Log(Y)

     Case Else
       MyBestLRFy = Y
       
    End Select
End Function

Function MyLRFx (X As Double) As Double
  MyLRFx = X
End Function

Function MyLRFy (Y As Double) As Double
  MyLRFy = Y
End Function

Function MyLRInvFy (Y As Double) As Double
  MyLRInvFy = Y
End Function

