��, (  P    Form1� 
Derivation8  �  �  �   B #����$Form158  6�  7�  8�  �   TestMnu�Test��0   	FBDiffMnu�Forward/Backward Difference �(   
CentralMnu�Central Difference �4   ExtCentralMnu�Extended Central Difference �   n1�-��   ExitMnu�&Exit�	  �  __	 
�р�-   � � � �� �� � � �Z �7�Qx �g �� U  �K�� {                  -  �  	Form_Load� ExitMnu_Click� n1_Click GaussJordanMnu_ClickB row�  numRows� col= MAT   A� B printf�  GaussJordan� numColso matErr_None� forrowf forcol myErr� Statc� GaussSeidelMnu_Clickk B2�xp GaussSeidel� LUDecompMnu_Click� dFlag� LUDecompC Indexc LUBackSubst� BisectionMnu_ClickO low� high� params [0]C 	Bisection  	tolerance? rootD f1+ NewtonsApproxMnu_Click n� NewtonApproxi maxIterR NewtonExactMnu_Clicko NewtonExact� NewtonApproxMnu_Click� RichmondApproxMnu_Click� RichmondApprox� RichmondExactMnu_Click� RichmondExact� CombinedMnu_Click3 Combined� BrentMnu_ClickL Brent� MultiRootaMnu_Click� rootsU NewtonMultiRoots numRoots-�iy putchar� PolyDefMnu_Click� 	polyOrder, coeff= [1]� [2]� [3]D fabs� SMALL� DeflatePolyRoots� LinBairstowMnu_Click LBPolyRoots� 	polyRoots� 	isComplex   real   imagu puts� polyRoot8 TwoRootsMnu_Click� paramsX� paramsY~ rootX� rootY   Newton2Functions   fx� fYZ SneMnu_ClickF numEqns4 VEC{ Xguess� NewtonSimNLE ff lg5 ff1   ff2   newVectS MultiRootsMnu_Click� Command1_Clickp Form& LagrangianMnu_Click   xarrG []��yArr� table� derivZ h   sizeofh yf LagrangeQ BareycentricMnu_Click� Barycentric   Statib doubel� EDBarycentricMnu_Click� ED_Barycentric� NewtonDivDiffMnu_Click NewtonDivDiff   NewtonDiffMnu_Click 
NewtonDiff   CubicSplineMnu_Clickr getDerivatives   Spline� dummy��
ARRAY_SIZE  �deltaX  �y0Idx   MyFx� FBDiffMnu_Click� FBDiffDeriv1� FBDiffDeriv2� FBDiffDeriv3 FBDiffDeriv45 FBUdfDeriv1& FBUdfDeriv2   FBUdfDeriv3   FBUdfDeriv4   CentralMnu_Click� CDiffDeriv1   CDiffDeriv2   CDiffDeriv3   CDiffDeriv4   
CUdfDeriv1   
CUdfDeriv2   
CUdfDeriv3   
CUdfDeriv4   ExtCentralMnu_Click� XCDiffDeriv1  XCDiffDeriv2 XCDiffDeriv3 XCDiffDeriv4   XCUdfDeriv1   XCUdfDeriv2   XCUdfDeriv3   XCUdfDeriv4l       �  Z   � *   Z     ,  � :     =  � J     h  � d     ?    r    	  ��������    CentralMnu_Click0J     X      w |  w f � & "   Using Central difference formula  u  w f r  � = � r w f �  l  Using diffence values u  w f �  �  y'( t  =t � 
 �  ) =   t    , h Z r  &� 
 �  ##.## u  w f � 
 �  y''(  t  =t � 
  ) =   t    , h Z r  5� 
 4 ##.## u  w f � 
 L y'''( t  =t � 
 b ) =   t    , h Z r  D� 
 � ##.## u  w f �  � y''''(  t  =t � 
 � ) =   t    , h Z r  S� 
 � ##.## u  w f r  w f �    Using user-defined equation u  w f �  > y(' t  =t � 
 R ) =   t  = h  b� 
 t ##.## u  w f � 
 � y(''  t  =t � 
 � ) =   t  = h  p� 
 � ##.## u  w f � 
 � y(''' t  =t � 
 � ) =   t  = h  ~� 
  ##.## u  w f �  , y''''(  t  =t � 
 D ) =   t  = h  �� 
 f ##.## u   9 	  ��������     ExitMnu_Click0      X  g     7   9 	  ��������     ExtCentralMnu_Click0R     X  �    w |  w f � . " ) Using Extended Central difference formula u  w f r  � = � r w f �  t  Using diffence values u  w f �  �  y'( t  =t � 
 �  ) =   t    , h Z r  �� 
 �  ##.## u  w f � 
 �  y''(  t  =t � 
  ) =   t    , h Z r  �� 
 < ##.## u  w f � 
 T y'''( t  =t � 
 j ) =   t    , h Z r  �� 
 � ##.## u  w f �  � y''''(  t  =t � 
 � ) =   t    , h Z r  �� 
 � ##.## u  w f r  w f �    Using user-defined equation u  w f �  F y(' t  =t � 
 Z ) =   t  = h  �� 
 | ##.## u  w f � 
 � y(''  t  =t � 
 � ) =   t  = h   � 
 � ##.## u  w f � 
 � y(''' t  =t � 
 � ) =   t  = h  � 
  ##.## u  w f �  4 y''''(  t  =t � 
 L ) =   t  = h  � 
 n ##.## u   9 	  ��������     FBDiffMnu_Click0R     X  �    w |  w f � . " ) Using Forward/Backward difference formula u  w f r  � = � r w f �  t  Using diffence values u  w f �  �  y'( t  =t � 
 �  ) =   t    , h Z r  �� 
 �  ##.## u  w f � 
 �  y''(  t  =t � 
  ) =   t    , h Z r  �� 
 < ##.## u  w f � 
 T y'''( t  =t � 
 j ) =   t    , h Z r  �� 
 � ##.## u  w f �  � y''''(  t  =t � 
 � ) =   t    , h Z r  �� 
 � ##.## u  w f r  w f �    Using user-defined equation u  w f �  F y(' t  =t � 
 Z ) =   t  = h  �� 
 | ##.## u  w f � 
 � y(''  t  =t � 
 � ) =   t  = h  �� 
 � ##.## u  w f � 
 � y(''' t  =t � 
 � ) =   t  = h  �� 
  ##.## u  w f �  4 y''''(  t  =t � 
 L ) =   t  = h  � 
 n ##.## u   9 	  ��������    	 Form_Load0�      X  Z     �       �? = � �������? h  ?�  Z�� > ��x   =  { ?  ,  = h�  =  ?N ��   9 	  ��������   �