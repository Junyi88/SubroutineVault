      PROGRAM Test2
	  
      real*8 ::  DSPIN(3), TRM0(3,3) , TERM(3,3)
      real*8 :: DROT(3,3) 
      DROT=reshape([
     1  0.903290522317386, -0.259014389274174, 0.342020143325669,
     2  0.414364123486001, 0.733346409125113, -0.538985544695756,
     3  -0.111214232269356, 0.628581411093448, 0.769751131320057], [3,3]
     4 )

      print *, 'DROT = '
      print *, DROT
      print *, '---------------------- '
      call DealWithRotation(DSPIN,DROT,TRM0,TERM)
      print *, 'DSPIN = '
      print *, DSPIN
      print *, '---------------------- '
      print *, 'TRM0 = '
      print *, TRM0
      print *, '---------------------- '
      print *, 'TRM0 (1,2) = '
      print *, TRM0(1,2)
	  
	  print *, '---------------------- '
      print *, 'TERM  = '
      print *, TERM
      END PROGRAM
	  
      INCLUDE 'DealWithRotation.f'   