      PROGRAM Test1
	  
      REAL, PARAMETER:: MYMAT(3,3)=reshape([
     1  11.0, 12.0, 13.0,
     2  21.0, 22.0, 23.0,
     3  31.0, 32.0, 33.0], [3,3])

      print *, 'MYMAT = '
      print *, MYMAT
      print *, '---------------------- '
	  
      print *, 'MYMAT(2,1)= ', MYMAT(2,1)
	  
	  
      END PROGRAM