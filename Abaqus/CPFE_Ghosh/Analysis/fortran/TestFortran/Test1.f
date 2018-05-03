      PROGRAM Test1
	  
      REAL, PARAMETER:: MYMAT(3,3)=reshape([
     1  11.0, 12.0, 13.0,
     2  21.0, 22.0, 23.0,
     3  31.0, 32.0, 33.0], [3,3])

	  REAL, PARAMETER:: MYMAT2(2,2,2)=reshape([
     1  1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0], [2,2,2])
	 
      print *, '---------------------- '
      print *, 'MYMAT = '
      print *, MYMAT
      
      print *, 'MYMAT(2,1)= ', MYMAT(2,1)
	  
      print *, '---------------------- '
      print *, 'MYMAT2 = '
      print *, MYMAT2
      
      print *, 'MYMAT(1,1,1)= ', MYMAT2(1,1,1)
      print *, 'MYMAT(1,1,2)= ', MYMAT2(1,1,2)
      print *, 'MYMAT(1,2,1)= ', MYMAT2(1,2,1)
      print *, 'MYMAT(1,2,2)= ', MYMAT2(1,2,2)
      print *, 'MYMAT(2,1,1)= ', MYMAT2(2,1,1)
      print *, 'MYMAT(2,1,2)= ', MYMAT2(2,1,2)
      print *, 'MYMAT(2,2,1)= ', MYMAT2(2,2,1)
      print *, 'MYMAT(2,2,2)= ', MYMAT2(2,2,2)
	  
      END PROGRAM