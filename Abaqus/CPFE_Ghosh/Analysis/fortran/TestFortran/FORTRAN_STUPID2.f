      PROGRAM FORTRAN_STUPID2
	  
      REAL*8, PARAMETER:: FCC_OHM(3,12)=reshape([
     1  1.1,2.1,3.1,
     2  1.2,2.2,3.2,
     3  1.3,2.3,3.3,
     4  1.4,2.4,3.4,
     5  1.5,2.5,3.5,
     6  1.6,2.6,3.6, 
     7  1.7,2.7,3.7,	 
     8  1.8,2.8,3.8,	 
     9  1.9,2.9,3.9, 
     1  1.10,2.10,3.10, 
     1  1.11,2.11,3.11, 
     2  1.12,2.12,3.12 
     1  ], [3,12])


      DO IX=1,3
      DO ISLIPS=1,12
	     print *, '(', IX,ISLIPS,')=',FCC_OHM(IX,ISLIPS)
      END DO
      END DO


	  
      END PROGRAM
	  
