SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
      
!
!
        implicit none
	!include 'ABA_PARAM.INC' 
	include 'params_xtal.inc' 


        integer lop, lrestart, kstep, kinc
!
	real*8 TIME(2), dtime

       !write(*,*) 'uexternaldb is called,', kstep, kinc, lop

	! Only call the externaldb at the start of the analysis
        if (lop==0 .or. lop==4) then 
!----------------------------------------------------------------------
        write(*,*) 'Start UEXTERNALDB'

	call SetUpCrystalProps( )
	write(*,*) 'uexternaldb done'
      endif

!----------------------------------------------------------------------
	RETURN


END



