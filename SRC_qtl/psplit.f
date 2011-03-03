      SUBROUTINE psplit(nemin,nemax,icase,iso)
      USE param
      USE case       
      USE abc          
      IMPLICIT REAL*8 (A-H,O-Z)
      
      CHARACTER *4     MODUS

      COMMON /CHAR/   MODUS

        nn=iso*(2*lcase(icase)+1)
        itap=29
           dummy=0.d0
        if (MODUS.eq.'SUMA') then
        DO 872 NUM=NEMIN,NEMAX
        WRITE(ITAP,204) NUM,E(NUM),dummy
	i0=0
        DO 100 II=1,isum(0,icase)
	sum=0d0
	do i=i0+1,isum(ii,icase)
	sum=sum+xden(i,num)
	end do
	i0=isum(ii,icase)
            WRITE(ITAP,31)sum
 100    continue
 872    CONTINUE

	elseif (MODUS.eq.'FULL') then
	DO 873 NUM=NEMIN,NEMAX
        WRITE(ITAP,204) NUM,E(NUM),dummy  
	DO II=1,NN
            WRITE(ITAP,31)xden(ii,num)
        END DO
 873    CONTINUE       

	elseif (MODUS.eq.'TOTA') then
	DO 874 NUM=NEMIN,NEMAX
        WRITE(ITAP,204) NUM,E(NUM),dummy
            WRITE(ITAP,31)1.d0
 874    CONTINUE

	elseif (MODUS.eq.'SPIN') then
	read(45,*)(xden(1,num),num=1,nemax)
        read(46,*)(xden(2,num),num=1,nemax)
        DO 875 NUM=NEMIN,NEMAX
        WRITE(ITAP,204) NUM,E(NUM),dummy
            WRITE(ITAP,31)xden(1,num)
            WRITE(ITAP,31)xden(2,num)
 875    CONTINUE

        else
        STOP 'MODUS UNKNOWN'

        end if

      RETURN
                            
   31 FORMAT(1X,'*****:',e17.8)                               
  204 FORMAT(1X,' BAND #',I3,'  E=',F9.5,'  WEIGHT=',F10.7)             
                                  
      END                                                               
