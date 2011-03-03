!BOP
! !ROUTINE: LatGen
! !INTERFACE:
      SUBROUTINE LatGen
! !USES:
      use constants, only : pi
	  use rotation_matrices, only : br1
      use struct,only : alat,alpha,nat,lattic
	  use program_control, only : verbosity
! !DESCRIPTION:
!     This routine sets up the Bravais matrix.
!                                                                       
!     LATGEN GENERATES A BRAVAIS MATRIX AND DEFINES THE VOLUME OF THE UNIT CELL                                    
!     BR1(3,3)  : TRANSFORMS INTEGER RECIPROCAL LATTICE VECTORS AS      
!                 GIVEN IN THE VECTORLIST OF LAPW1  ( GENERATED IN      
!                 COORS, TRANSFORMED IN BASISO, AND WRITTEN OUT IN      
!                 WFTAPE) INTO CARTESIAN SYSTEM                         

! !REVISION HISTORY:
!     Taken from SRC\_lapw2/latgen.f of wien2k\_04.10
!     Updated November 2004 (Kevin Jorissen)
!EOP

!
! KJ : This file was taken from SRC_lapw2 in the wien2k_04.10 version.
! KJ : I have changed the use statements, since I have other modules than in lapw2.
! KJ : The matrix br2, which was also generated here, has been removed.
! KJ : All output statements were added by me.
!
      IMPLICIT NONE
!  LOCAL VARIABLES :
	  real*8 aa,bb,cc,sqrt3,pia(3),vol,rvfac,sinab,sinbc,wurzel,cosab,cosac,cosbc,alfa(3)
      real*8 nul,een
	  integer i,j

      write(6,'(/,a)') 'Output from subroutine latgen :'
!---------------------------------------------------------------------  
      SQRT3=DSQRT(dble(3))

	  nul=dble(0)
	  een=dble(1)
      ALFA(1)=ALPHA(1)*PI/dble(180)                                             
      ALFA(2)=ALPHA(2)*PI/dble(180)                                             
      ALFA(3)=ALPHA(3)*PI/dble(180)                                             
      aa=alat(1)
	  bb=alat(2)
	  cc=alat(3)
      PIA(1)=dble(2)*PI/AA                                                 
      PIA(2)=dble(2)*PI/BB                                                 
      PIA(3)=dble(2)*PI/CC                                                 
      IF(LATTIC(1:1).EQ.'H') GOTO 10                                    
      IF(LATTIC(1:1).EQ.'S') GOTO 20                                    
      IF(LATTIC(1:1).EQ.'P') GOTO 20                                    
      IF(LATTIC(1:1).EQ.'F') GOTO 30                                    
      IF(LATTIC(1:1).EQ.'B') GOTO 40                                    
      IF(LATTIC(1:1).EQ.'C') GOTO 50                                    
      IF(LATTIC(1:1).EQ.'R') GOTO 60                                    

!        Error: wrong lattice, stop execution

         GOTO 900
                                                                       
!.....HEXAGONAL LATTICE                                                 
 10   CONTINUE                                                          
      BR1(1,1)=dble(2)/SQRT3*PIA(1)                                        
      BR1(1,2)=een/SQRT3*PIA(1)                                        
      BR1(1,3)=nul
      BR1(2,1)=nul
      BR1(2,2)=PIA(2)                                                   
      BR1(2,3)=nul
      BR1(3,1)=nul
      BR1(3,2)=nul
      BR1(3,3)=PIA(3)                                                   
                                                                       
      RVFAC=dble(2)/DSQRT(dble(3))                                             
      GOTO 100                                                          
!                                                                       
!.....RHOMBOHEDRAL CASE                                                    
 60   BR1(1,1)=een/DSQRT(dble(3))*PIA(1)                                          
      BR1(1,2)=een/DSQRT(dble(3))*PIA(1)                                          
      BR1(1,3)=-dble(2)/DSQRT(dble(3))*PIA(1)                                         
      BR1(2,1)=-PIA(2)                                                  
      BR1(2,2)= PIA(2)                                                    
      BR1(2,3)=nul                                                    
      BR1(3,1)=PIA(3)                                                    
      BR1(3,2)=PIA(3)                                                    
      BR1(3,3)=PIA(3)                                                    
      RVFAC=dble(6)/DSQRT(dble(3))
      GOTO 100                                                          
!                                                                       
!.....PRIMITIVE LATTICE                                                 
!                                                                       
  20  CONTINUE
      SINBC=DSIN(ALFA(1))
      COSAB=DCOS(ALFA(3))
      COSAC=DCOS(ALFA(2))
      COSBC=DCOS(ALFA(1))
      WURZEL=DSQRT(SINBC**2-COSAC**2-COSAB**2+2*COSBC*COSAC*COSAB)
      BR1(1,1)= SINBC/WURZEL*PIA(1)
      BR1(1,2)= (-COSAB+COSBC*COSAC)/(SINBC*WURZEL)*PIA(2)
      BR1(1,3)= (COSBC*COSAB-COSAC)/(SINBC*WURZEL)*PIA(3)
      BR1(2,1)= nul
      BR1(2,2)= PIA(2)/SINBC
      BR1(2,3)= -PIA(3)*COSBC/SINBC
      BR1(3,1)= nul
      BR1(3,2)= nul
      BR1(3,3)= PIA(3)
      RVFAC= een/WURZEL
      GOTO 100
!                                                                       
!.....FC LATTICE                                                        
 30   CONTINUE                                                          
      BR1(1,1)=PIA(1)                                                   
      BR1(1,2)=nul                                                    
      BR1(1,3)=nul                                                    
      BR1(2,1)=nul
      BR1(2,2)=PIA(2)                                                   
      BR1(2,3)=nul
      BR1(3,2)=nul
      BR1(3,1)=nul
      BR1(3,3)=PIA(3)                                                   
                                                                       
      RVFAC=dble(4)                                                        
      GOTO 100                                                          
!                                                                       
!.....BC LATTICE                                                        
 40   CONTINUE                                                          
      BR1(1,1)=PIA(1)                                                   
      BR1(1,2)=nul
      BR1(1,3)=nul
      BR1(2,1)=nul
      BR1(2,2)=PIA(2)                                                   
      BR1(2,3)=nul
      BR1(3,1)=nul
      BR1(3,2)=nul
      BR1(3,3)=PIA(3)                                                   
      RVFAC=dble(2)
      GOTO 100                                                    
!             
 50   CONTINUE                                                          
      IF(LATTIC(2:3).EQ.'XZ') GOTO 51                                    
      IF(LATTIC(2:3).EQ.'YZ') GOTO 52                                    
!.....CXY LATTICE                                                          
      BR1(1,1)=PIA(1)                                                   
      BR1(1,2)=nul
      BR1(1,3)=nul
      BR1(2,1)=nul
      BR1(2,2)=PIA(2)                                                   
      BR1(2,3)=nul
      BR1(3,1)=nul
      BR1(3,2)=nul
      BR1(3,3)=PIA(3)                                                   
      RVFAC=dble(2)
      GOTO 100                                                          
!                                                                       
!.....CXZ CASE (CXZ LATTICE BUILD UP)                                     
 51   CONTINUE                                     
!.....CXZ ORTHOROMBIC CASE 
      IF(ABS(ALPHA(3)-PI/2.0D0).LT.0.0001) then
         BR1(1,1)=PIA(1)                                                   
         BR1(1,2)=nul
         BR1(1,3)=nul
         BR1(2,1)=nul
         BR1(2,2)=PIA(2)                                                   
         BR1(2,3)=nul
         BR1(3,1)=nul
         BR1(3,2)=nul
         BR1(3,3)=PIA(3)                                                   
         RVFAC=dble(2)
         GOTO 100                                                          
      ELSE
!.....CXZ MONOCLINIC CASE 
         write(6,*) '  gamma not equal 90'
         SINAB=DSIN(ALFA(3))
         COSAB=DCOS(ALFA(3))
!                                                                       
         BR1(1,1)= PIA(1)/SINAB 
         BR1(1,2)= -PIA(2)*COSAB/SINAB
         BR1(1,3)= nul
         BR1(2,1)= nul
         BR1(2,2)= PIA(2)                                                     
         BR1(2,3)= nul
         BR1(3,1)= nul
         BR1(3,2)= nul
         BR1(3,3)= PIA(3)                                                     
         RVFAC=dble(2)/SINAB                                                   
         GOTO 100                                                          
      ENDIF
!                                                                       
!.....CYZ CASE (CYZ LATTICE BUILD UP)                                     
 52   CONTINUE                                     
      BR1(1,1)=PIA(1)                                                   
      BR1(1,2)=nul
      BR1(1,3)=nul
      BR1(2,1)=nul
      BR1(2,2)=PIA(2)                                                   
      BR1(2,3)=nul
      BR1(3,1)=nul
      BR1(3,2)=nul
      BR1(3,3)=PIA(3)                                                   
      RVFAC=dble(2)
      GOTO 100                                                          
!                                                                       
!.....DEFINE VOLUME OF UNIT CELL                                        
 100  CONTINUE                                                          
      VOL=AA*BB*CC/RVFAC                                                

      if (verbosity.ge.1) then
	    write(6,'(a,x,a)') 'lattice type ',lattic
		write(6,'(a,4x,3(F10.5,2x))') 'lattice parameters in atomic units ',aa,bb,cc
		write(6,'(a,4x,3(F10.3,2x))') 'lattice angles in degrees ',(alpha(j),j=1,3)
		write(6,'(a,4x,F10.5)') 'volume of (primitive?) unit cell in atomic units ^3 ',vol
		write(6,'(a)') 'Bravais matrix '
		do j=1,3
		  write(6,'(4x,3(F10.5,2x))') (br1(j,i),i=1,3)
		enddo
	  endif


                                                                       
      RETURN

!        Error messages

  900 CALL OUTERR('LATGEN','wrong lattice.')
      STOP 'LATGEN - Error'
      END                                                               
