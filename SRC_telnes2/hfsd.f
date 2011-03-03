!BOP
! !ROUTINE: Hfsd
! !INTERFACE:
      SUBROUTINE HFSD(NATOM,NPRIM,LPRIM,KAPPA,denxx) 
! !USES:
      use initialstate1
	  use program_control,only : verbosity,calcsplit
	  use dimension_constants,only : nrad
	  use struct, only : rela => rel
! !INPUT/OUTPUT PARAMETERS:
!   denxx    :   energy of the core state in Ry  
!   natom    :   atom for which core state is calculated
!   nprim    :   main quantum number of core state 
!   lprim    :   orbital quantum number of core state
!   kappa    :   relativistic quantum number of core state
! !DESCRIPTION:
!   Main routine for the calculation of core electron wave functions.
!   Large and small component is calculated.
!   'Atomic eigenvalues of a spherical total potential are calculated.'
!   hfsd stands for Hartree Fock Dirac Slater
! !REVISION HISTORY:
!     J P DESCLAUX   CEA PARIS 1969
!     Taken from SRC_lcore.
!     Heavily modified November 2004 (Kevin Jorissen)
!EOP
      implicit none
!  INPUT VARIABLES
      INTEGER,intent(in) :: NATOM,NPRIM,LPRIM,KAPPA
	  real*8,intent(out) :: denxx
!  LOCAL VARIABLES
      CHARACTER *10 BAR
	  character titre                             
	  real*8 DEN,DQ1,DFL,DEXV,DEXE,DCOP,TEST,TESTE, &
	    TESTY,TESTV,DVN(NRAD+41),DVF(NRAD+41),D(NRAD+41),DC(NRAD+41), &
                    DGC(NRAD+41),DPC(NRAD+41),ion
	  integer NQN,NQL,NK,NMAX,NEL,NITER
	  integer i,natt,jatom,npxx,iter,j,n,imax,kk,np40,im,l,k
	  real*8 vd(nrad+41),ymax,vmax,emax,de,val,y,yn,dval
	  real*8,external :: dalp
	  real*8,allocatable :: dpsave(:),dqsave(:)
	  	                                                    
      COMMON DEN,DQ1,DFL,NQN,NQL,NK,NMAX,NEL                                                         
      COMMON/PS2/DEXV,DEXE,DCOP,TEST,TESTE,TESTY,TESTV,ION,NITER,TITRE
      COMMON /CHAR/ BAR                                                 
      COMMON /DEUX/ DVN,DVF,D,DC,DGC,DPC 



      write(6,'(/,a)') 'Output from subroutine hfsd:'

      NSTOP=1                                             
      CALL INSLD(RELA,npxx,natom,nprim,lprim,kappa)

      ITER=1
      DGC(1:NP)=dble(0)
      DPC(1:NP)=dble(0)
	  if(verbosity.ge.2) WRITE(6,1002)
      N=-(int(ION)+1)  ! int added since ion is now real


      do iter=1,niter

        vd(1:NP)=dble(0)
        D(1:NP)=dble(0)
        TETS=TEST
        YMAX=dble(0)
        VMAX=dble(0)
        EMAX=dble(0)
! RESOLUTION DE L EQUATION DE DIRAC POUR CHAQUE ORBITALE                
! **********************************************************************
        DE=DEN

 16     CALL RESLD (NQN,NQL,NK,IMAX,DEN,DFL,DQ1,1,RELA) 
        IF (NSTOP.EQ.0) GO TO 18                                          
        IF (NSTOP.NE.362.OR.ITER.GE.10.OR.TETS.GT.TEST) GO TO 17          
        TETS=TESTV                                                        
        GO TO 16                                                          
 17     if(verbosity.ge.2) WRITE (6,1010) NSTOP,NQN,TITRE                              
        stop 'error in lcore'
 18     VAL=ABS((DEN-DE)/DE)                                           
        IF(VAL.GT.EMAX) EMAX=VAL                                          
        NMAX=IMAX                                                      
        DO I=1,NP                                                     
          VAL=DGC(I)-DP(I)                                                
          IF (ABS(DP(I)).GT.1.) VAL=VAL/DP(I)                               
          IF (ABS(VAL).LT.ABS(YMAX)) GO TO 21                               
          YMAX=VAL                                                          
          Y=DP(I)                                                           
          YN=DGC(I)                                                       
 21       VAL=DPC(I)-DQ(I)                                                
          IF (ABS(DQ(I)).GT.1.) VAL=VAL/DQ(I)                               
          IF (ABS(VAL).LT.ABS(YMAX)) GO TO 23                               
          YMAX=VAL                                                          
          Y=DQ(I)                                                           
          YN=DPC(I)                                                       
 23       DGC(I)=DP(I)                                                    
          DPC(I)=DQ(I)                                                    
          D(I)=D(I)+NEL*(DP(I)*DP(I)+DQ(I)*DQ(I))
!     to get r*P(core,n'l') 
!
!      In DP und DQ are r*r*
          DP(I)=DR(I)*DP(I)
          DQ(I)=DR(I)*DQ(I)

        ENDDO ! i=1,np

        if(.not.calcsplit) return !  This is a succesful exit.
!        Apparently, hfsd reuses the dp and dq arrays for something, so I save them first :
        if(iter.eq.1) then
	      allocate(dqsave(nrad+41),dpsave(nrad+41))
          dpsave=dp
	      dqsave=dq
	    endif

      
        CALL POTSL (DC,D,DP,DR,DPAS,DEXV,Z,NP,ION)                   
        DO I=1,NP   ! used to be do 33
          DVAL=ABS(DC(I)-DV(I))                                             
          IF ((DR(I)*DC(I)).LE.N) DVAL=-DVAL/DC(I)                          
          IF (DVAL.LE.VMAX) cycle !GO TO 33                                        
          VMAX=DVAL                                                         
          J=I                                                               
        enddo          ! used to be label 33 continue
	    if(verbosity.ge.2) WRITE(6,1003) ITER,VMAX,DR(J),DV(J),DC(J),EMAX,YMAX,YN,Y

        IF(TETS.LE.TEST.AND.EMAX.LE.TESTE.AND.VMAX.LE.TESTV.AND.YMAX.LE.TESTY) exit  ! calculation is converged

! no convergence : prepare next iteration
! POTENTIEL POUR L ITERATION SUIVANTE                                   
        if(ITER.EQ.2) then 
          DVAL=dble(1)-DCOP
          do I=1,NP
            DVN(I)=DV(I)
            DVF(I)=DC(I)
            DV(I)=DVAL*DV(I)+DCOP*DC(I)
          enddo
        else
          DO  I=1,NP
            DVAL=DALP(DVN(I),DVF(I),DV(I),DC(I))
            DVN(I)=DV(I)
            DVF(I)=DC(I)
            DV(I)=DVAL*DV(I)+(1.-DVAL)*DC(I)
          enddo
        endif

      enddo  ! iter=1,niter 

! calculation finished.


      if(iter.gt.niter) WRITE(6,1004) NITER                                               
      NSTOP=2                                                           

      JATOM=JATOM+1
      NP40=NP-npxx
! VALEURS MOYENNES DE R
! **********************************************************************
      DVF(1:NP)=DC(1:NP)
      DQ(1:NP)=dble(0)
      DVAL=dble(0)                                                           
        IM=NMAX                                                        
        DVAL=DVAL+NEL*DEN                                           
        do  J=1,IM                                                     
          DC(J)=DGC(J)*DGC(J)+DPC(J)*DPC(J)
        enddo
        L=5                                                               
        IF (IABS(NK).EQ.1) L=L-1                                       
        do J=1,L                                                       
          DP(J)=DFL+DFL
		  if (j.lt.2) then
		    n=4
		  elseif (j.eq.2) then
		    n=2
		  elseif (j.eq.3) then
		    n=1
		  elseif(j.eq.4) then
		    n=-1
		  elseif (j.gt.4) then
		    n=-3
		  endif
          CALL SOMM (DR,DC,DQ,DPAS,DP(J),N,IM)
        enddo
!.... CONVERT TO RYD                                                    
        DENXX=DEN*dble(2)
		if(verbosity.ge.2) write (6,1006) NQN,TITRE,DEN,(DP(J),J=1,L)                

! ENERGIE TOTALE EN MOYENNE SPHERIQUE                                   
! **********************************************************************
      DC(1)=1                                                           
      do I=1,NP
        DP(I)=D(I)/DR(I)
      enddo
      CALL SOMM (DR,DP,DQ,DPAS,DC(1),0,NP)                              
      do I=1,NP                                                    
        DP(I)=D(I)*DVF(I)                                                 
        D(I)=D(I)*((D(I)*DR(I))**(dble(1)/dble(3)))
      enddo                                 
      DC(2)=3                                                           
      DC(3)=1                                                           
      CALL SOMM (DR,DP,DQ,DPAS,DC(3),0,NP)                              
      CALL SOMM (DR,D,DQ,DPAS,DC(2),-1,NP)                              
      DC(2)=-3.*DC(2)/(105.27578**(1./3.))                              
      DC(1)=-Z*DC(1)                                                    
      DC(4)=DVAL-DC(3)                                                  
      DVAL=DVAL+(DC(1)-DC(3)+(DEXE-DEXV)*DC(2))/dble(2)
      DC(3)=(DC(3)-DC(1)-DEXV*DC(2))/2.                                 
      DC(2)=DC(2)*DEXE/2.                                               
      if(verbosity.ge.2) WRITE (6,1007) DVAL,DC(4),DC(3),DC(2),DC(1)                       

      dq=dqsave
	  dp=dpsave
	  deallocate(dqsave,dpsave)
      RETURN  !  This is a succesful exit.



 1002 FORMAT (' ITER',4X,'DVMAX',10X,'R',14X,'VN-1',13X,'VN',10X,'DEMAX' &
      ,6X,'DPMAX',9X,'PN-1',13X,'PN')                                   
 1003 FORMAT (I5,1PE11.2,3(1PE16.6),2(1PE11.2),2(1PE16.6))              
 1004 FORMAT (' NOMBRE D ITERATIONS SUPERIEUR A',I4)                    
 1005 FORMAT (/,12X,'ENERGIE',12X,'(R4)',14X,'(R2)',14X,'(R)',15X,'(R-1)', &
      13X,'(R-3)')                                                     
 1006 FORMAT (I3,A2,6(1PE18.7))                                         
 1007 FORMAT (6X,'ET=',1PE14.7,5X,'EC=',1PE14.7,5X,'EE=',1PE14.7,5X, &
      'EX=',1PE14.7,5X,'EN=',1PE14.7)                                   
 1010 FORMAT ('  NSTOP=',I4,'  POUR L ORBITALE',I3,A2)                  

      END                                                               
