!BOP
! !ROUTINE: ProductMatMat
! !INTERFACE:
      SUBROUTINE POTSL (DV,D,DP,DR,DPAS,DEXV,Z,NP,ION)             
! !INPUT/OUTPUT PARAMETERS:
!   dv   : potential
!   d    : density
!   dp   : work array
!   dr   : radial mesh
!   dpas : step size of radial mesh
!   dexv : multiplicative constant for exchange
!   z    : atomic number
!   np   : number of radial points
!   ion  : z-number of electrons
! !DESCRIPTION:
! POTENTIEL INTEGRATION PAR UNE METHODE A 4 POINTS
! !REVISION HISTORY:
!     Taken from SRC_lcore
!     Updated November 2004 (Kevin Jorissen)
!EOP


      IMPLICIT REAL*8 (A-H,O-Z)
	  real*8 ion
      DIMENSION DV(*),D(*),DP(*),DR(*)                                  
      DAS=DPAS/24.                                                      
      DO I=1,NP                                                       
        DV(I)=D(I)*DR(I)
      enddo                                                  
      DLO=EXP(DPAS)                                                     
      DLO2=DLO*DLO                                                      
      DP(2)=DR(1)*(D(2)-D(1)*DLO2)/(12.*(DLO-1.))                       
      DP(1)=DV(1)/3.-DP(2)/DLO2                                         
      DP(2)=DV(2)/3.-DP(2)*DLO2                                         
      J=NP-1                                                            
      DO I=3,J                                                        
        DP(I)=DP(I-1)+DAS*(13.*(DV(I)+DV(I-1))-(DV(I-2)+DV(I+1)))         
      enddo
      DP(NP)=DP(J)                                                      
      DV(J)=DP(J)                                                       
      DV(NP)=DP(J)                                                      
      DO I=3,J                                                        
        K=NP+1-I                                                          
        DV(K)=DV(K+1)/DLO+DAS*(13.*(DP(K+1)/DLO+DP(K))-(DP(K+2)/DLO2+DP(K-1)*DLO))
	  enddo
      DV(1)=DV(3)/DLO2+DPAS*(DP(1)+4.*DP(2)/DLO+DP(3)/DLO2)/3.          
      DLO=-(ION+1)                                                      
      DO I=1,NP                                                       
        DV(I)=DV(I)-(Z+3.*DEXV*((DR(I)*D(I)/105.27578)**(1./3.)))         
      enddo                                                 
      RETURN                                                            
      END                                                               
