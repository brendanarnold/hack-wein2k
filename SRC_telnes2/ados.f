!BOP
! !ROUTINE: ADos
! !INTERFACE:
      SUBROUTINE ADOS (NS,NFU,NEMIN,NEMAX,ebs,fc,nst,RNUMB,DENSTY,SNUMB,nlast)
! !USES:
       use tetra_params  
! !INPUT/OUTPUT PARAMETERS:
!     NFU: NUMBER OF EIGENFUNCTIONS PER K-POINT       
!     NNOC: NUMBER OF ALLREADY FILLED TETRAHEDRA      
!     NEMAX: NUMBER OF LAST BAND                      
!     NEMIN: NUMBER OF FIRST BAND-1                   
!     NUMK: NUMBER OF K-POINTS                        
!     NTT  NUMBER OF DIFFERENT TETRAHEDRA             
! !DESCRIPTION:
! !REVISION HISTORY:
!     Taken from SRC\_tetra.
!     Updated November 2004 (Kevin Jorissen)
!     Small changes February 2005 (Kevin Jorissen)
!EOP

      implicit none   !KJ
      COMMON /NCOUNT/ NNOC                                              

!KJ IN/OUTPUT
      integer ns,nfu,nemin,nemax,nst,nlast  !KJ
      real*8 EBS(Nst,*), FC(MG,Nst,*) !KJ
      real*8 RNUMB(nlast,mg), DENSTY(nlast,mg) !KJ
      real*8 SNUMB(nlast,mg)                       !KJ                  
!KJ LOCALS      

! KJ I destroyed the common because it objected to my modules :
!KJ      COMMON /EMICRO/ D(4), F(MF,4), V                                  
      real*8 d(4),f(mg,4),v  ! this line added KJ
	  real*8 v1  !KJ
	  integer nitt,ittfl,iwtet,nnoc,nwrit,ndim,ntt,nrec,nttw,n,kmax,k,i,jj,kp,kpp !KJ
! KJ end of my changes
      PARAMETER (NITT = 505)
      DIMENSION ITTFL(nitt)                                              
      DATA IWTET /0/  

      NNOC=0                                                            
      NWRIT=100                                                         
      READ(3,1234)Ndim,NTT,V1,nWRIT,NREC                               
 1234 format(2i10,e20.12,2i10) 
      WRITE(6,'(a,i7)')  'NUMBER OF K-POINTS:',NDIM
      write(6,'(a,i10)') 'NUMBER OF TETRAHEDRA:',NTT                        
      if((ntt/nwrit)*nwrit.eq.ntt) then
        Nttw=NTT/NWRIT
      else
        Nttw=NTT/NWRIT+1
      end if

      DO N=1,NTTW         
        IF (N.EQ.NTTW) THEN                                               
          KMAX=NTT-(N-1)*NWRIT                                            
        ELSE                                                              
          KMAX=NWRIT                                                      
        ENDIF                                                             

        if(5*kmax.gt.NITT) then
          write(*,*) 'increase parameter NITT in ados!'
          write(*,*) 'the value should be at least ',5*kmax
          STOP 'NITT in ados too small'
        endif

        READ(3,1235) (ITTFL(K),K=1,KMAX*5)                               
 1235   format(6i10) 
        DO K=1,KMAX          
          V=FLOAT(ITTFL((K-1)*5+1))*V1                                      
          DO JJ=NEMIN,NEMAX                                              
            DO KP=1,4                                                      
              KPP=ITTFL((K-1)*5+KP+1)                                           
              DO I=1,NFU 
                F(I,KP)=FC(I,KPP,JJ)                                              
              enddo
              D(KP)=EBS(KPP,JJ) 
            enddo
!KJ         old call:   CALL NOCULC (NS, RNUMB,DENSTY,SNUMB,nlast)  
            CALL NOCULC (NS, RNUMB,DENSTY,SNUMB,nlast,d,f,v)  ! KJ I added the 3 last arguments  
	      enddo
	    enddo
	  enddo

      RETURN                                                            
      END                                                               
