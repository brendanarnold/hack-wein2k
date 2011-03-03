!.......................................................................
      SUBROUTINE ADOS (NS,NFU,NEMIN,NEMAX,ebs,fc,nst,RNUMB,DENSTY,SNUMB,nlast)
!     **   ARBDOS READS EIGENVALUES AND EIGENFUNCTIONS FROM FILE 8    **
!     **                                                              **
!     **   SUBROUTINES USED:                                          **
!     **    NOCULC                                                    **
!     **                                                              **
!     **  INPUT AND OUTPUT:                                           **
!     **    NFU: NUMBER OF EIGENFUNCTIONS PER K-POINT                 * 
!     **    NNOC: NUMBER OF ALLREADY FILLED TETRAHEDRA                **
!     **    NEMAX: NUMBER OF LAST BAND                                **
!     **    NEMIN: NUMBER OF FIRST BAND-1                             **
!     **    NUMK: NUMBER OF K-POINTS                                  **
!     **   NTT  NUMBER OF DIFFERENT TETRAHEDRA                        **
!     **                                                              **
      INCLUDE 'param.inc'
      COMMON /NCOUNT/ NNOC                                              
      real*4 EBS(Nst,*), FC(MG,Nst,*)
      real*4 RNUMB(nlast,MF), DENSTY(nlast,MF)
      real*4 SNUMB(nlast,MF)                                        
      COMMON /EMICRO/ D(4), F(MF,4), V                                  
      PARAMETER (NITT = 505)
      DIMENSION ITTFL(nitt)                                              
      DATA IWTET /0/  
      NNOC=0                                                            
      NWRIT=100                                                         
      READ(3,1234)Ndim,NTT,V1,nWRIT,NREC                               
 1234 format(2i10,e20.12,2i10) 
!     READ(12,1101) NDIM,NTT,V1                                         
      WRITE (6,*) 'NUMBER OF K-POINTS:',NDIM
      write(6,*)  'NUMBER OF TETRAHEDRONS:',NTT                        
!      NTTW=(NTT/NWRIT)+1                                                
!.....correction by Nelhiebel
      if((ntt/nwrit)*nwrit.eq.ntt) then
        Nttw=NTT/NWRIT
      else
        Nttw=NTT/NWRIT+1
      end if
!.....
      DO 30 N=1,NTTW         
      IF (N.EQ.NTTW) THEN                                               
        KMAX=NTT-(N-1)*NWRIT                                            
      ELSE                                                              
        KMAX=NWRIT                                                      
      ENDIF                                                             
!cad
      if(5*kmax.gt.NITT) then
      write(*,*) 'increase parameter NITT in ados!'
      write(*,*) 'the value should be at least ',5*kmax
      STOP 'NITT in ados too small'
      endif
!cad
!     READ(12,1102) (ITTFL(K),K=1,KMAX*5)                               
      READ(3,1235) (ITTFL(K),K=1,KMAX*5)                               
 1235 format(6i10) 
      DO 30 K=1,KMAX          
      V=FLOAT(ITTFL((K-1)*5+1))*V1                                      
      DO 30 JJ=NEMIN,NEMAX                                              
      DO 20 KP=1,4                                                      
      KPP=ITTFL((K-1)*5+KP+1)                                           
      DO 10 I=1,NFU 
   10 F(I,KP)=FC(I,KPP,JJ)                                              
   20 D(KP)=EBS(KPP,JJ) 
      CALL NOCULC (NS, RNUMB,DENSTY,SNUMB,nlast)  
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               
