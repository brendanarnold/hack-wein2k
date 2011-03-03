SUBROUTINE qmix5(clmnew,clmold,roknew,rokold,jspin,nq3,nat,lmmax,jri,qmx, &
                 ndm,jatom1,ll,alx,aly,alz,dmat,alx_old,aly_old,alz_old,dmat_old, &
                 mult,inversion,scl1,scl2,nbroyd,nmax,LM)
  IMPLICIT REAL*8 (A-H,O-Z)
  INCLUDE 'param.inc'
  !                                                                       
  INTEGER, INTENT(IN)       :: ndm,jatom1(nat,4,3),ll(nat,4,3)
  REAL*8, INTENT(IN)        :: alx_old(nat,4,3),aly_old(nat,4,3),alz_old(nat,4,3)
  COMPLEX*16, INTENT(IN)    :: dmat_old(nat,4,-3:3,-3:3,3)
  REAL*8, INTENT(INOUT)     :: alx(nat,4,3),aly(nat,4,3),alz(nat,4,3)
  COMPLEX*16, INTENT(INOUT) :: dmat(nat,4,-3:3,-3:3,3)
  COMPLEX*16 ROKOLD,ROKNEW
  INTEGER LM(2,NCOM,NAT)                                          
  !                                                                       
  real*8,allocatable :: YB(:),PS(:,:),SB(:),FN1(:),T1(:),VT(:),UI(:)
!  real*8,allocatable :: SCALE(:),SCOLD(:)
  DIMENSION CLMOLD(NRAD,NCOM,NAT,2),CLMNEW(NRAD,NCOM,NAT,2),       &
       ROKNEW(nq3,2),ROKOLD(nq3,2)
  INTEGER mult(*)                                 
  DIMENSION LMMAX(NAT),JRI(NAT)
  logical inversion                                   
  !                                                                       
  !    CHARGE DENSITY MIXING USING BROYDENS METHOD                        
  !    PROGRAM MODIFIED BY D. SINGH  DECEMBER 1985                        
  !    PROGRAM MODIFIED BY P. BLAHA  JANUARY  1987                        
  !    PROGRAM MODIFIED BY G. MADSEN MAY      2002
  !                                                                       
  !    Modifications by L. D. Marks
  !     July 2004
  !             Added seperate scalings for plane waves & CLM
  !             Added inversion trapping (does it matter?)
  !             Added experimental minimum cycle routine - not good
  !             Corrected code for first cycle
  !             Corrected error in density matrix mixing
  !     August 2004
  !             Added Harmonic scaling -- and removed it (not viable)
  !             Added monitor of Broyd/Pratt step size and relative angle
  !             Added autorestart of Broyd if Angle > Acut
  !             Added debug output of scales to mixer.scale, not currently active
  DATA ZERO/0.D0/,ONE/1.D0/,ACUT/138.6D0/       
  MAXMQ=NQ3
  idmat=0
  DO idm=1,ndm
     DO  JATOM=1,NAT
        IF(jatom1(jatom,1,idm).EQ.0) EXIT
        DO iorb=1,4
           IF(jatom1(jatom,iorb,idm).EQ.0) EXIT      
           l=ll(jatom,iorb,idm)
           idmat=idmat+3+(l*l+1)*(l*l+1)*2
        ENDDO
     ENDDO
  ENDDO
!  NBSIZE=2*nq3*2+NRAD*NCOM*NAT*2+idmat
  NBSIZE=jspin*nq3*2+NRAD*SUM(LMMAX(1:NAT))*jspin+idmat
  allocate ( YB(NBSIZE),PS(NBSIZE,2),                              &
       SB(NBSIZE),FN1(NBSIZE),T1(NBSIZE),VT(NBSIZE),UI(NBSIZE))
!  allocate (SCALE(NBSIZE),SCOLD(NBSIZE))
!
! These are the plane wave coefficients
! They are small, hence badly scaled against the CLM !
  scl_plane=scl1
  Nplane=0
  Splane=0
  DO ISPIN=1,JSPIN                                                
     K0=2*MAXMQ*(ISPIN-1)                                              
     DO K=1,MAXMQ                                                    
        PS(2*K-1+K0,2)=dble(ROKNEW(K,ISPIN))*scl_plane                              
        PS(2*K-1+K0,1)=dble(ROKOLD(K,ISPIN))*scl_plane                              
        Nplane=Nplane+1
!        Splane=Splane+abs(PS(2*K-1+K0,2)-PS(2*K-1+K0,1))
        tmp=PS(2*K-1+K0,2)-PS(2*K-1+K0,1)
        Splane=Splane+tmp*tmp
        if(inversion)then
          PS(2*K+K0,1)=0
          PS(2*K+K0,2)=0
        else
          Nplane=Nplane+1
          PS(2*K+K0,2)=DIMAG(ROKNEW(K,ISPIN))*scl_plane  
          PS(2*K+K0,1)=DIMAG(ROKOLD(K,ISPIN))*scl_plane                               
!          Splane=Splane+abs(PS(2*K+K0,2)-PS(2*K+K0,1))
        tmp=PS(2*K+K0,2)-PS(2*K+K0,1)
        Splane=Splane+tmp*tmp
        endif
!        scale(2*K+K0-1)=scl_plane
!        scale(2*K+K0)=scl_plane
     ENDDO
  ENDDO
!  I believe Splane is the mean-square distance in real space
!  The factor of 100 puts it on a comparable scale to :DIS, :ENE
   if(Splane.gt.1d-20)Splane=sqrt(Splane)*100
   write(21,210)':PLANE:  INTERSTITAL DISTANCE  ',Splane
!  Write(21,210)':PLANE:  INTERSTITAL DISTANCE  ',Splane/Nplane*1d5
!  Write(21,*)'     SCALES',scl1,scl2
  210 FORMAT(A,F9.7)
  IS=2*MAXMQ*JSPIN                                                  
  !
  PI=2.d0*asin(1.D0)                                                                       
!  sqfp=SQRT(4.d0*pi)
  sqfp=1.D0
! The L=0 term is divided by sqfp -- multiply others ?
  scl_clm=scl2
!  lmult=1
  DO ISPIN=1,JSPIN                                              
     DO N=1,NAT                                                     
        JRIN=JRI(N)                                                       
        NH=LMMAX(N)
        sclm=scl_clm*dble(Mult(N))
!        write(6,*)'Norm ',n,sclm                                                       
        DO LH=1,NH
!           Lmult=abs(LM(1,LH,N))
!           if(LH.gt.1)sclm=scl_clm*sqfp*Lmult*dble(Mult(N))
           DO J=1,JRIN                                                    
              IS=IS+1
!              scale(is)=sclm/mult(n)
              PS(IS,1)=CLMOLD(J,LH,N,ISPIN)*sclm
              PS(IS,2)=CLMNEW(J,LH,N,ISPIN)*sclm
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  DO idm=1,ndm
     DO  JATOM=1,NAT
        IF(jatom1(jatom,1,idm).EQ.0) EXIT
        sclm=dble(mult(jatom))
!        write(6,*)'Norm 2 ',sclm
        DO iorb=1,4
           IF(jatom1(jatom,iorb,idm).EQ.0) EXIT      
           is=is+3
           ps(is-2,1)=alx_old(jatom,iorb,idm)*sclm
           ps(is-1,1)=aly_old(jatom,iorb,idm)*sclm
           ps(is,1)  =alz_old(jatom,iorb,idm)*sclm
           ps(is-2,2)=alx(jatom,iorb,idm)*sclm
           ps(is-1,2)=aly(jatom,iorb,idm)*sclm
           ps(is,2)  =alz(jatom,iorb,idm)*sclm
           l=ll(jatom,iorb,idm)
!           scale(is-2)=sclm/mult(jatom)
!           scale(is-1)=sclm/mult(jatom)
!           scale(is)=sclm/mult(jatom)
           DO m=-l,l
              DO mp=-l,l
                 is=is+2
                 ps(is-1,2)=DBLE(dmat(jatom,iorb,m,mp,idm))*sclm
                 ps(is-1,1)=DBLE(dmat_old(jatom,iorb,m,mp,idm))*sclm
!       Should these be zero for inversion centering ???
!                 if(inversion)then
!                  ps(is,2)=0
!                  ps(is,1)=0
!                 else
                  ps(is,2)  =DIMAG(dmat(jatom,iorb,m,mp,idm))*sclm
                  ps(is,1)  =DIMAG(dmat_old(jatom,iorb,m,mp,idm))*sclm
!                 endif
!                  scale(is-1)=sclm
!                  scale(is)=sclm/mult(jatom)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !                                                                       
  MAXMIX=IS
  !  Scale -- for debug
!  open(unit=88,file='mixer.scale')
!  do j=1,maxmix
!     scold(j)=0
!  enddo
!  igots=0
!  do j=1,maxmix
!     read(88,*,err=88,end=88)kk,scold(j)
!  enddo
!  igots=1
!  88 rewind 88
  WRITE( 6,1001)MAXMQ,MAXMIX                                        
  IF(MAXMIX.GT.NBSIZE) STOP 'QMIX5:NBSIZE'                            
  !                                                                       
  !     PERFORM MIXING ON CHARGE DENSITIES QPW AND FITTED RHO             
  !     USE BROYDENS METHOD                                               
  !                                                                       
!  REWIND(74)
!  READ(74,*,err=119,end=119) idum
! PRATT mixing for first or more than nbroyd iterations
  READ(31,END=119) DMIX,LASTIT
  if(lastit.gt.nbroyd) goto 119
!                                       
!  DMIX=QMX                                                          
  if(qmx.lt.dmix) then
    dmix=qmx
  else
! Smaller of the average and twice the prior step
! Rather a classic "Trust Region" approach
    DMIX=Min(2.d0*dmix,(dmix+qmx)*0.5d0,qmx)
  endif
  QMX=DMIX                                                         
  READ(31)(FN1(K),K=1,MAXMIX)                                       
  READ(31)(SB(K),K=1,MAXMIX)                                        
  YNORM=ZERO                                                        
  DO K=1,MAXMIX
  ! Change in position
     SB(K)=PS(K,1)-SB(K)                                               
     TEMP=PS(K,2)-PS(K,1)
  ! Change in distance
     YB(K)=TEMP-FN1(K)
  ! Scaling stuff
  !   tt1=SB(k)*YB(K)/(YB(K)*YB(K)+1d-24)
  !   tt1=abs(tt1)
  !   tt2=scale(k)
  !   if(igots.eq.0)then
  !       scale(k)=tt1
  !       scold(k)=tt1
  !   else
  !       scale(k)=(tt1+2.d0*scold(k))/3.d0
  !   endif
  !   write(88,881)k,scale(k),scold(k),tt1,tt2,PS(K,1),PS(K,2)
  !   881 format(i6,6d11.2)
     YNORM=YNORM+YB(K)**2                                              
     FN1(K)=TEMP                                                       
  ENDDO
  LASTIT=LASTIT+1                                                   
  REWIND 31                                                         
  REWIND 32                                                         
  WRITE(31)DMIX,LASTIT                                              
  WRITE(31)(FN1(K),K=1,MAXMIX)                                      
  WRITE(31)(PS(K,1),K=1,MAXMIX)
  DO K=1,MAXMIX                                                 
     UI(K)=DMIX*YB(K)+SB(K)                                            
     VT(K)=YB(K)/YNORM                                                 
  ENDDO
  ! Only use at most nmax cycles
  nskip=LASTIT-1-nmax
  nuse=min(nmax,lastit-1)
  IF(LASTIT.NE.2)THEN                                               
     JMAX=LASTIT-1                                                    
  if(nskip.gt.0)then
        do j=1,nskip
                READ(32)(SB(K),K=1,MAXMIX)
                READ(32)(T1(K),K=1,MAXMIX)
        enddo
  endif
  DO J=2,NUSE
  !     DO J=2,JMAX
        READ(32)(SB(K),K=1,MAXMIX)                                       
        READ(32)(T1(K),K=1,MAXMIX)                                       
        AKJ=ZERO                                                         
        DO K=1,MAXMIX                                                 
           AKJ=AKJ+T1(K)*YB(K)                                              
        ENDDO
        DO  K=1,MAXMIX                                                 
           UI(K)=UI(K)-AKJ*SB(K)                                            
        ENDDO
     ENDDO
  ENDIF
  WRITE(32)(UI(K),K=1,MAXMIX)                                       
  WRITE(32)(VT(K),K=1,MAXMIX)                                       
  REWIND 32                                                         
  XX=(ONE-DMIX)                                                     
  YY=DMIX                                                          
  DO K=1,MAXMIX                                                  
     SB(K)=XX*PS(K,1)+YY*PS(K,2)                                       
  ENDDO
  if(nskip.gt.0)then
        do j=1,nskip
                READ(32)(UI(K),K=1,MAXMIX)                                        
                READ(32)(VT(K),K=1,MAXMIX)     
        enddo
  endif
  DO J=1,NUSE                                                  
     READ(32)(UI(K),K=1,MAXMIX)                                        
     READ(32)(VT(K),K=1,MAXMIX)                                        
     CKI=ZERO                                                          
     DO K=1,MAXMIX                                                  
        CKI=CKI+VT(K)*FN1(K)                                              
     ENDDO
     DO K=1,MAXMIX                                                  
        SB(K)=SB(K)-CKI*UI(K)                                             
     ENDDO
  ENDDO
!
! Consistency check
        DOT=0
        BMOD=0
        PMOD=0
!        open(unit=99,file='mixer.debug')
        DO K=1,MAXMIX
                SPRATT=(PS(K,2)-PS(K,1))*DMIX
                SBROY=SB(K)-PS(K,1)
!                write(99,99)K,SPRATT,SBROY,SB(K),PS(K,1),PS(K,2)
!99              format(i6,5D13.5)
                DOT=SPRATT*SBROY+DOT
                BMOD=BMOD+SBROY*SBROY
                PMOD=PMOD+SPRATT*SPRATT
        ENDDO
!        IF(DOT.LT.0.0)THEN
!               Moving opposite to Pratt direction
!               Broyden matrix may have collapsed
                BMOD=sqrt(BMOD)
                PMOD=sqrt(PMOD)
                DTT=acos(DOT/(BMOD*PMOD))*90./asin(1.D0)
                Write(21,98)BMOD,PMOD,DTT
                IF(abs(DTT).GT.ACUT.and.lastit.ge.6)Then
                        QMX=DMIX*MIN(BMOD,PMOD*0.5D0)/PMOD
                        QMX=max(QMX,0.01D0)
                        Goto 119
                endif
98      format(':DIRB :  |BROYD|=',D10.3,' |PRATT|=',D10.3, &
        ' ANGLE=',f6.1,' DEGREES')   
!                GOTO 119
!        ENDIF
100  DO ISPIN=1,JSPIN                                              
     K0=MAXMQ*2*(ISPIN-1)                                              
     DO K=1,MAXMQ
        if(inversion)then
        ROKNEW(K,ISPIN)=SB(2*K-1+K0)/scl_plane
        else                                                   
        ROKNEW(K,ISPIN)=CMPLX(SB(2*K-1+K0),SB(2*K+K0))/scl_plane               
        endif
     ENDDO
  ENDDO
  IS=2*MAXMQ*JSPIN                                                  
  DO ISPIN=1,JSPIN                                              
     DO N=1,NAT                                                     
        JRIN=JRI(N)                                                       
        NH=LMMAX(N)
        sclm=1.D0/(scl_clm*mult(n))
!        write(6,*)'Renorm ',N,sclm                                                       
        DO LH=1,NH
!           Lmult=abs(LM(1,LH,N))
!           if(LH.gt.1)sclm=1.D0/(scl_clm*sqfp*Lmult*mult(n))

           DO J=1,JRIN                                                    
              IS=IS+1
              CLMNEW(J,LH,N,ISPIN)=SB(IS)*sclm
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  DO idm=1,ndm
     DO  JATOM=1,NAT
        IF(jatom1(jatom,1,idm).EQ.0) EXIT
        sclm=1.d0/dble(mult(jatom))
!        write(6,*)'Renorm ',sclm
        DO iorb=1,4
           IF(jatom1(jatom,iorb,idm).EQ.0) EXIT      
           is=is+3
           alx(jatom,iorb,idm)=sb(is-2)*sclm
           aly(jatom,iorb,idm)=sb(is-1)*sclm
           alz(jatom,iorb,idm)=sb(is)*sclm
           l=ll(jatom,iorb,idm)
           DO m=-l,l
              DO mp=-l,l
                 is=is+2
!       Is this right?
!                 if(inversion)then
!                  dmat(jatom,iorb,m,mp,idm)=sb(is-1)*sclm
!                 else
                  dmat(jatom,iorb,m,mp,idm)=CMPLX(sb(is-1),sb(is))*sclm
!                 endif
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  REWIND 31                                                         
  REWIND 32                                                         
  GOTO 120                                                          
! ... PRATT MIXING FOR FIRST ITERATION                                  
119 REWIND 31
  LASTIT=1 
  DMIX=min(QMX,0.1d0) 
  WRITE(31)DMIX,LASTIT                                                         
  WRITE(21,*)':Pratt iteration over-ride, DMIX ',DMIX 
  QMX=DMIX
!  WRITE(31)DMIX,LASTIT                                           
  DO K=1,MAXMIX                                                 
     FN1(K)=PS(K,2)-PS(K,1)
     SB(K)=PS(K,1)+DMIX*FN1(K)                                            
  ENDDO
  WRITE(31)(FN1(K),K=1,MAXMIX)                                      
  WRITE(31)(PS(K,1),K=1,MAXMIX)
  GOTO 100                                     

120 CONTINUE                                                          
  deallocate ( YB,PS,SB,FN1,T1,VT,UI )
  RETURN                                                            
!                                                                       
1001 FORMAT(/,' QMIX5:',/,8X,'MAXMQ =',I8,'   MAXMIX =',I10,/)          
1022 FORMAT(10X,'DENSITY FOR ITERATION',I4,' PREPARED')                
END SUBROUTINE qmix5
