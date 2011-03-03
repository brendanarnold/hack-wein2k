     SUBROUTINE ATPAR (REL,NAT,JATOM,LATOM,cform,ZZ)                            
      use param
      use defs
      use struk
      use lo; USE atspdt, only: e=>el,p,dp,pe,dpe,pei
      USE parallel
      IMPLICIT REAL*8 (A-H,O-Z)
#ifdef Parallel
  INCLUDE 'mpif.h'
#endif
      CHARACTER*4      cform                                           
      LOGICAL          REL
      dimension        emist2(0:lomax,nloat)
!                                                                       
      COMMON /RADFU/   RRAD1(NRAD,0:LMAX2),RADE1(NRAD,0:LMAX2), &
                       RRAD2(NRAD,0:LMAX2),RADE2(NRAD,0:LMAX2)                
      COMMON /POTNLC/  VR(NRAD)           
      COMMON /UHELP/   A(NRAD),B(NRAD),AP(NRAD),BP(NRAD),AE(NRAD),       &
                       BE(NRAD)                                         
      common /normal/  pu1u1(0:lxdos,0:lxdos),pu1ue(0:lxdos,0:lxdos) &
                      ,pu1u2(0:lxdos,0:lxdos,nloat),pueue(0:lxdos,0:lxdos) &
                      ,pueu2(0:lxdos,0:lxdos,nloat),pu2u2(0:lxdos,nloat,0:lxdos,nloat)
!---------------------------------------------------------------------  
!.....READ TOTAL SPHERICAL POTENTIAL V(0,0) OF TAPE18=VSP               
!     NORM OF V(0,0)=V00(R)*R/SQRT(4.D0*PI)                             

      READ(18,1980)                                                     
      READ(18,2000) IDUMMY                                              
      READ(18,2031)                                                     
      READ(18,2022) ( VR(J), J=1,JRI(JATOM) )                           
      READ(18,2031)                                                     
      READ(18,2030)                                                     
      DO J=1,JRI(JATOM)                                             
         VR(J)=VR(J)/2.0D0                                                 
      enddo
!                                                                       
      nlo=0
      nlov=0
      nlon=0
! nlo #LO on this atom, nlov #LO up til now, nlon #LO left
      ilo(0:lmax2)=0
      DO i=1,jatom                                                   
         IF(myid.EQ.0) THEN
            READ(10) e
            READ(10) elo
         ENDIF
#ifdef Parallel
        CALL MPI_BCAST(e,lmax2+1,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(elo,(lmax2+1)*nloat,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif
        IF(i.EQ.jatom) THEN
           DO l=0,lmax2
              lapw(l)=.TRUE.
              IF(e(l).GT.150.) THEN
                 e(l)=e(l)-200.d+0
                 lapw(l)=.FALSE.
              ENDIF
           ENDDO
        ENDIF
        DO l = 0,lomax
           DO k=1,nloat
              loor(k,l)=.FALSE.
              rlo(k,l)=.FALSE.
              IF (i.EQ.jatom) THEN 
                 IF (elo(l,k).LT.(995.d+0)) THEN
                    ilo(l)=ilo(l)+1
                    nlo=nlo+((2*l+1))*mult(i)
                    IF(.NOT.lapw(l).AND.k.EQ.1) GOTO 666
                    IF(k.EQ.nloat) THEN
                          rlo(ilo(l),l)=.TRUE.
                    ENDIF
                    loor(ilo(l),l)=.TRUE.
666                 CONTINUE
                 ENDIF
              ELSE
                 IF (elo(l,k).LT.(995.d+0)) nlov=nlov+((2*l+1))*mult(i)
              ENDIF
           ENDDO
        ENDDO
      ENDDO
      IF(jatom.EQ.nat) GOTO 30                                          
      DO i=jatom+1,nat  
         READ(10) emist                                                    
         READ(10) emist2
#ifdef Parallel
        CALL MPI_BCAST(emist2,(lmax2+1)*nloat,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif         
         DO l = 0,lomax
            DO k=1,nloat
               IF (emist2(l,k).LT.(995.0d+0))  &
                    nlon=nlon+((2*l+1))*mult(i)
            ENDDO
         ENDDO
      ENDDO
  30  CONTINUE     
      IF(myid.EQ.0) THEN
         WRITE(6,13)  JATOM,IATNR(JATOM),(POS(I,LATOM),I=1,3) &
              ,MULT(JATOM)        
         WRITE(21,13) JATOM,IATNR(JATOM),(POS(I,LATOM),I=1,3) &
              ,MULT(JATOM)        
         WRITE(6,1060) ANAME(JATOM)                                    
         DO JROW=1,3                                               
            WRITE(6,1070) ( ROTLOC(JCOL,JROW,JATOM),JCOL=1,3 )         
         ENDDO
         INDEX=LATOM-1                                                 
         DO M=1,MULT(JATOM)                                        
            INDEX=INDEX+1                                               
            WRITE(6,1080) M,( POS(JC,INDEX),JC=1,3 )                    
            DO JR=1,3                                               
               WRITE(6,1090)(ROTIJ(JC,JR,INDEX),JC=1,3),TAUIJ(JR,INDEX)   
            ENDDO
         ENDDO
         WRITE(6,7) ANAME(JATOM)                                           
         WRITE(6,5) E                                                      
         WRITE(6,14)                                                       
      ENDIF

      DO 70 l=0,LMAX2                                                  
         DELE=2.0D-3                                                       
         DELEI=0.25D0/DELE                                                 
         FL=L                                                              
         EI=E(l)/2.0d0                                                         
!                                                                       
!     CALCULATE ENERGY-DERIVATIVE BY FINITE DIFFERENCE                  
!     DELE IS THE UPWARD AND DOWNWARD ENERGY SHIFT IN HARTREES          
!                                                                       
         E1=EI-DELE                                                        
         CALL OUTWIN(REL,VR,R0(JATOM),DX(JATOM),JRI(JATOM),E1,            &
              FL,UVB,DUVB,NODEL,ZZ) 
         CALL RINT13(REL,A,B,A,B,OVLP,JATOM)                               
         TRX=1.0D0/SQRT(OVLP)                                              
         IMAX=JRI(JATOM)                                                   
         DO M=1,IMAX                                                    
            AE(M)=TRX*A(M)                                                    
            BE(M)=TRX*B(M)                                                    
         ENDDO
         UVB=TRX*UVB                                                       
         DUVB=TRX*DUVB                                                     
         E1=EI+DELE                                                        
         CALL OUTWIN(REL,VR,R0(JATOM),DX(JATOM),JRI(JATOM),E1,            &
              FL,UVE,DUVE,NODE,ZZ)
         CALL RINT13(REL,A,B,A,B,OVLP,JATOM)                               
         TRX=1.0d0/SQRT(OVLP)                                                
         UVE=DELEI*(TRX*UVE-UVB)                                           
         DUVE=DELEI*(TRX*DUVE-DUVB)                                        
         IMAX=JRI(JATOM)                                                   
         DO M=1,IMAX                                                    
            AE(M)=DELEI*(TRX*A(M)-AE(M))                                      
            BE(M)=DELEI*(TRX*B(M)-BE(M))                                      
         ENDDO
!                                                                       
!     CALCULATE FUNCTION AT EI                                          
!                                                                       
         CALL OUTWIN(REL,VR(1),R0(JATOM),DX(JATOM),JRI(JATOM),EI,         &
              FL,UV,DUV,NODES,ZZ)
         CALL RINT13(REL,A,B,A,B,OVLP,JATOM)                               
         TRX=1.0d0/SQRT(OVLP)                                                
         P(l)=TRX*UV                                                       
         DP(l)=TRX*DUV                                                     
         IMAX=JRI(JATOM)                                                   
         DO M=1,IMAX                                                    
            A(M)=TRX*A(M)
            B(M)=TRX*B(M)                                                     
         ENDDO
!                                                                       
!     INSURE ORTHOGONALIZATION                                          
!                                                                       
         CALL RINT13(REL,A,B,AE,BE,CROSS,JATOM)                            
         TRY=-CROSS                                                        
         IF(TRY.LT.(-0.05d0).AND.myid.EQ.0) WRITE(6,9) L,TRY,OVLP
         IMAX=JRI(JATOM)                                                   
         DO M=1,IMAX                                                    
            AE(M)=(AE(M)+TRY*A(M))
            BE(M)=(BE(M)+TRY*B(M))                                            
         ENDDO
         IMAX=JRI(JATOM)                                                   
        DO 80 I=1,IMAX                                                    
            RRAD1(I,l)=A(I)                                                   
            RRAD2(I,l)=B(I)                                                   
            RADE1(I,l)=AE(I)                                                  
            RADE2(I,l)=BE(I)                                                  
 80      CONTINUE
         PE(l)=UVE+TRY*P(l)                                                
         DPE(l)=DUVE+TRY*DP(l)                                             
         CALL RINT13(REL,AE,BE,AE,BE,PEI(l),JATOM)                         
         IF(myid.EQ.0) WRITE(6,8) L,P(l),DP(l),PE(l),DPE(l),PEI(l),NODEL,NODES,NODE                         
 70   continue
!                         
! nun fur lo
!
      pr12lo(1:nloat,1:nloat,0:lomax)=0.0d0
      pi12lo(1:nloat,0:lomax)=0.0d0
      pe12lo(1:nloat,0:lomax)=0.0d0
      a1lo(1:nrad,1:nloat,0:lomax)=0.0d0
      b1lo(1:nrad,1:nloat,0:lomax)=0.0d0

      lo_loop: DO l=0,lomax 
         DO jlo=1,ilo(l)
            if (.not.loor(jlo,l)) CYCLE
            DELE=2.0D-3                                 
            DELEI=0.25D0/DELE 
            FL=L
            EI=elo(l,jlo)/2.d0 
!     
!     CALCULATE FUNCTION AT EI                                          
!
            IF(rlo(jlo,l)) THEN
               ei=elo(l,nloat)/2.d0
               kappa=l
!               kappa=-2
               CALL diracout(rel,vr(1),r0(jatom),dx(jatom),jri(jatom),   &
                    ei,fl,kappa,uv,duv,nodes,zz)
!!$               CALL diracout_old(rel,vr(1),r0(jatom),dx(jatom),jri(jatom),   &
!!$                    ei,fl,kappa,uv,duv,nodes,zz)
!               do m=1,jri(jatom)
!                  write(78,'(2f13.7)') a(m),b(m)
!               enddo
               CALL dergl(a,b,r0(jatom),dx(jatom),jri(jatom))
               DO m = 1, jri(jatom)
                  r_m = r0(jatom)*exp(dx(jatom)*(m-1))
                  b(m) = b(m)*r_m / (2.d0*clight+(elo(l,jlo)-2.d0*vr(m)/r_m)/(2.d0*clight))
!                  write(79,'(2f13.7)') r_m,b(m)
!                  b(m)=b(m)*clight
!                   b(m)=0.0d0
               ENDDO
            ELSE
               CALL outwin(rel,vr(1),r0(jatom),dx(jatom),jri(jatom),  &
                    ei,fl,uv,duv,nodes,zz) 
            ENDIF
            CALL rint13(rel,a,b,a,b,ovlp,jatom)
            TRX=1.0d0/SQRT(OVLP)
            plo(jlo,l)=trx*uv 
            dplo(jlo,l)=trx*duv 
            IMAX=JRI(JATOM)
            DO M=1,IMAX 
               a1lo(M,jlo,l)=TRX*A(M) 
               b1lo(M,jlo,l)=TRX*B(M) 
            ENDDO
            IF(l.EQ.1) THEN
               DO m=1,imax
                  r_m = r0(jatom)*EXP(dx(jatom)*(m-1))
                  WRITE(94,'(4e14.6)') r_m,a1lo(m,jlo,l),rrad1(m,l),rade1(m,l)
               ENDDO
            ENDIF
           CALL RINT13(REL,rrad1(1,l),rrad2(1,l),a1lo(1,jlo,l),b1lo(1,jlo,l), &
                 pi12lo(jlo,l),JATOM)
            CALL RINT13(REL,rade1(1,l),rade2(1,l),a1lo(1,jlo,l),b1lo(1,jlo,l), &
                 pe12lo(jlo,l),JATOM)         
            
            IF(myid.EQ.0) WRITE(6,800) L,Plo(jlo,l),DPlo(jlo,l),NODEL,NODES,NODE
         ENDDO

         DO jlo=1,ilo(l)
            CALL ABC(l,jlo,jatom)
         ENDDO

         DO jlo=1,ilo(l)
            IF (.NOT.loor(jlo,l)) CYCLE
            DO jlop=1,ilo(l)
               IF (.NOT.loor(jlop,l)) CYCLE
               CALL RINT13(REL,a1lo(1,jlop,l),b1lo(1,jlop,l), &
                    a1lo(1,jlo,l),b1lo(1,jlo,l), &
                    pr12lo(jlop,jlo,l),JATOM)
            ENDDO
         ENDDO
      ENDDO lo_loop

!
!....overlap integrals for x-dos
            pu1u2=0.d0
            pueu2=0.d0
            pu2u2=0.d0
      do l=0,lxdos
         do lp=0,lxdos
            CALL RINT13(REL,rrad1(1,l),rrad2(1,l),rrad1(1,lp) &
                 ,rrad2(1,lp),pu1u1(l,lp),JATOM)         
            CALL RINT13(REL,rrad1(1,l),rrad2(1,l),rade1(1,lp) &
                 ,rade2(1,lp),pu1ue(l,lp),JATOM)         
            CALL RINT13(REL,rade1(1,l),rade2(1,l),rade1(1,lp) &
                 ,rade2(1,lp),pueue(l,lp),JATOM)         

            do jlo=1,ilo(lp)
!            if(jlo.eq.3) jlo=2   !!! exclude rlo
            if(lp.le.lmax2.and.loor(jlo,lp)) CALL RINT13(REL,rrad1(1,l),rrad2(1,l) &
                 ,a1lo(1,jlo,lp),b1lo(1,jlo,lp),pu1u2(l,lp,jlo),JATOM)         
            if(lp.le.lmax2.and.loor(jlo,lp)) CALL RINT13(REL,rade1(1,l),rade2(1,l) &
                 ,a1lo(1,jlo,lp),b1lo(1,jlo,lp),pueu2(l,lp,jlo),JATOM)         
            do jlop=1,ilo(l)
!            if(jlop.eq.3) jlop=2   !!! exclude rlo
            if((lp.le.lmax2).and.(l.le.lmax2).and.loor(jlo,lp).and.loor(jlop,l)) CALL RINT13(REL,a1lo(1,jlop,l) &
                 ,b1lo(1,jlop,l),a1lo(1,jlo,lp),b1lo(1,jlo,lp),pu2u2(l,JLOP,lp,JLO),JATOM)
            enddo
            enddo
         enddo
      enddo
!
      RETURN                                                            
!                                                                       
    2 FORMAT(5E14.7)                                                    
    3 FORMAT(16I5)                                                      
    4 FORMAT(I4,E16.7)                                                  
    5 FORMAT(10X,' ENERGY PARAMETERS ARE',7F7.2)                        
    6 FORMAT(10X,'E(',I2,')=',F10.4)                                    
    7 FORMAT(/10X,'ATOMIC PARAMETERS FOR ',A10/)                        
   14 FORMAT(/11X,'L',5X,'U(R)',10X,                                     &
             'U''(R)',9X,'DU/DE',8X,'DU''/DE',6X,'NORM-U''')               
    8 FORMAT(10X,I2,5E14.6,5X,3I2)                                      
 800   FORMAT(10X,I2,2E14.6,42X,5X,3I2)                                      
    9 FORMAT(10X,'FOR L=',I2,' CORRECTION=',E14.5,' OVERLAP=',E14.5)    
   11 FORMAT(7F10.5)                                                    
   13 FORMAT(////,':POS',i3.3,':',1x,'AT.NR.',I4,1X,'POSITION =', &
             3F8.5,2X,'MULTIPLICITY =',I3)
 1040 FORMAT(8I2)                                                       
 1050 FORMAT(A10,I5,5X,2F10.9,I5,F5.2)                                  
 1060 FORMAT(//,3X,'NOT EQUIV ATOM ',A10,'  LOCAL ROTATION MATRIX')     
 1070 FORMAT(30X,3F10.5)                                                
 1080 FORMAT(13X,'EQUIV ATOM ',I3,3X,'POSITION: ',3F8.3)                
 1090 FORMAT(30X,3F10.5,5X,F10.5)                                       
 1980 FORMAT(3X)                                                        
 2000 FORMAT(16X,I2//)                                                  
 2022 FORMAT(3X,4E19.12)                                       
 2030 FORMAT(///)
 2031 FORMAT(/)
      END
