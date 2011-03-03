       SUBROUTINE ENERGY (VCOUL0,QEL,LM,LMMAX,ekvol1,iex,etotal)
!                                                                       
!     ENERGY CALCULATES ALL COMPONENTS OF THE   TOTAL ENERGY            
!     EXCEPT THE SUM OF EIGENVALUES.                                    
!     THIS ALGORITHM FOLLOWS THE PAPER OF M.WEINERT,E.WIMMER,           
!     PHYS.REV.B,26,4571-4578 (1982)                                    
!                                                                       
!     FIVE ENERGY VARIABLES ARE DEFINED  :                              
!     1. ECONST                                                         
!        ECONST DESCRIBES THE CONSTANT CONTRIBUTION OF EACH SPHERE,     
!        CONSISTING OF NUCLEAR CHARGE ZZ, ELECTRONIC CHARGE QEL AND     
!        SPHERICAL AVERAGE OF POTENTIAL VCOUL0                          
!     2. EKVOL                                                          
!        EKVOL IS THE INTEGRAL (PW-CHARGE RHOK * PW-POTENTIAL VKSUM)    
!        OVER CELL VOLUME OMEGA                                         
!     3. EKSPHER                                                        
!        DESCRIBES THE VALUE OF EKVOL IN SPHERES AND COMPENSATES IT     
!        IN THE FINAL SUM (TO ACHIEVE INTEGRATION ONLY OVER PW.         
!        BOTH ARGUMENTS ARE SPLIT UP INTO (CUBIC) HARMONICS  BY THE     
!        RAYLEIGH-EXPANSION OF EXP(IKR),THE RADIAL COEFFICIENTS ARE     
!        CALLED RADRHOK AND RADPOTK                                     
!     4. ESPHERE                                                        
!        ESPHERE IS THE INTEGRAL (     SPHERE-CHARGE CLM * SPHERE -     
!        POTENTIAL VRADSUM + TERM(ZZ) ) OVER RADIUS                     
!     5. ETOTAL                                                         
!        ETOTAL IS THE SUM OF TERMS 1-4                                 
!                                                                       
!     POTENTIAL (VKSUM,VRADSUM) MEANS SUM OF  ( 0.5*COULPOT + XCPOT     
!     -EPSILON ) IN THE ABOVE GIVEN TERMS.                              
!     ALL INPUT-DATAS ARE GIVEN IN RYDBERG,   ALL ENERGY COMPONENTS     
!     ARE CALCULATED IN HARTREES                                        
!                                                                       
!
!     Addition by LDM, Jan 2007, extra nuclear energy term for Efield
!                                                                       
!                                                                       
!                                                                       
      use struct
      use densit
      use work     !, only : v,r,vcoul2,value
      use vresp1
      use parallel
      use efeld
      IMPLICIT NONE 
!
      INCLUDE 'param.inc'
#ifdef Parallel
      INCLUDE 'mpif.h'
#endif
!
      REAL*8   ::  c_kub(0:10,0:10)
      COMMON /norm_kub/ c_kub

      LOGICAL  ::  HEXA,CUBIC                                       
      INTEGER  ::  bufnelec,NELEC,INDEX,JATOM,MU,J,L,M,LMX,ISPIN
      INTEGER  ::  IEX,LMMAX(NAT),LM(2,NCOM+3,NAT), jat
!      REAL*8   ::  exc(1:nrad,1:ncom,1:nat,1:2)
!      REAL*8   ::  vclmsp(1:nrad,1:ncom,1:nat,1:2)
!      REAL*8   ::  vxc(1:nrad,1:ncom,1:nat,1:2)
      REAL*8   ::  RHO(NRAD)
      REAL*8   ::  ECONST(NAT),ESP(NCOM,NAT)            
      REAL*8   ::  EKINSP(NCOM,NAT),RHO1(NRAD)         
      REAL*8   ::  C(NCOM),BES(0:20)                
      REAL*8   ::  VCOUL0(1:nat),QEL(NAT),CINV(NCOM)
      REAL*8   ::  PI,SQ4PII,EKVOL,ETOTAL,EKVOL1,EINT,FIELDE
      REAL*8   ::  SUM,SUMKIN,bufetotal
      COMPLEX*16   ::  RRHOK,RPOTK,IMAG,IMAG1,CFAC
      COMPLEX*16   ::  RKINK,imagl,temp                

      INTEGER :: imode, ifourier
      REAL*8  :: z, tmp1, backget
      
      PI    =ACOS(-1.D0)                                                
      IMAG  =(0.D0,1.D0)                                                
      sq4pii=1.d0/sqrt(4.d0*pi)
      EKVOL =0.d0                                                         
      ETOTAL=0.d0                                                         
      bufetotal=0.d0                                                         
      bufNELEC =0.d0                                                         
      NELEC =0.d0                                                         
!                                                                       
!     START LOOP OVER ALL ATOMS                                         
!     ---------------------------------------------------               
      DO 20 JATOM=NSTART(myid),NSTOP(myid)
      DO 20 MU=1,MULT(JATOM)                                            
         HEXA  =.FALSE.                                                 
         CUBIC =.FALSE.                                                 
         IF(IATNR(JATOM).LT.0) HEXA =.TRUE.                             
         IF(IATNR(JATOM).GT.0) CUBIC=.TRUE.                             
         DO 15 J =1,JRI(JATOM)                                          
 15      R(J)=R0(JATOM)*EXP( DX(JATOM)*(J-1) )                        
!                                                                       
!........CALCULATE ECONST                                               
         ECONST(JATOM)=-0.5D0*ZZ(JATOM)/RMT(JATOM)*( ZZ(JATOM)+          &
                       RMT(JATOM)*VCOUL0(JATOM)/2.D0-QEL(JATOM) )*       &
                       MULT(JATOM)                                      

!........CALCULATE ESPHERE                                              
         DO 40 LMX=1,LMMAX(JATOM)                                       
            L=IABS( LM(1,LMX,JATOM) )                                   
            M=      LM(2,LMX,JATOM)                                     
         ESP(LMX,JATOM)=0.0                                             
         EKINSP(LMX,JATOM)=0.0                                          
         DO 40 ISPIN=1,2  
            IF(LMX.EQ.1) THEN                                                 
            DO 45 J=1,JRI(JATOM)                                        
               RHO(J)=CLMSP(J,LMX,JATOM,ISPIN)*sq4pii           
               VALUE(J)=-RHO(J)*(VXC(J,LMX,JATOM,ISPIN)+ &
               VCOUL2(J,JATOM)-EXC(J,LMX,JATOM,ISPIN))*0.5d0            
               RHO1(J)=-RHO(J)*(V(J,LMX,JATOM)-4.D0* &
               EXC(J,LMX,JATOM,ISPIN)+2.D0*VXC(J,LMX,JATOM,ISPIN)) &
               *0.5d0                             
 45         CONTINUE                                                    
!
      if(iex.eq.30.or.iex.eq.31.or.iex.eq.32.or.iex.eq.33.or. &
      iex.eq.34.or.iex.eq.41.or.iex.eq.42) then
      write(6,*) ' Exchange energy only in energy'
            DO 345 J=1,JRI(JATOM)                                        
                VALUE(J)=-RHO(J)*(-EXC(J,LMX,JATOM,ISPIN))*0.5d0            
 345         CONTINUE                                                    
      else if(iex.eq.40) then
            DO  J=1,JRI(JATOM)                                        
             VALUE(J)=(-RHO(J)*(VXC(J,LMX,JATOM,ISPIN)+V(J,lmx,JATOM))+ &
                      EXC(J,LMX,JATOM,ISPIN)*0.25d0*r(j)*r(j)- &
                      vclmsp(j,lmx,jatom,ispin))*0.5d0
            ENDDO                                                    
           write(6,*) (vxc(1,1,jatom,1)+V(1,1,JATOM))*r(1)
      endif  
!
      ELSE IF(CUBIC.AND.(l.ne.9.and.m.eq.4)) then
               goto 40
      else if (cubic.and.(m.EQ.6.OR.m.EQ.8.OR.m.EQ.10)) THEN 
               GOTO 40
            ELSE IF(CUBIC.AND.(l.EQ.4.OR.l.EQ.6 &
               .OR.l.EQ.7.OR.l.EQ.8.OR.l.EQ.9.OR.l.EQ.10)) THEN
            DO 145 J=1,JRI(JATOM)
                  RHO(J)= &
                  CLMSP(J,LMX,JATOM,ISPIN)*C_kub(l,m)+CLMSP(J,LMX+1, &
                  JATOM,ISPIN)*C_kub(l,m+4)
         IF(L.EQ.8.OR.l.EQ.10)  &
             rho(j)=rho(j)+CLMSP(J,LMX+2,JATOM,ISPIN)*C_kub(l,m+8)
                  rho1(j)= &
                  (VXC(J,LMX,JATOM,ISPIN)+V(J,LMX,JATOM)*0.5d0- &
                    EXC(J,LMX,JATOM,ISPIN))*C_kub(l,m) + &
                  (VXC(J,LMX+1,JATOM,ISPIN)+V(J,LMX+1,JATOM)*0.5d0- &
                    EXC(J,LMX+1,JATOM,ISPIN))*C_kub(l,m+4)
         IF(L.EQ.8.OR.l.EQ.10) rho1(j)=rho1(j)+ &
                  (VXC(J,LMX+2,JATOM,ISPIN)+V(J,LMX+2,JATOM)*0.5d0- &
                    EXC(J,LMX+2,JATOM,ISPIN))*C_kub(l,m+8)  
                  VALUE(J)=-RHO(J)*rho1(j)*0.5d0
      RHO1(J)=-RHO(J)*(V(J,LMX,JATOM)-4.D0*EXC(J,LMX,JATOM,ISPIN) &
        +2.D0*VXC(J,LMX,JATOM,ISPIN))*0.5d0
145         CONTINUE
!
      if(iex.eq.30.or.iex.eq.31.or.iex.eq.32.or.iex.eq.33.or. &
      iex.eq.34.or.iex.eq.41.or.iex.eq.42) then
      write(6,*) ' Exchange energy only in energy'
            DO 445 J=1,JRI(JATOM)                                        
                  RHO(J)= &
                  CLMSP(J,LMX,JATOM,ISPIN)*C_kub(l,m)+CLMSP(J,LMX+1, &
                  JATOM,ISPIN)*C_kub(l,m+4)
               IF(L.EQ.8.OR.l.EQ.10)  &
                rho(j)=rho(j)+CLMSP(J,LMX+2,JATOM,ISPIN)*C_kub(l,m+8)
                  rho1(j)= &
                  (-EXC(J,LMX,JATOM,ISPIN))*C_kub(l,m) + &
                  (-EXC(J,LMX+1,JATOM,ISPIN))*C_kub(l,m+4)
         IF(L.EQ.8.OR.l.EQ.10) rho1(j)=rho1(j)+ &
                  (-EXC(J,LMX+2,JATOM,ISPIN))*C_kub(l,m+8)
                  VALUE(J)=-RHO(J)*rho1(j)*0.5d0
 445         CONTINUE                                                    
      else if(iex.eq.40) then
            DO  j=1,jri(jatom)
                  rho(j)= &
                  clmsp(j,lmx,jatom,ispin)*c_kub(l,m)+clmsp(j,lmx+1, &
                  jatom,ispin)*c_kub(l,m+4)
         IF(l.EQ.8.OR.l.EQ.10) rho(j)=rho(j)+ &
                   clmsp(j,lmx+2,jatom,ispin)*c_kub(l,m+8)
             value(j)=(-rho(j)*((vxc(j,lmx,jatom,ispin)+v(j,lmx,jatom)) &
                 *c_kub(l,m)+(vxc(j,lmx+1,jatom,ispin)+v(j,lmx+1,jatom)) &
                 *c_kub(l,m+4))+ &
               (exc(j,lmx,jatom,ispin)*0.25d0*r(j)*r(j)- &
                vclmsp(j,lmx,jatom,ispin)) &
                 *c_kub(l,m)+(exc(j,lmx+1,jatom,ispin)*0.25d0*r(j)*r(j)- &
                  vclmsp(j,lmx+1,jatom,ispin))*c_kub(l,m+4))*0.5d0
         IF(l.EQ.8.OR.l.EQ.10) value(j)=value(j) &
                        -0.5d0*rho(j)*(c_kub(l,m+8)* &
                          (vxc(j,lmx+2,jatom,ispin)+v(j,lmx+2,jatom))+ &
                c_kub(l,m+8)*(exc(j,lmx+2,jatom,ispin)*0.25d0*r(j)*r(j)- &
                  vclmsp(j,lmx+2,jatom,ispin)))
            enddo                                                    
      endif  
!
            ELSE                          
            DO 245 J=1,JRI(JATOM)
               RHO(J)=CLMSP(J,LMX,JATOM,ISPIN)
         VALUE(J)=-RHO(J)*(VXC(J,LMX,JATOM,ISPIN)+V(J,LMX,JATOM)*0.5d0- &
                  EXC(J,LMX,JATOM,ISPIN))*0.5d0
      RHO1(J)=-RHO(J)*(V(J,LMX,JATOM)-4.D0*EXC(J,LMX,JATOM,ISPIN) &
        +2.D0*VXC(J,LMX,JATOM,ISPIN))*0.5d0
245         CONTINUE
!
      if(iex.eq.30.or.iex.eq.31.or.iex.eq.32.or.iex.eq.33.or. &
      iex.eq.34.or.iex.eq.41.or.iex.eq.42) then
      write(6,*) ' Exchange energy only in energy'
            DO 545 J=1,JRI(JATOM)                                        
                VALUE(J)=-RHO(J)*(-EXC(J,LMX,JATOM,ISPIN))*0.5d0            
 545         CONTINUE                                                    
      else if(iex.eq.40) then
            DO  J=1,JRI(JATOM)                                        
             VALUE(J)=(-RHO(J)*(VXC(J,LMX,JATOM,ISPIN)+V(J,lmx,JATOM))+ &
                      EXC(J,LMX,JATOM,ISPIN)*0.25d0*r(j)*r(j)- &
                      vclmsp(j,lmx,jatom,ispin))*0.5d0
!             VALUE(J)=(RHO(J)*(VXC(J,LMX,JATOM,ISPIN)+V(J,lmx,JATOM))
!     &                )*0.5d0
            enddo                                                    
      endif  
!
           END IF
!           INTEGRATE VALUE OVER RADIUS                                 
            CALL CHARGE (R,DX(JATOM),VALUE,1,JRI(JATOM),EINT)           
            ESP(LMX,JATOM)=ESP(LMX,JATOM) + EINT                        
            CALL CHARGE (R,DX(JATOM),RHO1,1,JRI(JATOM),EINT)            
            EKINSP(LMX,JATOM)=EKINSP(LMX,JATOM) + EINT                  
 40      CONTINUE                                                       
!                                                                       
!     END LOOP OVER ATOMS                                               
 20   CONTINUE                                                          
!                                                                       
!                                                                       
!.....WRITE OUT RESULTS TO TAPE6 AND TAPE21=SCFDATA                     
      WRITE(6,1000)                                                     
!     WRITE OUT DENSITY-POTENTIAL INTEGRALS                             
      INDEX=0                                                           
      DO jatom=1,nstart(myid)-1
         index=index+mult(jatom)
      ENDDO
      DO 50 jatom=nstart(myid),nstop(myid)
      DO 50 MU=1,MULT(JATOM)                                            
         INDEX=INDEX+1                                                  
         HEXA  =.FALSE.                                                 
         CUBIC =.FALSE.                                                 
         IF(IATNR(JATOM).LT.0) HEXA =.TRUE.                             
         IF(IATNR(JATOM).GT.0) CUBIC=.TRUE.                             
         SUM=0.0                                                        
         SUMKIN=0.
         IF(jatom.EQ.nstart(myid).AND.mu.EQ.1) WRITE(6,1040)
         WRITE(6,1050)  INDEX,ANAME(JATOM)                              
         DO 55 LMX=1,LMMAX(JATOM)                                       
            L=IABS( LM(1,LMX,JATOM) )                                   
            M=      LM(2,LMX,JATOM)                                     
            IF(CUBIC) THEN                                              
               IF(m.EQ.4.OR.m.EQ.6.OR.m.EQ.8.OR.m.EQ.10) GOTO 55
            ENDIF                                                       
            WRITE(6,1060) L,ESP(LMX,JATOM),EKINSP(LMX,JATOM)           
            SUM=SUM + ESP(LMX,JATOM)
 55      CONTINUE                                                       
         bufETOTAL=bufETOTAL+SUM                                              
         WRITE(6,1070)  SUM                                      
 50   CONTINUE                                                          
!                                                                       
!     WRITE OUT EKVOL, CONSTANT TERMS                                   
      WRITE(6,1080)  -EKVOL1*0.5d0
!ccccc
      ETOTAL=ETOTAL+EKVOL1*(-0.5d0)                                           
      if(iex.eq.30.or.iex.eq.31.or.iex.eq.32.or.iex.eq.33.or. &
      iex.eq.34.or.iex.eq.40.or.iex.eq.41.or.iex.eq.42) goto 61
      WRITE(6,1090)                                                     
      DO 60 JATOM=NSTART(myid),NSTOP(myid)
         bufNELEC=bufNELEC + ZZ(JATOM)*MULT(JATOM)                            
         WRITE(6,1100)  JATOM,ANAME(JATOM),ECONST(JATOM)                
         bufETOTAL=bufETOTAL+ECONST(JATOM)                                    
!        Include here the core term for an Efield
         if(abs(iefeld) .gt. 0)then
                INDEX=1
                DO jat=1,jatom-1
                   index=index+mult(jat)
                ENDDO
                z=pos(3,index)
!               LDM general case, core energy in Hartrees
                imode=abs(iefeld)/1000
                ifourier=1
                tmp1= -backget(z,imode,wefeld,ifourier)*0.5D0*refeld &
                        *ZZ(JATOM)*MULT(JATOM)
                write(6,1101)tmp1
                bufETOTAL=bufETOTAL+tmp1
         endif
 60   CONTINUE                                                          
!                                                                       
!     WRITE OUT MUFFIN-TIN ZERO * NUMBER OF ELECTRONS,  ETOTAL          

 61   CONTINUE

#ifdef Parallel
         call MPI_ALLREDUCE(bufnelec,nelec,1,MPI_INTEGER,            &
              MPI_SUM,comm_coul,ierr)
         call MPI_ALLREDUCE(bufETOTAL,etotal,1,MPI_DOUBLE_PRECISION, &
              MPI_SUM,comm_coul,ierr)
#else
               etotal=bufetotal
               nelec=bufnelec
#endif
      ETOTAL=ETOTAL+EKVOL1*(-0.5d0)
       WRITE(6,1120)  ETOTAL
      ETOTAL=ETOTAL*2.D0                                                
      WRITE(6,1121)  ETOTAL
      RETURN                                                            
!                                                                       
!                                                                       
 1000 FORMAT(1H1,/,30X,30(1H-),/,33X,'T O T A L   E N E R G Y',          &
             /   ,30X,30(1H-),/,/)                                       
 1040 FORMAT(3X,'1. DENSITY-POTENTIAL INTEGRALS',5X,'(H)',/)            
 1050 FORMAT(9X,'SPHERE=',I2,3X,A10,/)                                  
 1060 FORMAT(12X,'L=',I2,'   :',2F13.6,5X,'SPHERE CHARGE',/,             &
             20X,2F13.6,5X,'PLANE WAVE CHARGE  EXPANDED',/ )            
 1070 FORMAT(12X,'TOTAL=  ',F13.6,/,/)                                  
 1080 FORMAT(9X,'PLANE WAVE  CONTRIBUTION',/,/,12X,'TOTAL=  ',F13.6,/,/) 
 1090 FORMAT(3X,'2. CONSTANT TERMS',5X,'(H)',/)                         
 1100 FORMAT(9X,'SPHERE',I3,3X,A10,'  =',F13.6)                        
 1101 FORMAT(21X,'EField      =',F13.6)
 1120 FORMAT(3X,80(1H-),/,/,3X,'TOTAL VALUE =',F14.6,5X,'(H)')          
 1121 FORMAT(':DEN  :',' DENSITY INTEGRALS =',3x,2F20.6,3X,'(Ry)')           
!                                                                       
      END                                                               
