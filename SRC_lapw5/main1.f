      SUBROUTINE MAIN1                                                   
      use atpos
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*4 LATTIC, IUNIT,SWITCH,addsub,cform,cnorm,deb        
      CHARACTER*10 TITEL                                                
      character*1 noortho,noorth1
      INCLUDE 'param.inc'
      COMMON /SYM2/ IZ(3,3,NSYM),TAU(3,NSYM),IORD                       
      real*8,allocatable ::  RM(:,:),CLM(:,:,:),ROTLOC(:,:,:) 
      real*8,allocatable :: RNOT(:),RMT(:),DX(:),z(:)
      integer,allocatable :: JRI(:),LM(:,:,:),LMMAX(:),mult(:)  
      DIMENSION BR1(3,3),VEC(3,3),A(3)                 
      dimension TITEL(8),VX(3),VY(3)      
      DIMENSION XI(3),BR2(3,3),br3(3,3),vt1(3)                      
      DIMENSION V(3),VT(3)
      integer,allocatable :: IATNR(:),IOP(:)    
      real*8,allocatable :: POS(:,:)    
      logical ortho
      DATA BR1/9*0.0/,BR2/9*0.0/                                        
!                                                                       
!                      
      noortho='O'                                               
      READ(8,102) TITEL                                                 
      READ(8,103) LATTIC,NAT,cform
      nato=nat              
      ndif=48*nato                              
      allocate ( RM(NRAD,NATO),CLM(NRAD,NCOM,NATO),ROTLOC(3,3,NATO)) 
      allocate ( RNOT(NATO),RMT(NATO),DX(NATO),z(nato))
      allocate ( JRI(NATO),LM(2,NCOM,NATO),LMMAX(NATO),mult(nato) ) 
      allocate ( IATNR(NDIF),IOP(NDIF))    
      allocate ( POS(3,NDIF))    
      READ(8,100) A,alpha,beta,gamma
      if(gamma.eq.0.d0) gamma=90.d0
      gamma1=gamma
      alpha1=alpha
      beta1=beta
      alpha=alpha/180.d0*acos(-1.d0)                                          
      beta=beta/180.d0*acos(-1.d0)                                          
      gamma=gamma/180.d0*acos(-1.d0)                                          
      WRITE(6,119) TITEL                                                
      WRITE(6,104) LATTIC,NAT                                           
      WRITE(6,105) A                                                    
!      if(cform.eq.'NEW ') then
        assign 2021 to iform1
!      else
!        assign 2020 to iform1
!      end if                                    
 2020 FORMAT(3X,5E13.7)                                                 
 2021 FORMAT(3X,4E19.12)                                                 
!                                                                       
!....DEFINE DIRECT space BRAVAIS MATRIX  BR2(i,*))                            
!....BR1    DIRECT space (non-primitive cell)     
      IF(LATTIC(1:1).EQ.'S'.OR.LATTIC(1:1).EQ.'P') THEN
      cosg1=(cos(gamma)-cos(alpha)*cos(beta))/sin(alpha)/sin(beta)
      gamma0=acos(cosg1)
      BR2(1,1)=A(1)*1.0d0*sin(gamma0)*sin(beta)
      BR2(1,2)=A(1)*1.0d0*cos(gamma0)*sin(beta)
      BR2(1,3)=A(1)*1.0d0*cos(beta)    
      BR2(2,1)=0.0d0              
      BR2(2,2)=A(2)*1.0d0*sin(alpha)
      BR2(2,3)=A(2)*1.0d0*cos(alpha)
      BR2(3,1)=0.0d0              
      BR2(3,2)=0.0d0              
      BR2(3,3)=A(3)*1.0d0              
      BR1(1,1)=A(1)*1.0d0*sin(gamma0)*sin(beta)
      BR1(1,2)=A(1)*1.0d0*cos(gamma0)*sin(beta)
      BR1(1,3)=A(1)*1.0d0*cos(beta)    
      BR1(2,1)=0.0d0              
      BR1(2,2)=A(2)*1.0d0*sin(alpha)
      BR1(2,3)=A(2)*1.0d0*cos(alpha)
      BR1(3,1)=0.0d0              
      BR1(3,2)=0.0d0              
      BR1(3,3)=A(3)*1.0d0              
!  old version with monoclinic lattic only
!      BR2(1,1)=A(1)*sin(gamma)
!      br2(1,2)=a(1)*cos(gamma)                                               
!      BR2(2,2)=A(2)                                                     
!      BR2(3,3)=A(3)                                                     
!      BR1(1,1)=1.*A(1)*sin(gamma)
!      br1(1,2)=1.*a(1)*cos(gamma)                                             
!      BR1(2,2)=1.*A(2)                                                  
!      BR1(3,3)=1.*A(3) 
      ortho=.true.
        if(gamma1.ne.90.d0) ortho=.false.
        if(alpha1.ne.90.d0) ortho=.false.
        if(beta1.ne.90.d0) ortho=.false.
      ELSE IF(LATTIC(1:1).EQ.'F') THEN                                  
      BR2(1,1)=0.5*A(1)                                                 
      BR2(2,1)=0.5*A(1)                                                 
      BR2(2,2)=0.5*A(2)                                                 
      BR2(3,2)=0.5*A(2)                                                 
      BR2(1,3)=0.5*A(3)                                                 
      BR2(3,3)=0.5*A(3)                                                 
      BR1(1,1)=A(1)                                                     
      BR1(2,2)=A(2)                                                     
      BR1(3,3)=A(3)                                                     
      ortho=.true.                                                 
      ELSE IF(LATTIC(1:1).EQ.'B') THEN                                  
      BR2(1,1)=-0.5*A(1)                                                
      BR2(2,1)=0.5*A(1)                                                 
      BR2(3,1)=0.5*A(1)                                                 
      BR2(1,2)=0.5*A(2)                                                 
      BR2(2,2)=-0.5*A(2)                                                
      BR2(3,2)=0.5*A(2)                                                 
      BR2(1,3)=0.5*A(3)                                                 
      BR2(2,3)=0.5*A(3)                                                 
      BR2(3,3)=-0.5*A(3)                                                
      BR1(1,1)=A(1)                                                     
      BR1(2,2)=A(2)                                                     
      BR1(3,3)=A(3)                                                     
      ortho=.true.                                                 
      ELSE IF(LATTIC(1:1).EQ.'H') THEN                                  
      BR1(1,1)=SQRT(3.)/2.*A(1)                                         
      BR1(1,2)=-0.5*A(2)                                                
      BR1(2,2)=1.*A(2)                                                  
      BR1(3,3)=1.*A(3)                                                  
      BR2(1,1)=SQRT(3.)/2.*A(1)                                         
      BR2(1,2)=-0.5*A(2)                                                
      BR2(2,2)=1.*A(2)                                                  
      BR2(3,3)=1.*A(3)                                                  
      ortho=.false.                                                 
      ELSE IF(LATTIC(1:1).EQ.'R') THEN                                  
      BR1(1,1)=1./SQRT(3.)/2.*A(1)                                         
      BR1(1,2)=-0.5*A(2)                                                
      BR1(1,3)=1./3.*A(3)                                                
      BR1(2,1)=1./SQRT(3.)/2.*A(1)                                         
      BR1(2,2)=0.5*A(2)                                                
      BR1(2,3)=1./3.*A(3)                                                
      BR1(3,1)=-1./SQRT(3.)*A(1)                                         
      BR1(3,2)=0.                                         
      BR1(3,3)=1./3.*A(3)                                                  
      BR2(1,1)=1./SQRT(3.)/2.*A(1)                                         
      BR2(1,2)=-0.5*A(2)                                                
      BR2(1,3)=1./3.*A(3)                                                
      BR2(2,1)=1./SQRT(3.)/2.*A(1)                                         
      BR2(2,2)=0.5*A(2)                                                
      BR2(2,3)=1./3.*A(3)                                                
      BR2(3,1)=-1./SQRT(3.)*A(1)                                         
      BR2(3,2)=0.                                         
      BR2(3,3)=1./3.*A(3)                                                  
      ortho=.false.                                                 
      ELSE IF(LATTIC(1:3).EQ.'CXY') THEN                                  
      BR2(1,1)=0.5*A(1)                                                
      BR2(2,1)=0.5*A(1)                                                 
      BR2(3,1)=0.                                                 
      BR2(1,2)=0.5*A(2)                                                 
      BR2(2,2)=-0.5*A(2)                                                
      BR2(3,2)=0.                                                 
      BR2(1,3)=0.                                                 
      BR2(2,3)=0.                                                 
      BR2(3,3)=A(3)                                                
      BR1(1,1)=A(1)                                                     
      BR1(2,2)=A(2)                                                     
      BR1(3,3)=A(3)                                                     
      ortho=.true.                                                 
      ELSE IF(LATTIC(1:3).EQ.'CYZ') THEN                                  
      BR2(1,1)=A(1)                                                
      BR2(2,1)=0.                                                 
      BR2(3,1)=0.                                                 
      BR2(1,2)=0.                                                 
      BR2(2,2)=-0.5*A(2)                                                
      BR2(3,2)=0.5*A(2)                                                 
      BR2(1,3)=0.                                                 
      BR2(2,3)=0.5*A(3)                                                 
      BR2(3,3)=0.5*A(3)                                                
      BR1(1,1)=A(1)                                                     
      BR1(2,2)=A(2)                                                     
      BR1(3,3)=A(3)                                                     
      ortho=.true.                                                 
!      ELSE IF(LATTIC(1:3).EQ.'CXZ'.and.gamma1.ne.90.d0) THEN                  
      ELSE IF(LATTIC(1:3).EQ.'CXZ') THEN                  
      BR2(1,1)=0.5d0*A(1)*sin(gamma)                                          
      BR2(1,2)=0.5d0*A(1)*cos(gamma)                                         
      BR2(1,3)=-0.5d0*A(3)                                                 
      BR2(2,2)=A(2)                                                
      BR2(3,1)=0.5d0*A(1)*sin(gamma)                                          
      BR2(3,2)=0.5d0*A(1)*cos(gamma)                                         
      BR2(3,3)=0.5d0*A(3)                                                
      BR1(1,1)=1.*A(1)*sin(gamma)
      br1(1,2)=1.*a(1)*cos(gamma)                                             
      BR1(2,2)=1.*A(2)                                                  
      BR1(3,3)=1.*A(3) 
      ortho=.false.                                                 
      ELSE                                                              
      STOP 'LATTIC NOT DEFINED'                                         
      END IF                                                            
      call gbass(br1,br3)
!                                                                       
!.... READ NUMBER OF DIFFERENT ATOMS(MAX 17)                            
!.... READ AT.NR (EQ.STRUCT) AND POSITION                               
!                                                                       
      WRITE(6,106)((BR1(I,J),I=1,3),J=1,3)                              
      WRITE(6,106)((BR3(I,J),I=1,3),J=1,3)                              
      WRITE(6,108)                                                      
      INDEX=0                                                           
      DO 2 JATOM=1,NAT                                                  
      INDEX=INDEX+1                                                     
      READ(8,1012) IATNR(INDEX),POS(1,INDEX),POS(2,INDEX),POS(3,INDEX),  &
       MULT(JATOM)                                                      
!                                                                       
!    POS IATNR: MEANS CUB SYM                                           
!    NEG IATNR: MEANS NO CUB SYM                                        
      DO 51 MU=1,MULT(JATOM)-1                                          
      INDEX=INDEX+1                                                     
      READ(8,1013) IATNR(INDEX),POS(1,INDEX),POS(2,INDEX),POS(3,INDEX)  
 51   CONTINUE                                                          
!.... READ(RADIAL PARAMETERS, Z, SYMMETRY OPERATIONS                    
!                                                                       
      READ(8,113) ANAME,JRI(JATOM),RNOT(JATOM),RMTT,Z(JATOM)            
      DX(JATOM)=LOG(RMTT/RNOT(JATOM)) / (JRI(JATOM)-1)                 
      READ(8,1051) ((ROTLOC(I1,J1,JATOM),I1=1,3),J1=1,3)                
 2    CONTINUE                                                          
      NDAT=INDEX                                                        
      READ(8,114) IORD                                                  
      DO 11 I=1,IORD                                                    
   11 READ(8,115) ((IZ(I1,I2,I),I1=1,3),TAU(I2,I),I2=1,3)               
      CALL ROTDEF(NAT,MULT,IOP,POS,lattic)                                     
      DO 53 INDEX=1,NDAT                                                
      AIX=POS(1,INDEX)                                                  
      AIY=POS(2,INDEX)                                                  
      AIZ=POS(3,INDEX)                                                  
      DO 3 J=1,3                                                        
    3 POS(J,INDEX)=BR1(1,J)*AIX+BR1(2,J)*AIY+BR1(3,J)*AIZ               
 53   WRITE(6,107) IATNR(INDEX),IOP(INDEX),(POS(J,INDEX),J=1,3)         
!                                                                       
!.... READ VECTORS OF ORIGIN AND ENDPOINTS OF THE PLOT                  
!                                                                       
      WRITE(6,109)                                                      
      DO 4 I=1,3                                                        
      READ(5,*) IX,IY,IIZ,IDV                                         
      DO 5 J=1,3                                                        
    5 VEC(J,I)=BR1(1,J)*IX/IDV+BR1(2,J)*IY/IDV+BR1(3,J)*IIZ/IDV         
    4 WRITE(6,100) (VEC(J,I),J=1,3)                                     
!                                                                       
!..... READ NSHELL-PARAMETERS IN X,Y AND Z- DIRECTION                   
!                                                                       
      READ(5,*) NSA,NSB,NSC                                           
      WRITE(6,120) NSA,NSB,NSC                                          
  120 FORMAT('0NSHELL PARAMETERS IN X, Y AND Z-DIRECTION:',3I5)         
!                                                                       
!.... READ NUMBER OF PLOTTING AND PRINT - POINTS IN X AND Y -DIR        
!.... READ SWITCH: RHO ,DIFF,OVER, ....N; IUNIT: ANG OR ATU             
!      SWITCH1:  NUM,CON,3D                                             
      READ(5,*) NPX,NPY
!      if(npx.gt.nxdim) stop 'decrease NPX or increase PARAMETER NXDIM'
!      if(npy.gt.nydim) stop 'decrease NPY or increase PARAMETER NYDIM'
        NPXO=1
        NPYO=1                                     
      WRITE(6,110) NPX,NPY,NPXO,NPYO                                    
      READ(5,112) SWITCH,addsub,IUNIT,cnorm,deb                           
      read(5,'(a1)',end=444,err=444) noorth1
      noortho=noorth1
 444  continue
      WRITE(6,111) SWITCH,IUNIT 
!                                                                       
!.... GENERATE DIRECTION OF BASISLINES                                  
!                                                                       
      CALL VDIF(VEC,VEC(1,2),VX)                                        
      CALL VDIF(VEC,VEC(1,3),VY)                                        
      TEST1=VX(1)*VY(1)+VX(2)*VY(2)+VX(3)*VY(3)
      test2=test1/sqrt(VX(1)**2+VX(2)**2+VX(3)**2)/sqrt(Vy(1)**2+Vy(2)**2+Vy(3)**2)                           
      test2=acos(test2)*180.d0/acos(-1.d0)
      IF(abs(TEST1).GT.0.001) then
       write(6,*) ' directions not orthogonal '
       write(6,*) ' vector 1:',vx
       write(6,*) ' vector 2:',vy
       write(6,*) ' cos(alpha), alpha:',test1,test2
      if(npy.gt.1.and.noortho.eq.'O') STOP 'DIR WRONG'            
      end if
! .... TEST ORTHOGONALITY OF BASIS LINES                                
      pi=acos(-1.d0)
      sqfp=sqrt(4.d0*pi)                                        
      FAK=1.
      fadd=0.d0             
      if(addsub.eq.'ADD ') fadd=1.d0
      if(addsub.eq.'SUB ') fadd=-1.d0
      if(addsub.eq.' ADD') fadd=1.d0
      if(addsub.eq.' SUB') fadd=-1.d0
                                  
      IF(IUNIT.EQ.'ANG ') FAK=6.7551d0                                    
      R=VNORM(VX)                                                       
      RELX=R                                                            
      R=VNORM(VY)                                                       
      RELY=R                                                            
      CMIN=1000.                                                        
      CMAX=0.0                                                          
      WRITE(10) NPX,NPY,RELX,RELY                                       
      WRITE(21,'(2i5,2f10.5)') NPX,NPY,RELX,RELY
!                                                                       
!.... GENERATE ORIGINS OF NEXT UNIT CELLS                               
!                                                                       
      CALL GENER(NSA,NSB,NSC,BR2)                                       
!                                                                       
      IF(SWITCH.NE.'OVER'.AND.SWITCH.NE.'DIFF') GOTO 12
      write(6,*) 'YOU HAVE SELECTED SWITCH:',switch                 
!                                                                       
!.... READ SIGMA OBTAINED BY "ATOM" AND INTERPOLATE ON UNIFORM MESH     
!                                                                       
      CALL HSLGRS(Z,NAT)                                                     
   12 CONTINUE                                                          
      IF(SWITCH.EQ.'OVER') GOTO 16                                      
!                                                                       
!.... READ LMMAX OF ALL ATOMS                                           
!.... READ CLMs, GENERATE R-MESH                                       
!                                                                       
      READ(9,2032)                                                      
      READ(11,2032,END=38,IOSTAT=ICLM)              
      if(fadd.gt.0.d0) then
          write(6,*) 'YOU ARE ADDING THE DENSITIES FROM UNITs 9 and 11'
      else if(fadd.lt.0.d0) then
          write(6,*) 'YOU ARE SUBTRACTING THE DENSITY OF UNIT 11 from 9'
      else
          write(6,*) 'YOU USE ONLY THE DENSITY OF UNIT 9'
      endif                    
 38   DO 13 JATOM=1,NAT                                                 
      JRJ=JRI(JATOM)                                                    
      READ(9,118) LL         
      if(ll.gt.ncom) stop 'increase ncom'                                    
      IF(ICLM.EQ.0) READ(11,118) LL                                     
      LMMAX(JATOM)=LL                                                   
      DO 15 L=1,LL                                                      
      READ(9,2010) LM(1,L,JATOM),LM(2,L,JATOM)                          
       IF(ICLM.NE.0) GOTO 37                                            
      READ(11,2010,END=37) LM(1,L,JATOM),LM(2,L,JATOM)                  
      READ(11,iform1,END=37) (RM(I,1),I=1,JRJ)                             
      READ(9,iform1) (CLM(I,L,JATOM),I=1,JRJ)                              
      READ(9,2031)                                                      
      READ(11,2031)                                                     
      DO 36 I=1,JRJ                                                     
  36  CLM(I,L,JATOM)=CLM(I,L,JATOM)+RM(I,1)*fadd 
      if(deb.eq.'DEBU')  &
        write(6,8687) jatom,l, (CLM(I,L,JATOM),I=jrj-5,JRJ)               
8687   format(2i3,8e15.6)                            
      GOTO 15                                                           
  37  READ(9,iform1) (CLM(I,L,JATOM),I=1,JRJ)                              
      READ(9,2031)                                                      
  15  CONTINUE    
       if(cnorm.eq.'TOT ') then
        do 8484 i=1,jrj
 8484   clm(i,1,jatom)=clm(i,1,jatom)/sqfp 
        end if                                                     
      READ(9,2033)                                                      
      IF(ICLM.EQ.0) READ(11,2033)                                       
      WRITE(6,117) JATOM,LL,(LM(1,L,JATOM),LM(2,L,JATOM),L=1,LL)        
   13 CONTINUE                                                          
!
      do 8686 jatom=1,nat
      jrj=jri(jatom)
      DO 14 I=1,JRJ
      RM(I,JATOM)=RNOT(JATOM)*EXP((I-1)*DX(JATOM))                      
  14  CONTINUE
 8686 RMT(JATOM)=RM(JRJ,JATOM)                                          
!                                                                       
!.... READ FOURIER COEFFICIENTS, REC.LATTIC VECTORS(IN 2*PI/A)          
!.... GENERATE STAR                                                     
!                                                                       
      CALL OUTIN(ICLM,LATTIC,cform,fadd,ortho,br3)                           
  16  CONTINUE                                                          
!                                                                       
!.... GENERATE RECTANGULAR MESH  X = X1 + LAMBDA*A                      
!                                                                       
      NPPY=NPY                                                          
      IF(NPY.EQ.1) NPPY=2                                               
      DO 20 I=1,NPX                                                     
      DO 21 II=1,3                                                      
   21 XI(II)=VEC(II,1)+VX(II)*(I-1)/(NPX-1)                             
      DO 20 J=1,NPY                                                     
      DO 22 II=1,3                                                      
   22 V(II)=VY(II)*(J-1)/(NPPY-1)+XI(II)                                
      IF(SWITCH.EQ.'OVER') GOTO 29                                      
!                                                                       
!.... TEST, IN WHICH SPHERE V FALLS                                     
!                                                                       
!.... LOOP OVER CELL ORIGINS                                            
      DO 25 IPOS=1,NPOS                                                 
!.... LOOP OVER DIFFERENT ATOMS                                         
      DO 25 IAT=1,NDAT                                                  
      DO 26 II=1,3                                                      
   26 VT(II)=V(II)-(ATP(II,IPOS)+POS(II,IAT))                           
      R=VNORM(VT)                                                       
      JATOM=IABS(IATNR(IAT))                                            
   25 IF(R.LT.RMT(JATOM)) GOTO 27                                       
!                                                                       
!.... POINT IN INTERSTITIAL                                             
!                                                                       
!ccc
        DO 28 II=1,3                                                      
      if(ortho) then
        VT(II)=V(II)/A(II)
      else
        VT(II)=V(II)
      endif
 28   continue
      CALL RHOOUT(VT,CHG)  
      if(deb.eq.'DEBU') write(6,*) 'interstitial:',vt,chg                                             
      GOTO 30                                                           
!                                                                       
!.... POINT IN SPHERE IATNR(IAT)                                        
!                                                                       
   27 JATOM=IABS(IATNR(IAT))                                            
!                                                                       
      if(deb.eq.'DEBU') WRITE(6,*) 'JATOM,IAT,POINT, DIFFVECTOR', &
        JATOM,IAT,V,VT           
      VT(1)=V(1)                                                        
      VT(2)=V(2)                                                        
      VT(3)=V(3)                                                        
      if(ortho) then
!     TRANSFER COORDINATES FROM SPHERE IAT TO JATOM                     
         CALL ROTATE(VT,IZ(1,1,IOP(IAT)),TAU(1,IOP(IAT)),A)                
      else
      DO 945 J1=1,3                                                        
 945  Vt1(J1)=BR3(J1,1)*vt(1)+BR3(J1,2)*vt(2)+BR3(J1,3)*vt(3)         
      if(deb.eq.'DEBU') WRITE(6,*) ' VECTOR in hex coord',VT1                 
         CALL ROTATO(VT1,IZ(1,1,IOP(IAT)),TAU(1,IOP(IAT)))                
      if(deb.eq.'DEBU') WRITE(6,*) ' rotated VECTOR in hex coord',VT1
      DO 946 J1=1,3                                                        
 946  Vt(J1)=BR1(1,J1)*vt1(1)+BR1(2,J1)*vt1(2)+BR1(3,J1)*vt1(3)         
      end if
      if(deb.eq.'DEBU') WRITE(6,*) 'ROTATED VECTOR',VT                        
!     REDUCE VECTOR TO SMALLEST POSSIBLE ONE                            
      INDEX=1                                                           
      JATOM1=1                                                          
 31   IF(JATOM1.NE.JATOM) THEN                                          
      INDEX=INDEX+MULT(JATOM1)                                          
      JATOM1=JATOM1+1                                                   
      GOTO 31                                                           
      END IF                                                            
      CALL REDUC(VT,ATP,NPOS,POS,INDEX,rmt(jatom))                                 
      if(deb.eq.'DEBU') WRITE(6,*) ' REDUCED VECTOR',VT                                   
!     TRANSFER COORDINATES INTO LOCAL COORDINATE SYSTEM                 
      CALL ROTAT(VT,ROTLOC(1,1,JATOM))                                  
      if(deb.eq.'DEBU') WRITE(6,*) 'LOCAL VECTOR', VT       
!                                                                       
      IF(R.LT.RNOT(JATOM)) R=RNOT(JATOM)                                      
      IR=1+LOG(R/RNOT(JATOM))/DX(JATOM)                                
      IF(IR.LT.1) IR=1                                                  
      CALL CHARGE(CHG,IR,R,JATOM,VT,IATNR(IAT),clm,lm,lmmax,nat,rm,jri)
      if(deb.eq.'DEBU') write(6,*) 'sphere',chg
!                                                                       
   30 IF(SWITCH.EQ.'RHO ') GOTO 35                                      
   29  CONTINUE                                                         
      CALL CHARG(CHG1,V,NDAT,POS,IATNR,Z,nat)                               
      IF(SWITCH.EQ.'OVER') CHG=CHG1                                     
      IF(SWITCH.EQ.'DIFF') CHG=CHG-CHG1                                 
!                                                                       
   35 CHG=CHG*FAK                                                       
      WRITE(10) CHG                                                     
      IF(NPY.EQ.1) THEN                                                 
      DO 6432 II=1,3                                                    
 6432 V(II)=V(II)-VEC(II,1)                                             
      R=VNORM(V)                                                        
      WRITE(20,'(2f15.7)') R,CHG                                               
      END IF                                                            
      IF(CHG.LT.CMIN) CMIN=CHG                                          
      IF(CHG.GT.CMAX) CMAX=CHG                                          
  20  CONTINUE                                                          
      WRITE(6,121) CMIN,CMAX                                            
  121 FORMAT('0 CMIN=',F10.5,'    CMAX=',F10.5)                         
      RETURN                                                            
  100 FORMAT(6F10.5)                                                    
  101 FORMAT(6I5)                                                       
  102 FORMAT(8A10)                                                      
  103 FORMAT(A4,23X,I3,1x,a4,/,4X,A4)                                         
  104 FORMAT(1X,A10,' LATTIC',10X,I5,' - ATOMS')                        
  105 FORMAT(' LATTIC CONSTANTS:',F10.7,2x,F10.7,2x,F10.7)                               
  106 FORMAT(' BRAVAIS MATRIX:',/,3(10X,3F10.5,/))                      
  107 FORMAT(2I5,3F10.3)                                                
  108 FORMAT('0ATOMNR,OPERATION  POS X    POS Y     POS Z')             
  109 FORMAT('0ENDPOINTS OF PLOT (X,Y,Z):'/)                            
  110 FORMAT('0NUMBER OF PLOTTING POINTS IN X AND Y - DIRECTION:',2I5    &
      ,/' NUMBER OF POINTS IN PRINTOUT:',2I5)                           
  111 FORMAT(1H0,A4,10X,'CHARGE IN:',A4,/)                              
  112 FORMAT(A4,a4,/,A4,a4,a4)                                                   
  113 FORMAT(A10,5X,I5,5X,F10.5,5X,F10.5,5X,F5.2)                       
  114 FORMAT(I4)                                                        
  115 FORMAT(3(3I2,F10.5,/))                                             
  116 FORMAT(3X,5E13.7)                                                 
  117 FORMAT('0ATOM ',I4,'  LMMAX:',I4,'  LM:',11(I3,I2,1X))               
  118 FORMAT(/,15X,I3,//)                                               
  119 FORMAT(1H1,8A10)                                                  
 1012 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/15X,I2)                  
 1013 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
 1051 FORMAT(20X,3F10.8)                                                
 2010 FORMAT(15X,I3,5X,I2,/)                                            
 2031 FORMAT(/)                                                         
 2032 FORMAT(//)                                                        
 2033 FORMAT(///)                                                       
      END                                                               
