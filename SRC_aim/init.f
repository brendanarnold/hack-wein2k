      SUBROUTINE INIT()

!
!     AIM (1999-2000) by Javier D. Fuhr and Jorge O. Sofo
!
!     Instituto Balseiro and Centro Atomico Bariloche
!     S. C. de Bariloche - Rio Negro
!     Argentina
!     e-mail: fuhr@cab.cnea.gov.ar
!             sofo@cab.cnea.gov.ar
!
!     Reads case.struct and case.clmsum files
!
      use rad
      use sphe
      use blank
      use crit
      use zsphe
      use atpos
      implicit none
      INCLUDE 'param.inc'

      real*8 tau,otau,oiz,invoiz,br1,br2,br3,br4
      real*8 a,alpha,beta,gamma
      real*8 gamma1
      real*8 aix,aiy,aiz,rmtt,fadd,sqfp
      real*8 pi,tpi,fpi,sqpi,sqtpi
      real*8 mat,y,yp1,ypn
      real*8 x1,x2,x3,f1,f2,f3,tmp

      integer iz,iord,iszero
      integer i,j,jatom,index,jrj,nat,iform1,mu,i1,j1,ll
      integer l,iclm,i2,k,itmp
      integer ipiv,id,info,nato,ndif

      CHARACTER*4 LATTIC,cform
      CHARACTER*10 TITEL,aname

      logical ortho,deb

      COMMON /DEBUG/ deb
      COMMON /SYM2/ tau(3,NSYM),otau(3,NSYM),oiz(3,3,NSYM), &
         invoiz(3,3,NSYM),iz(3,3,NSYM),iord
!      COMMON /RAD/ RM(NRAD,NATO),JRI(NATO)                       
!      COMMON // CLM(NRAD,NCOM,NATO),CLM2(NRAD,NCOM,NATO), &
!         LM(2,NCOM,NATO),LMMAX(NATO)      
!      COMMON /ATPOS/ ATP(3,nnpos),NPOS,NSA,NSB,NSC
      COMMON /BRAV/ BR1(3,3),BR2(3,3),br3(3,3),br4(3,3)
      COMMON /LATT/ A(3),alpha,beta,gamma,gamma1,ortho
!      COMMON /SPHE/ POS(3,NDIF),RMT(NATO),ROTLOC(3,3,NATO),DX(NATO), &
!         RNOT(NATO),ndat,IATNR(NDIF),IOP(NDIF),MULT(NATO)
      COMMON /CTES/ pi,tpi,fpi,sqpi,sqtpi
!      COMMON /CRIT/ pc(3,NNPOS*NDIF),evpc(3,NNPOS*NDIF), &
!         zpc(3,3,NNPOS*NDIF),pcrb(3,NNPOS*NDIF),npc,icpc(NNPOS*NDIF), &
!         npcm3
!      COMMON /ZSPHE/ ZZ(NDIF)
      DIMENSION TITEL(8)   
      DIMENSION mat(3,3),y(3,3),ipiv(3)
      real*8,allocatable :: z(:)    ! Z(NATO)
      pi=acos(-1.d0)
      tpi=2.d0*pi
      fpi=4.d0*pi
      sqpi=sqrt(pi)
      sqtpi=sqrt(tpi)
      sqfp=2.d0*sqrt(pi)
      npc=0
!     
!     
      READ(8,102) TITEL                                              
      READ(8,103) LATTIC,NAT,cform
      READ(8,100) A,alpha,beta,gamma
      WRITE(6,119) TITEL                                             
      WRITE(6,104) LATTIC,NAT                                        
      WRITE(6,105) A       
      nato=nat
      ndif=nat*48
      allocate( RM(NRAD,NATO),JRI(NATO),z(nato) )                            
      allocate( CLM(NRAD,NCOM,NATO),CLM2(NRAD,NCOM,NATO),LM(2,NCOM,NATO),LMMAX(NATO) ) 
      allocate( POS(3,NDIF),RMT(NATO),ROTLOC(3,3,NATO),DX(NATO),RNOT(NATO) )
      allocate( IATNR(NDIF),IOP(NDIF),MULT(NATO),ZZ(NDIF) )
      assign 2021 to iform1
 2021 FORMAT(3X,4E19.12)

      call gen_brav(LATTIC)
!     
!.... READ NUMBER OF DIFFERENT ATOMS(MAX 17)                         
!.... READ AT.NR (EQ.STRUCT) AND POSITION                            
!     
      WRITE(6,106)((BR1(I,J),I=1,3),J=1,3)                           
      WRITE(6,106)((BR3(I,J),I=1,3),J=1,3)                           
      WRITE(6,106)((BR4(I,J),I=1,3),J=1,3)                           
      WRITE(6,108)                                                   
      INDEX=0                                                        
      DO JATOM=1,NAT                                                  
        INDEX=INDEX+1                                
!          if(index.gt. ndif) stop 'NDIF too small'
        READ(8,1012) IATNR(INDEX),POS(1,INDEX),POS(2,INDEX), &
           POS(3,INDEX),MULT(JATOM)
!       Since later parts of the code use IATNR to indicate which atom is which
!       Force a valid value
        if(iatnr(index).lt.0)then
           iatnr(index)=-jatom
        else
           iatnr(index)=jatom
        endif
!     
!     POS IATNR: MEANS CUB SYM
!     
!     NEG IATNR: MEANS NO CUB SYM
!     
        DO MU=1,MULT(JATOM)-1                                       
          INDEX=INDEX+1
!          if(index.gt. ndif) stop 'NDIF too small'
          READ(8,1013) IATNR(INDEX),POS(1,INDEX), &
             POS(2,INDEX),POS(3,INDEX)  
!       Fix, similar to above
!       Force a valid value
        if(iatnr(index).lt.0)then
           iatnr(index)=-jatom
        else
           iatnr(index)=jatom
        endif
        enddo
!.... READ(RADIAL PARAMETERS, Z, SYMMETRY OPERATIONS                 
!     
        READ(8,113) ANAME,JRI(JATOM),RNOT(JATOM),RMTT,Z(JATOM)
        DX(JATOM)=LOG(RMTT/RNOT(JATOM)) / (JRI(JATOM)-1)
        READ(8,1051) ((ROTLOC(I1,J1,JATOM),I1=1,3),J1=1,3)
           DO MU=0,MULT(JATOM)-1                                      
           zz(index-mu)=z(jatom)
           enddo
      enddo
      NDAT=INDEX 
      READ(8,114) IORD                                               
      DO I=1,IORD                                                    
        READ(8,115) ((IZ(I1,I2,I),I1=1,3),TAU(I2,I),I2=1,3)
      enddo
!     Let's ensure that the first symmetry element is the identity if possible
      iszero=-1
      do I=1,Iord
!       No translations
        if( (abs(TAU(1,I)).gt.1d-12) .or. (abs(tau(2,i)).gt.1d-12) &
                .or. (abs(tau(3,i)).gt.1d-12) ) goto 900
!       Unity on the axes
        if( (IZ(1,1,I) .ne. 1) .or. (IZ(2,2,I) .ne. 1) &
                .or. (IZ(3,3,I) .ne. 1) ) goto 900
!       Angst check
        DO I1=1,3
           DO J1=1,3
                IF( (I1.NE.J1) .and. (IZ(I1,J1,I) .ne. 0) ) goto 900
           enddo
        enddo
!       Jump out
        iszero = I
        goto 910
900     continue
      enddo
!     Not found, abort check
      goto 920
!     Switch so identity is first
910   do I1=1,3
                tmp=tau(I1,1)
                tau(I1,1)       = tau(I1,iszero)
                tau(I1,iszero)  = tmp
                DO J1=1,3
                        itmp            = IZ(I1,J1,1)
                        IZ(I1,J1,1)     = IZ(I1,J1,iszero)
                        IZ(I1,J1,iszero)= itmp
                enddo
      enddo
920   continue
!
!     CREATES MATRICES THAT TRANSFORM "index-ATOM" ONTO "main-ATOM"
!     in orthogonal coordinates
!     OIZ = Tr[BR3] * IZ * BR1
!
!     OTAU
!
      do i1=1,iord
        if (ortho) then
          do i=1,3
            do j=1,3
              oiz(i,j,i1)=1.*iz(i,j,i1)
              invoiz(i,j,i1)=oiz(i,j,i1)
            enddo
          enddo
        else
          do i=1,3
            do j=1,3
              mat(i,j)=0.d0
              do k=1,3
                mat(i,j)=mat(i,j)+iz(i,k,i1)*br1(k,j)
              enddo
            enddo
          enddo
          do i=1,3
            do j=1,3
              oiz(i,j,i1)=0.d0
              do k=1,3
                oiz(i,j,i1)=oiz(i,j,i1)+br4(i,k)*mat(k,j)
                invoiz(i,j,i1)=oiz(i,j,i1)
              enddo
            enddo
          enddo
        endif
        call ludcmp(invoiz(1,1,i1),3,3,ipiv,id,info)
        if(info.ne.0) then
          write(6,*) 'Error inverting oiz:'
          do i=1,3
            write(6,*) (invoiz(i,j,i1),j=1,3)
          enddo
          stop 'ERROR INVERTING OIZ'
        endif
        do  i=1,3
          do j=1,3
            y(i,j)=0.
          enddo 
          y(i,i)=1.
        enddo 
        do j=1,3
          call lubksb(invoiz(1,1,i1),3,3,ipiv,y(1,j))
        enddo         
        do i=1,3
          do j=1,3
            invoiz(i,j,i1)=y(i,j)
          enddo
        enddo
        do i=1,3
          otau(i,i1)=0.d0
          do k=1,3
            otau(i,i1)=otau(i,i1)+tau(k,i1)*br1(k,i)
          enddo
        enddo
      enddo

      CALL ROTDEF(NAT,MULT,IOP,POS,LATTIC)                                  
      DO INDEX=1,NDAT                                                
        AIX=POS(1,INDEX)
        AIY=POS(2,INDEX)
        AIZ=POS(3,INDEX)
        DO J=1,3
          POS(J,INDEX)=BR1(1,J)*AIX+BR1(2,J)*AIY+BR1(3,J)*AIZ
        enddo
        WRITE(6,107) IATNR(INDEX),IOP(INDEX),(POS(J,INDEX),J=1,3)
      enddo
      write (6,*)

!     
!.... READ LMMAX OF ALL ATOMS
!.... READ CLM'S, GENERATE R-MESH
!     
      READ(9,2032)
      DO JATOM=1,NAT
        JRJ=JRI(JATOM)
        do i=1,JRJ
          RM(I,JATOM)=RNOT(JATOM)*EXP((I-1)*DX(JATOM))
        enddo
        RMT(JATOM)=RM(JRJ,JATOM)
        READ(9,118) LL
        LMMAX(JATOM)=LL
        DO L=1,LL
          READ(9,2010) LM(1,L,JATOM),LM(2,L,JATOM)
          if(lm(1,l,jatom).gt. lmax2) stop 'lmax2 too small'
          READ(9,iform1) (CLM(I,L,JATOM),I=1,JRJ)
          READ(9,2031)
        end do    
!        cnorm='TOT'
!        if(cnorm.eq.'TOT ') then
        do i=1,jrj
          clm(i,1,jatom)=clm(i,1,jatom)/sqfp
        enddo
!        end if

        do l=1,ll
           yp1=(clm(2,l,jatom)-clm(1,l,jatom))/(rm(2,jatom)-rm(1,jatom))
           x1=rm(jri(jatom)-2,jatom)
           x2=rm(jri(jatom)-1,jatom)
           x3=rm(jri(jatom),jatom)
           f1=clm(jri(jatom)-2,l,jatom)
           f2=clm(jri(jatom)-1,l,jatom)
           f3=clm(jri(jatom),l,jatom)
           ypn= (-(f3*(x1-x2)*(x1+x2-2*x3))+f2*(x1-x3)*(x1-x3)- &
                f1*(x2-x3)*(x2-x3))/((x1-x2)*(x1-x3)*(x2-x3))
!$$$           ypn=(clm(jri(jatom),l,jatom)-clm(jri(jatom)-1,l,jatom))/
!$$$     $          (rm(jri(jatom),jatom)-rm(jri(jatom)-1,jatom))
!           write(6,*) 'l= ',l,' yp1=',yp1,' ypn=',ypn
          call spline(rm(1,jatom),clm(1,l,jatom),jri(jatom),yp1,ypn, &
             clm2(1,l,jatom))
        enddo
        READ(9,2033)
        WRITE(6,117) JATOM,LL,(LM(1,L,JATOM),LM(2,L ,JATOM),L=1,LL)
      enddo                          
!     
!.... READ FOURIER COEFFICIENTS, REC.LATTIC VECTORS(IN 2*PI/A)
!.... GENERATE STAR
!     
      call OUTIN(ICLM,LATTIC,cform,fadd,ortho,br3,br4)

!.... Endo of initialization
     
      return


 121  FORMAT('0 CMIN=',F12.5,'    CMAX=',F12.5)
 100  FORMAT(6F10.5)
 101  FORMAT(6I5)
 102  FORMAT(8A10)
 103  FORMAT(A4,23X,I3,1x,a4,/,4X,A4)
 104  FORMAT(1X,A10,' LATTIC',10X,I5,' - ATOMS')
 105  FORMAT(' LATTIC CONSTANTS:',F10.7,2x,F10.7,2x,F10.7)
 106  FORMAT(' BRAVAIS MATRIX:',/,3(10X,3F10.5,/))
 107  FORMAT(2I5,3F10.3)
 108  FORMAT('0ATOMNR,OPERATION  POS X    POS Y     POS Z')
 109  FORMAT('0ENDPOINTS OF PLOT (X,Y,Z):'/)
 110  FORMAT &
         ('0NUMBER OF PLOTTING POINTS IN X AND Y - DIRECTION:',2I5, &
          /' NUMBER OF POINTS IN PRINTOUT:' &
         ,2I5)                           
 111  FORMAT(1H0,A4,10X,'CHARGE IN:',A4,/)
 112  FORMAT(A4,a4,/,A4,a4,a4)
 113  FORMAT(A10,5X,I5,5X,F10.5,5X,F10.5,5X,F5.2)
 114  FORMAT(I4)
 115  FORMAT(3(3I2,F10.5,/))
 116  FORMAT(3X,5E13.7)
 117  FORMAT('0ATOM ',I4,'  LMMAX:',I4,'  LM:',11(i3,2X))
 118  FORMAT(/,15X,I3,//)
 119  FORMAT(1H1,8A10)
 1012 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/15X,I2)
 1013 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)
 1051 FORMAT(20X,3F10.8)
 2010 FORMAT(15X,I3,5X,I2,/)
 2031 FORMAT(/)
 2032 FORMAT(//)
 2033 FORMAT(///)
      end

