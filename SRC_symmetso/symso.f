      subroutine symso(tm,ipot,nA,nB,invers,indso1)
! symso finds out which symmetry operations of the lattice (read in the
! original case.struct file) are preserved in the presence of s-o 
! coupling and  magnetization. M is along (tm1,tm2,tm3) (expressed in
! units of lattice vectors). symso also finds new distribution of
! the atoms into magnetically equivalent groups and then writes new (20)
! structure file.
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
!                                                                       
      CHARACTER*4       LATTIC,IREL,cform,type 
      CHARACTER*10, allocatable::      NAME(:) 
      CHARACTER*40      TITLE 
      character*69      help2
      COMMON /SYM2/ IZ(3,3,noper),TAU(3,noper),IORD                       
      COMMON /SYM2SO/ IZSO(3,3,noper),TAUSO(3,noper),IORDSO             
      COMMON /CHAR   /  LATTIC                           
      INTEGER, allocatable::          IATNR(:), MULT(:)
      real*8   VI
      real*8   A(3), ALPHA(3), PIA(3), br1(3,3)
      real*8, allocatable::   RMT(:), V(:), POS(:,:) 
!      COMMON  /STRUK/    POS, A, ALPHA, RMT, V, PIA, VI, IATNR, MULT
      COMMON  /STRUK/    A, ALPHA, PIA, VI
      SAVE    /STRUK/
!-----------------------------------------------------------------------
      dimension tm(3),rotlso(3,3),rotlo(3,3),roti(3,3)
      INTEGER indso(noper),invmat(3,3),indso1(noper)
      real*8, allocatable:: rotloc(:,:,:),rotij(:,:,:),tauij(:,:)
      INTEGER, allocatable:: multso(:),ninso(:),jrj(:),ind(:)
      INTEGER, allocatable:: ind1(:,:),ind2(:,:),isplit(:),ieq(:)
      real*8, allocatable:: ro(:),zz(:)
! ksym=0 for operation not compatible wit s-o, ksym=1 for compatible op.
! ieq is number of s-o sets originated from non s-o set
      integer ksym(noper)
      common/honza/ SYM(4,3,48),isigk(48)
      dimension sym1(3,3,48)
!-----------------------------------------------------------------------
!     set the inversion matrix
      do i=1,3
      do j=1,3
      if(i.eq.j)then
      invmat(i,j)=-1
      else
      invmat(i,j)=0
      endif
      enddo
      enddo
!     read original structure file
      READ(22,1510) TITLE                                               
      READ(22,3511) LATTIC,NAT,cform,IREL                              
      WRITE(6,1510) TITLE                                               
      write(6,1511) LATTIC,NAT,cform,IREL                              
 17   FORMAT(6F10.6)                                                    
 1011 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
 3011 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
 1012 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)     
 3012 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)     
 1020 FORMAT(///,3X,'ERROR IN EBND  : MULT(JATOM)=0 ...',                &
             /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)           
 1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5)                               
 1510 FORMAT(A40)                                                       
 3511 FORMAT(A4,23X,I3,1x,a4,/,13X,A4,18X,A4)                                 
 1511 FORMAT(A4,23X,I3,1x,a4,/,13X,A4,18X,A4)                                 
!     READ IN LATTICE CONSTANTS                                         
      READ(22,17) A(1),A(2),A(3),(alpha(i),i=1,3)    
      write(6,17) A(1),A(2),A(3),(alpha(i),i=1,3)  

      allocate(NAME(nat),IATNR(NAT),MULT(NAT),RMT(NAT),V(NAT),POS(3,Nat*48))
      allocate(rotloc(3,3,nat),rotij(3,3,Nat*48),tauij(3,Nat*48))
      allocate(multso(Nat*48),ninso(Nat*48),jrj(nat),ind(Nat*48))
      allocate(ind1(nat*48,Nat*48),ind2(nat*48,Nat*48),isplit(nat),ieq(nat))
      allocate(ro(nat),zz(nat))
      ninso=0
      INDEX=0                                                           
      DO 50 JATOM = 1,NAT                                               
         INDEX=INDEX+1
       ind(index)=jatom
         READ(22,3012) IATNR(JATOM),( POS(J,INDEX),J=1,3 ), &
                       MULT(JATOM),ISPLIT(JATOM)
         write(6,1012) IATNR(JATOM),( POS(J,INDEX),J=1,3 ), &
                       MULT(JATOM),ISPLIT(JATOM)
            DO 55 M=1,MULT(JATOM)-1                                     
               INDEX=INDEX+1                                            
       ind(index)=jatom
               READ(22,3011) IATNR(JATOM),( POS(J,INDEX),J=1,3)
               write(6,1011) IATNR(JATOM),( POS(J,INDEX),J=1,3)
  55        CONTINUE
         READ (22,5050) NAME(jatom),JRJ(jatom),RO(jatom),RMT(jatom) &
                       ,ZZ(jatom)
         write(6 ,5051) NAME(jatom),JRJ(jatom),RO(jatom),RMT(jatom) &
                       ,ZZ(jatom)
 5050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)
 5051 FORMAT(A10,' NPT=',I5,'  R0=',F10.9,' RMT=',F10.5,'   Z:',F10.5)
         READ (22,5060) ((ROTLOC(I1,I2,JATOM),I1=1,3),I2=1,3)
         write(6 ,5060) ((ROTLOC(I1,I2,JATOM),I1=1,3),I2=1,3)
 5060 FORMAT(20X,3F10.6)
 50   CONTINUE
      ndif=index 
      READ(22,*) IORD                                                
      DO 29 I=1,IORD                                                    
      do i2=1,3
      READ(22,115) (IZ(I1,I2,I),I1=1,3),TAU(I2,I)               
      enddo
      read(22,*)iii
29    continue
  115 FORMAT(3(3I2,F10.5))                                             
! ***** STRUCT file read, find which symm. oper. invert spin ****
      call angle(tm,a,alpha,lattic,theta,phi)
      call br1dm(br1,lattic,a,alpha)
      call symho(lattic,a,alpha,theta,phi,iord,iz,br1,isigk)
      call latgen(nat,ndif,rotij,tauij,rotloc,pos,mult)
! read in2 file (34), begin to write modified in2 file (54)
      read(34,778,err=779)help2
      write(54,778)help2
      read(34,778,err=779)help2
      write(54,778)help2
      read(34,778,err=779)help2
      write(54,778)help2
      do i=1,nat
      read(34,778,err=779)help2
      enddo
 778  format(a69)
 779  continue
      ncount=index
      iordso=0
      DO 11 I=1,IORD                                                    
      write(6,*)' symmetry operation',i
      do 39 i1=1,3
39    write(6,107) (IZ(I1,I2,I),I2=1,3)               
107   format(3i4)
      if(isigk(i).eq.0)then
      write(6,101)
      ksym(i)=0
      else
      if(isigk(i).eq.1)then
      write(6,102)
      else
      write(6,103)
      endif
101   format('  operation type C does not preserve M direction'/)
102   format('  operation type A preserve M direction'/)
103   format('  operation type B *(time inv.) preserve M direction'/)
      ksym(i)=1
      iordso=iordso+1
      indso(iordso)=i
      do i1=1,3
       do i2=1,3
       sym1(i1,i2,iordso)=sym(i1,i2,i)
       enddo
      enddo
      endif
11    continue
! **************************************************************
! ******** find groups of magnetically equivalent atoms ********
      one=1.d0+0
      toler=0.00001
      natso=0  
      do 47 i=1,ncount
47    multso(i)=0
      DO 30 index=1,ncount
            ni=0
         DO 31 ib=1,natso
         index1=ninso(ib)
            DO 32 K = 1, IORDSO
            i=indso(k)
               X = 0.0D+0
               DO 35 J = 1, 3
                  X = X + iz(J,1,I)*POS(J,INDEX1)
   35          CONTINUE
               X = X + TAU(1,I) + 1.0D+0
               X = MOD(X,ONE)
               Y = 0.0D+0
               DO 34 J = 1, 3
                  Y = Y + iz(J,2,I)*POS(J,INDEX1)
   34          CONTINUE
               Y = Y + TAU(2,I) + 1.0D+0
               Y = MOD(Y,ONE)
               Z = 0.0D+0
               DO 33 J = 1, 3
                  Z = Z + iz(J,3,I)*POS(J,INDEX1)
   33          CONTINUE
               Z = Z + TAU(3,I) + 1.0D+0
               Z = MOD(Z,ONE)
               X1 = ABS(X-POS(1,INDEX))
               Y1 = ABS(Y-POS(2,INDEX))
               Z1 = ABS(Z-POS(3,INDEX))
               IF ((X1 .LT. TOLER) .AND. (Y1 .LT. TOLER) .AND. &
                   (Z1 .LT. TOLER))then
               ni=1 
               goto 37
               endif
!....check positions for centered lattices
            if(lattic(1:1).eq.'B') then
              x1=mod(x1+0.500001d0,one)
              y1=mod(y1+0.500001d0,one)
              z1=mod(z1+0.500001d0,one)
              IF(X1.LT.TOLER.AND.Y1.LT.TOLER.AND.Z1.LT.TOLER) THEN        
              ni=1
              GOTO 37                                                     
              END IF 
            endif                                    
            if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CXY') then
              x1=mod(x1+0.500001d0,one)
              y1=mod(y1+0.500001d0,one)
              IF(X1.LT.TOLER.AND.Y1.LT.TOLER.AND.Z1.LT.TOLER) THEN        
              ni=1
              GOTO 37                                                     
              END IF                                     
              x1=mod(x1+0.500001d0,one)
              y1=mod(y1+0.500001d0,one)
            endif                 
            if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CXZ') then
              x1=mod(x1+0.500001d0,one)
              z1=mod(z1+0.500001d0,one)
              IF(X1.LT.TOLER.AND.Y1.LT.TOLER.AND.Z1.LT.TOLER) THEN        
              ni=1
              GOTO 37                                                     
              END IF                                                      
              x1=mod(x1+0.500001d0,one)
              z1=mod(z1+0.500001d0,one)
            endif                 
            if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CYZ') then
              y1=mod(y1+0.500001d0,one)
              z1=mod(z1+0.500001d0,one)
              IF(X1.LT.TOLER.AND.Y1.LT.TOLER.AND.Z1.LT.TOLER) THEN        
              ni=1
              GOTO 37                                                     
              END IF                                                      
            endif                 
   32       CONTINUE
   31    CONTINUE
   37    CONTINUE
       if(ni.eq.1)then
       write(6,200)index,ib,k
200    format(' atom',i4,' equiv. set',i3,' oper.',i3)
       multso(ib)=multso(ib)+1
! ind1 is index of atom, ind2 is index of cryst. equiv. set
! ind1 relates the atoms in new group to old atom
! ind2 relates the atoms in new group to old atom sets
       ind1(ib,multso(ib))=index
       ind2(ib,multso(ib))=ind(index)
       else
! atom nonequivalent to orig. groups encountered
       natso=natso+1
       multso(natso)=1
       ninso(natso)=index
       write(6,201)index,natso
201    format(' atom',i4,' not equivalent, new set ',i4)
       ind1(natso,1)=index
       ind2(natso,1)=ind(index)
       endif
  30  CONTINUE
! find how many new sets originate from old set
      do i=1,nat
      ieq(i)=0
      enddo
      do 38 j=1,natso
      i1=ind2(j,1)
      ieq(i1)=ieq(i1)+1
      write(6,568)j,i1,ieq(i1)
568   format(' j ind2,ieq',3i4)
38    continue
      write(6,555)natso,(multso(ii),ii=1,natso)
      write(6,556)(ieq(ii),ii=1,nat)
556   format(' ieq',12i4)
555   format(' natso=',i4,/,' multso',12i4,/,' ieq',12i4)
!  ************** nonequiv. atoms end ***************************
! *************** find new local rotation matrices ***************
      call latgen(nat,ndif,rotij,tauij,rotloc,pos,mult)
! ************** write new structure file ************************
 2010 FORMAT('ATOM',I4,': X=',F10.7,' Y=',F10.7,' Z=',F10.7)            
 2011 FORMAT(4X,I4,': X=',F10.7,' Y=',F10.7,' Z=',F10.7)                
 2012 FORMAT('          MULT=',i2,'          ISPLIT=',i2)
 2014 format('LOCAL ROT MATRIX:   ',3f10.7,/,20x,3f10.7,/,20x,3f10.7)
      WRITE(20,1509) TITLE
      WRITE(24,2509) TITLE,(tm(i),i=1,3)
2509  format(a40,' s-o calc. M||',3f6.2) 
1509  format(a40,' intermediate structure file for s-o calc.') 
!cccccccccccc
        pi=acos(-1.d0)
        alpha(1)=alpha(1)*180.d0/pi           
        alpha(2)=alpha(2)*180.d0/pi           
        alpha(3)=alpha(3)*180.d0/pi           
      write(20,1511) LATTIC,natso,cform,IREL                          
      write(20,17) A(1),A(2),A(3),(alpha(i),i=1,3)    
      write(24,1511) LATTIC,natso,cform,IREL                          
      write(24,17) A(1),A(2),A(3),(alpha(i),i=1,3)    
       jatom=0
       do 51 jat=1,nat
       do 51 kat=1,ieq(jat)
       jatom=jatom+1
       i1=ind1(jatom,1)
       i2=ind2(jatom,1)
       do 61 i=1,3
       do 61 j=1,3
       rotlo(i,j)=rotloc(i,j,i2)
       roti(i,j)=rotij(i,j,i1)
61     continue
       call matmm(rotlso,rotlo,roti)
       write(6,*)' jatom,ind1,ind2',jatom,i1,i2
       write(6,*)' rotloc             rotij              rotlso'
       write(6,514)((rotlo(j1,j2),j1=1,3),(roti(j1,j2),j1=1,3) &
                   ,(ROTLSO(j1,j2),j1=1,3),j2=1,3)
514    format(3(3f8.4,3x))
       do 52 jdif=1,multso(jatom)
       i1=ind1(jatom,jdif)
       i2=ind2(jatom,jdif)
! with M present, symmetry is always lower than cubic -> IATNR<0
! but for intermediate structure file IATNR must be the same as original
       if(iatnr(i2).gt.0)then
       iat=jatom
       else
       iat=-jatom
       endif
       if(jdif.eq.1)then
       write(24,2010)-JATOM,( POS(J,I1),J=1,3 )
       write(24,2012)multso(jatom),isplit(i2)
       write(20,2010)IAT,( POS(J,I1),J=1,3 )
       write(20,2012)multso(jatom),isplit(i2)
       else
       write(24,2011)-JATOM,( POS(J,I1),J=1,3 )  
       write(20,2011)IAT,( POS(J,I1),J=1,3 )  
       endif
52     continue
        write(24,5051) NAME(i2),JRJ(i2),RO(i2),RMT(i2) &
                       ,ZZ(i2)
        write(24,2014) ((ROTLSO(j1,j2),j1=1,3),j2=1,3)
        write(20,5051) NAME(i2),JRJ(i2),RO(i2),RMT(i2) &
                       ,ZZ(i2)
        write(20,2014) ((ROTLSO(j1,j2),j1=1,3),j2=1,3)
51     continue
! write the symmetry operations, check,whether the group contains
! center of inversion
      invers=0
 1014    format(i4,6x,'NUMBER OF SYMMETRY OPERATIONS')
      write(20,1014) iord
      do i=1,iord
        do i2=1,3
         write(20,115) (IZ(I1,I2,I),I1=1,3),TAU(I2,I)               
        enddo
      write(20,116)i
      enddo
      write(24,1014) iordso
      l=0
      DO k = 1,iordso
      idelta=0
      i=indso(k)
       if(isigk(i).eq.1)then
        type='   A'
        l=l+1
        indso1(l)=i
        write(6,118)l,i,type
        do i2=1,3
         write(24,115) (IZ(I1,I2,I),I1=1,3),TAU(I2,I)               
         tauso(i2,l)=tau(i2,i)
         do i1=1,3
          izso(i1,i2,l)=iz(i1,i2,i)
          idelta=idelta+(iz(i1,i2,i)-invmat(i1,i2))**2
         enddo
        enddo
        write(24,116)l
        if(idelta.eq.0)invers=1
       endif
      enddo   
      nA=l
      write(6,*)' number of A-type operations=',nA
      DO k = 1,iordso
      i=indso(k)
       if(isigk(i).eq.-1)then
        l=l+1
        type='   B'
        indso1(l)=i
        write(6,118)l,i,type
118     format(' so. oper.',i4,' orig.',i4,' type',a4)
        do i2=1,3
         write(24,115) (IZ(I1,I2,I),I1=1,3),TAU(I2,I)               
         tauso(i2,l)=tau(i2,i)
         do i1=1,3
          izso(i1,i2,l)=iz(i1,i2,i)
          idelta=idelta+(iz(i1,i2,i)-invmat(i1,i2))**2
         enddo
        enddo
        write(24,116)l
        if(idelta.eq.0)invers=1
       endif
      enddo   
      nB=l-nA
      write(6,*)' number of B-type operations=',nB
116   format(i4)
      if(invers.eq.1)then
      write(6,*)' System has center of inversion'
      else
      write(6,*)' System does not have center of inversion'
      endif
! write new VSP files (for i=1,2)   
! write new VNS files (for i=3,4)   
      if(ipot.gt.0)then
      ifr=25
      ifw=45
!....vspup/dn,vnsup/dn
      do 36 i=1,4
      call clmchange(ifr,ifw,nat,ieq,jrj,iord,iz,tau,ksym)
!     call vchange(ifr,ifw,nat,jrj,ieq,iprv)
      ifr=ifr+1
      ifw=ifw+1
36    continue
      ifr=29
      ifw=49
      call in1ch(ifr,ifw,nat,ieq)
      ifr=30
      ifw=50
      call incch(ifr,ifw,nat,ieq)
!...clmsum/up/dn
      ifr=35
      ifw=55
      call clmchange(ifr,ifw,nat,ieq,jrj,iord,iz,tau,ksym)
      ifr=36
      ifw=56
      call clmchange(ifr,ifw,nat,ieq,jrj,iord,iz,tau,ksym)
      ifr=37
      ifw=57
      call clmchange(ifr,ifw,nat,ieq,jrj,iord,iz,tau,ksym)

      endif
      return
      end
