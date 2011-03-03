      subroutine latsym(nsym,amat,atrans,nat)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
                                                                        
      LOGICAL           ORTHO
!                                                                       
      CHARACTER*4       LATTIC                                          
      CHARACTER*10      KNAME                                      
      CHARACTER*80      TITLE                                           
!
      COMMON /SYM2/ IZ(3,3,48),TAU(3,48),IORD                       
      COMMON /CHAR   /  LATTIC                                          
      COMMON /GENER  /  BR2(3,3)
                                     
!-----------------------------------------------------------------------
      Dimension SYM(4,3,48),a1(3,3),a2(3,3),b1(3,3),trans(3)
      INTEGER,allocatable :: Iatnr(:),INDXa(:),IAS(:,:),iatom(:)
      dimension amat(3,3,NOPER),a(3),atrans(3,NOPER)
      real*8,allocatable :: pos(:,:)
      CHARACTER*10,allocatable :: name(:)
      INTEGER,allocatable :: index1(:),index2(:),mult(:)
!-----------------------------------------------------------------------
      READ(20,1510) TITLE                                               
      WRITE(6,1510) TITLE                                               
      READ(20,1511) LATTIC,NAT                                          
 17   FORMAT(6F10.7)                                                    
 1011 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
 1012 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,7X,I2,8X,I2)     
 1020 FORMAT(///,3X,'ERROR IN EBND  : MULT(JATOM)=0 ...',                &
             /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)           
 1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5)                               
 1510 FORMAT(A80)                                                       
 1511 FORMAT(A4,23X,I3,/,13X,A4)         
      allocate (iatom(nat*48*16),Iatnr(nat*48*16),INDXa(Nat*48*16),IAS(48,Nat*48*16))
      allocate (pos(3,nat*48*16))
      allocate ( name(nat),index1(nat),index2(nat),mult(nat))
!     READ IN LATTICE CONSTANTS                                         
      READ(20,17) A(1),A(2),A(3),alpha,beta,gamma                          
      if(alpha.eq.0.d0) alpha=90.d0                                       
      if(beta .eq.0.d0) beta =90.d0                                       
      if(gamma.eq.0.d0) gamma=90.d0     
      if(lattic(1:1).eq.'H'.or.lattic(1:1).eq.'R') gamma=90.d0
      pi=acos(-1.d0)                                  
!.....SET UP LATTICE, AND IF REQUIRED ROTATION MATRICES                 
      CALL DIRLAT (NAT,ortho,alpha,beta,gamma)                                
!     INDEX COUNTS ALL ATOMS IN THE UNIT CELL, JATOM COUNTS ONLY THE    
!     NONEQUIVALENT ATOMS                                               
      INDEX=0                                                           
      DO 50 JATOM = 1,NAT                                               
         INDEX=INDEX+1
         READ(20,1012) IATNR(index),( POS(J,INDEX),J=1,3 ),MULT(JATOM)  
         IF (MULT(JATOM).EQ.0) THEN                                     
            WRITE(6,1020) JATOM,INDEX,MULT(JATOM)                       
            STOP ' NNN: MULT EQ 0'                                      
         ENDIF                     
            DO 55 M=1,MULT(JATOM)-1                                     
               INDEX=INDEX+1                                            
               READ(20,1011) IATNR(index),( POS(J,INDEX),J=1,3)
  55        CONTINUE
         READ(20,1050) NAME(JATOM) 
         read(20,*)           
         read(20,*)           
         read(20,*)           
 50   CONTINUE 
      READ(20,*) IORD                                                
      DO 11 I=1,IORD                                                    
   11 READ(20,115) ((IZ(I1,I2,I),I2=1,3),TAU(I1,I),I1=1,3)               
  115 FORMAT(3(3I2,F10.5,/))                                             
!     ALL OF POTE IS READ                                               
!....scaling and redefinition: R_i= BR2(*,i)
!....scaling special for all monoclinic lattices
      if(gamma.ne.90.d0) then
      br2(1,1)=br2(1,1)*a(1)
      br2(2,2)=br2(2,2)*a(2)
      br2(3,3)=br2(3,3)*a(3)
      br2(1,2)=br2(1,2)*a(1)
      br2(1,3)=br2(1,3)*a(3)
      br2(3,1)=br2(3,1)*a(1)
      br2(3,2)=br2(3,2)*a(1)
      else
      DO 140 I=1,3
      DO 140 J=1,3
         Br2(J,I)=Br2(J,I)*a(i)
  140 CONTINUE
      end if
      DO 141 I=1,3
      DO 141 J=1,i
         ahelp=br2(j,i)
         Br2(J,I)=Br2(i,J)
         br2(i,j)=ahelp
 141  CONTINUE
      write(6,*) ' ATOMIC POSITIONS:'
      IF (ortho) THEN
         DO 160 IA=1,index
         DO 161 J=1,3
 161        pos(J,IA)=pos(J,IA)*a(j)
             write(6,8181) ia,pos(1,ia),pos(2,ia),pos(3,ia)
  160    CONTINUE
      ELSE 
         DO 170 IA=1,index
            TAU1=pos(1,IA)
            TAU2=pos(2,IA)
            TAU3=pos(3,IA)
            if(lattic(1:3).eq.'CXZ') then 
!...........positions are in primitiv monoclinic coord.
              pos(1,IA)=tau1*sin(Gamma/180.d0*pi)*a(1)
              pos(2,IA)=tau1*cos(Gamma/180.d0*pi)*a(1)+tau2*a(2)
              pos(3,IA)=tau3*a(3)
            else
            pos(1,IA)=Br2(1,1)*TAU1+Br2(1,2)*TAU2+Br2(1,3)*TAU3
            pos(2,IA)=Br2(2,1)*TAU1+Br2(2,2)*TAU2+Br2(2,3)*TAU3
            pos(3,IA)=Br2(3,1)*TAU1+Br2(3,2)*TAU2+Br2(3,3)*TAU3
!            TAU(1,IA)=B(1,1)*TAU1+B(2,1)*TAU2+B(1,3)*TAU3
!            TAU(2,IA)=B(1,2)*TAU1+B(2,2)*TAU2+B(2,3)*TAU3
!            TAU(3,IA)=B(3,1)*TAU1+B(3,2)*TAU2+B(3,3)*TAU3
             endif
             write(6,8181) ia,pos(1,ia),pos(2,ia),pos(3,ia)
  170    CONTINUE
      END IF
 8181 format(i5,3f15.8)
! FIND THE SYMMETRY OF THE SYSTEM
!      write(6,*) br2
      CALL PGLSYM(Br2,SYM,NSYM)
      CALL PGBSYM(Br2,pos,Iatnr,index,SYM,NSYM,a)
      if(nsym.ne.iord.and.iord.ne.0) write(6,8182) nsym,iord
 8182 format(//,'        WARNING !!!!!!',/, &
      'nsym found by symmetry differs from iord read in struct',2i4)
!      if(iord.eq.0) then
!         backspace 20
!         write(20,'(i4,"    number of symmetry operations")') nsym
!      endif
!
!
             CALL RECLAT (br2,a2,0)
      DO 143 I=1,3
      DO 143 J=1,i
         ahelp=a2(j,i)
         a2(J,I)=a2(i,J)
         a2(i,j)=ahelp
 143  CONTINUE
             call matmm(a1,br2,a2)
!             write(6,111) (a2(1,j),j=1,3)
!             write(6,111) (a2(2,j),j=1,3)
!             write(6,111) (a2(3,j),j=1,3)
!             write(6,111) (br2(1,j),j=1,3)
!             write(6,111) (br2(2,j),j=1,3)
!             write(6,111) (br2(3,j),j=1,3)
!             write(6,111) (a1(1,j),j=1,3)
!             write(6,111) (a1(2,j),j=1,3)
!             write(6,111) (a1(3,j),j=1,3)
       do 887 i=1,nsym
       write(6,*) 'Symmetry operation',i
       do 888 i1=1,4
       write(6,8183) (sym(i1,i2,i),i2=1,3)
 8183  format(3f15.5)
 888   continue
         do 2 i1=1,3
         trans(i1)=sym(4,i1,i)/a(i1)
         if(trans(i1).lt.-1.d-4) trans(i1)=trans(i1)+1.d0
         do 2 j=1,3
!ccc  2      a1(i1,j)=sym(i1,j,i)
  2      a1(j,i1)=sym(i1,j,i)
!...symmetry op. for mon.CXZ lattice in simple mon.coord. 
             if(.not.ortho.and.lattic(1:3).ne.'CXZ') then
             call matmm(b1,a2,a1)
             call matmm(a1,b1,br2)
!
             do 765 i1=1,3
             trans(i1)=sym(4,1,i)*a2(i1,1)+sym(4,2,i)*a2(i1,2)+ &
                       sym(4,3,i)*a2(i1,3)
         if(trans(i1).lt.-1.d-4) trans(i1)=trans(i1)+1.d0
 765         continue
             end if
             write(6,111) (a1(1,j),j=1,3),trans(1)
             write(6,111) (a1(2,j),j=1,3),trans(2)
             write(6,111) (a1(3,j),j=1,3),trans(3)
!             write(6,111) (a1(j,1),j=1,3),trans(1)
!             write(6,111) (a1(j,2),j=1,3),trans(2)
!             write(6,111) (a1(j,3),j=1,3),trans(3)
111   format(4f10.4)
 112  format(3i2,f10.8)
      if(iord.gt.0) then
             call checks(a1,trans)
!.....write symmetry operations to struct file
!      else
!             write(20,112) (nint(a1(1,j)),j=1,3),trans(1)
!             write(20,112) (nint(a1(2,j)),j=1,3),trans(2)
!             write(20,112) (nint(a1(3,j)),j=1,3),trans(3)
!             write(20,'(i4)') i
      end if
       do j=1,3
       atrans(j,i)=trans(j)
         do j1=1,3
         amat(j1,j,i)=a1(j1,j)
         enddo
       enddo
 887   continue
!
      return
      end
