      SUBROUTINE RECPR (NWAVE,INDMAX,KXMAX,KYMAX,KZMAX,GMAX)            
!                                                                       
      USE param
      use defs
      use struk
      use reclat
      use reallocate
      USE parallel
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER    IM(3),NST,ISTM(3,NSYM)
      COMPLEX*16       TAUP(NSYM)                                     
      CHARACTER*67     ERRMSG
      LOGICAL          KDELTA
!                                                                       
      REAL*8, pointer             :: ABSK(:)
      DIMENSION        AM(3),M(3),KMAX(3)                    
      DATA DELTA/0.01D0/
      real*4,  allocatable :: radii(:),iradii(:)
      integer, allocatable :: ihkl1(:,:),ihkl2(:,:),icnt(:)

!---------------------------------------------------------------------  
!                                       
      ABMAX=-999
      if(myid.eq.0) WRITE(6,*) ' KXMAX,KYMAX,KZMAX',KXMAX,KYMAX,KZMAX
       KMAX(1)=1                                                     
       KMAX(2)=1                                                     
       KMAX(3)=1                                                     
      LL1=KXMAX+1                                                       
      LL2=KYMAX+1                                                       
      LL3=KZMAX+1                                                       
      TEST1=1.D-05    
       NWAV=11111                                                   
       NWAV1=NWAV                                                        
      allocate ( kzz(3,nwav), absk(nwav) )
      absk=0.0d0
!                                                                       
      ILL1=-LL1+2                                                           
      ILL2=-LL2+2                                                           
      ILL3=-LL3+2
!     Presort values
      irad=0
      NRADM=(LL1-ILL1+1)*(LL2-ILL2+1)*(LL3-ILL3+1)+1
      allocate (icnt(NRADM),radii(NRADM),ihkl1(3,NRADM),iradii(nradm))
!     Find all possible radii
      DO 199 I1=ILL1,LL1  
      M(1)=I1-1                                                         
      DO 199 I2=ILL2,LL2                                                     
      M(2)=I2-1                                                         
      DO 199 I3=ILL3,LL3                                                
      M(3)=I3-1                                                         
!                                                                       
      ABSM=0.0D0                                                        
      DO IND1=1,3                                                    
         AMIND1=0.0                                                   
         DO  IND2=1,3                                                 
             AMIND1=AMIND1+M(IND2)*BR2(IND1,IND2)
         ENDDO                       
         ABSM=ABSM+AMIND1*AMIND1                                                        
      ENDDO                                                          
                                                                       
      if(absm.gt.1d-8)ABSM=SQRT(ABSM)                                                   
      if(absm.le.gmax) then
                irad=irad+1
                ICNT(irad)=IRAD
                RADII(IRAD)=ABSM
                ihkl1(1,irad)=M(1)
                ihkl1(2,irad)=M(2)
                ihkl1(3,irad)=M(3)
      endif 
  199 CONTINUE                         
!     Sort the radii, preserving irad order
!     Mimic this
!     Quicksort
      call sortag(radii,irad,icnt)
!     Now re-order equal radii
      do j=1,irad
        iradii(j)=icnt(j)
      enddo
      ilow=2
      iup=ilow
      imcheck=-999
1000  rtest=radii(ilow)
      iup=iup+1
      if(iup.eq.irad)goto 1001
      if(radii(iup)-rtest.gt.test)then
!       radii are different, sort
        iuse=iup-ilow
        if(iuse.gt.1)then
                call sortag(iradii(ilow),iuse,icnt(ilow))
!               Keep track of how many we had to sort over for later
                imcheck=max(imcheck,iuse)
        endif
        ilow=iup
      endif
      goto 1000
1001  continue
!     Anything left over ?
      iup=irad
      iuse=iup-ilow
      if(iuse.gt.1)call sortag(iradii(ilow),iuse,icnt(ilow))
!     Keep track of how many we had to sort over for later
      imcheck=max(imcheck,iuse)

!     Copy a sorted hkl list from ihkl1 to ihkl2
      allocate (ihkl2(3,irad))
      do j=1,irad
        i=icnt(j)
        ihkl2(1,j)=ihkl1(1,i)
        ihkl2(2,j)=ihkl1(2,i)
        ihkl2(3,j)=ihkl1(3,i)
      enddo
      deallocate (ihkl1)
      deallocate (radii)
      deallocate (iradii)    
!
      J=1                                                               
      NK=0 
!      DO 1 I1=ILL1,LL1                                                 
!      M(1)=I1-1                                                         
!      DO 1 I2=ILL2,LL2                                                     
!      M(2)=I2-1                                                         
!      DO 1 I3=ILL3,LL3                                                
!      M(3)=I3-1                                                         
       DO 1 INEW=1,IRAD
        M(1)=ihkl2(1,INEW)
        M(2)=ihkl2(2,INEW)
        M(3)=ihkl2(3,INEW)
!                                                                       
!   THE NEGATIVE INDICES ARE NECESSARY IN ORDER TO OBTAIN K-VECTORS     
!   WITH 0 COMPONENTS AFTER MULTIPLICATION WITH THE BRAVAIS MATRIX.     
!
      ABSM=0.0D0                                                        
      DO 22 IND1=1,3                                                    
         AM(IND1)=0.0                                                   
         DO 23 IND2=1,3                                                 
  23     AM(IND1)=AM(IND1)+M(IND2)*BR2(IND1,IND2)                       
         ABSM=ABSM+AM(IND1)**2                                          
         AHELP=AM(IND1)/PIA(IND1)                                       
         IM(IND1)=AHELP+SIGN(DELTA,AHELP)                               
         IF(ORTHO.or.lattic(1:3).eq.'CXZ') GOTO 22                            
         IM(1)=M(1)                                                     
         IM(2)=M(2)                                                     
         IM(3)=M(3)                                                     
   22 CONTINUE                                                          
!  transformation into primitiv monoclinic basis
      IF(.not.ORTHO.and.lattic(1:3).eq.'CXZ') then
        im(1)=m(1)+m(3)
        im(2)=m(2)
        im(3)=-m(1)+m(3)
      endif 
!                                                                       
      ABSM=SQRT(ABSM)                                                   
      IF(NK.EQ.0) GOTO 44
!     Skip search if larger than any others
      if(ABSM.GT.ABMAX)goto 44                                                

      DO 2 J=max(NK-imcheck-1,1),NK                                                       
      ABST=ABSM-ABSK(J)                                                 
      IF(ABS(ABST).LT.TEST)GOTO 5                                       
   2  CONTINUE                                                          
! .... NEW VECTOR Different                                           
201   goto 44

! .... NEW VECTOR EQU. OLD ONE ?
  5   CALL STERN(IM,NST,ISTM,TAUP)
      DO 3 I=J,NK                                                       
       IF(KDELTA(KZZ(1,I),IM,NST,ISTM)) goto 1
  3   CONTINUE                                                          

!     New vector
 44   nk=nk+1
      j=nK
      IF(J.GE.NWAV1) then
        nwav1=nwav1*3
        call doreallocate(kzz, 3, nwav1)
        call doreallocate(absk, nwav1)
      endif                                      
      ABSK(J)=ABSM                                                      
      ABMAX=MAX(ABSM,ABMAX)+2.D0*TEST
! .... PUT NEW VECTOR IN LIST                                           
      DO 8 N=1,3                                                        
   8  KZZ(N,J)=IM(N)                         

   1  CONTINUE                                                          
      IND=0                                                             
      deallocate (ihkl2)
      deallocate (icnt)
!                                                                       
      if(myid.eq.0) WRITE(6,1010) Nk                                      
!     reduce dimensions to actual value
         if(myid.eq.0) write(6,*) 'nwav1,kn',nwav1,nk
	call doreallocate(kzz, 3, nk)
	call doreallocate(absk, nk)
      allocate ( inst(nk) )
      allocate ( tauk(nk*nsym) )
      kxmax=0
      kymax=0
      kzmax=0
      NWAVE=Nk                                                      
!      NWAVE0=Nk                                                      
      do i=1,nk
       do j=1,3
            IF (IABS(KZZ(J,I)).GT.KMAX(J)) kmax(j)=IABS(KZZ(J,I))
       enddo
      enddo                    
!      allocate ( krcval(-kmax(1):kmax(1),-kmax(2):kmax(2),-kmax(3):kmax(3)) )
!
      DO 30 I=1,Nk                                                  
         DO 31 J=1,3                                                    
!            IF (IABS(KZZ(J,I)).GT.KMAX(J)) THEN                         
!                WRITE(6,1000) I,KZZ(1,I),KZZ(2,I),KZZ(3,I),KMAX(J)      
!                NWAVE=I-1                                               
!                GO TO 900                                                
!            END IF                                                      
            IM(J)=KZZ(J,I)                                              
 31      CONTINUE                                                       
         CALL STERN(IM,NST,ISTM,TAUP)
         INST(I)=NST                                                    
         if(myid.eq.0) WRITE(6,1020) I,KZZ(1,I),KZZ(2,I),KZZ(3,I),absk(i),nst             
!         if(absk(i).gt.gmin) then
!           nwave0=i-1
!           gmin=1d10
!         endif
         DO 32 JJ=1,NST   
!            write(6,1021)istm(1,jj),istm(2,jj),istm(3,jj),taup(jj)          
         if(iabs(istm(1,jj)).gt.kxmax) kxmax=iabs(istm(1,jj))
         if(iabs(istm(2,jj)).gt.kymax) kymax=iabs(istm(2,jj))
         if(iabs(istm(3,jj)).gt.kzmax) kzmax=iabs(istm(3,jj))
            IND=IND+1                                                   
            TAUK(IND)=TAUP(JJ)                                          
!c            KRCVAL( ISTM(1,JJ),ISTM(2,JJ),ISTM(3,JJ) )=IND              
 32      CONTINUE                                                       
 30   CONTINUE                                                          
 33   INDMAX=IND                                                        
      if(myid.eq.0) WRITE(6,1030) INDMAX                                  
      nwave1=nwave
      deallocate (absk)
      RETURN                                                            
!  900 WRITE (6,*) 'KMAX1, KMAX2, KMAX3 too small'
!      WRITE(6,9000) I,KZZ(1,I),KZZ(2,I),KZZ(3,I)
!      STOP 'RECPR - Error'
!                                                                       
 100  FORMAT(I4)                                                        
! 1000 FORMAT(///,3X,'WARNING IN RECPR : ',I4,'TH VECTOR:',3I4,          
!     *  ' GT KMAX',                                                     
!     *       /,20X, 'KMAX=',I3,10X,' VECTOR LIST TRUNCATED')            
 1010 FORMAT(3X,I10,' PLANE WAVES GENERATED (INCLUDING FORBIDDEN ',       &
             'H,K,L)',/)                                                
 1020 FORMAT(7X,'KVEC(',I10,') = ',3I5,f10.4,i5)   
 1021 FORMAT(8X,'             ',3i5,2f10.5)
 1030 FORMAT(8X,'SIZE INCLUDING STAR MEMBERS = ',I10)                    
! 9000 FORMAT(I4,'TH VECTOR:',3I4)          
      END                                                               


