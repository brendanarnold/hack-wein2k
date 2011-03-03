      SUBROUTINE RECPR (INDMAX,kxmax,kymax,kzmax &
        ,gmin,gmax,nwave0)
!                                                                       
      use struct
      use waves
      use reallocate
      IMPLICIT REAL*8 (A-H,O-Z)
!
!      INCLUDE 'param.inc'
!      COMPLEX*16       TAUP,TAUK                                     
      COMPLEX*16       TAUP                                     
!      CHARACTER*4      LATTIC                                           
!      CHARACTER*5      MODUS                                            
!      CHARACTER*10     ANAME                                            
      CHARACTER*67     ERRMSG
!      CHARACTER*80     TITLE                                            
      LOGICAL          KDELTA,ORTHO
!                                                                       
      COMMON /ORTH/   ORTHO                                            
      COMMON /GENER/  BR1(3,3),BR2(3,3)                                 
!      COMMON /STRUK/  POS(3,NDIF),AA,BB,CC,ALPHA(3),RMT(NATO),V(NATO),
!     *                PIA(3),VOL,IATNR(NATO),MULT(NATO),ISPLIT(NATO)           
!      COMMON /CHAR/   TITLE,LATTIC,MODUS,ANAME(NATO)                    
!      COMMON /GITT/   KZZ(3,NWAV),ABSK(NWAV),NWAVE                            
      COMMON /FOUR/   IM(3),NST,ISTM(3,NSYM),TAUP(NSYM)                 
!      COMMON /SYM2/   TAU(3,NSYM),IORD,IZ(3,3,NSYM)                     
!      COMMON /INVERS/ TAUK(NWAV*NSYM),                                  
!     *                KRCVAL(-KMAX1:KMAX1,-KMAX2:KMAX2,-KMAX3:KMAX3)    
!      DIMENSION        AM(3),M(3),INST(NWAV),KMAX(3)                    
      DIMENSION        AM(3),M(3),KMAX(3)                    
      DATA DELTA/0.01D0/
      real*4,  allocatable :: radii(:),iradii(:)
      integer, allocatable :: ihkl1(:,:),ihkl2(:,:),icnt(:)
!---------------------------------------------------------------------  
!                                       
      WRITE(6,*) ' KXMAX,KYMAX,KZMAX',KXMAX,KYMAX,KZMAX
       KMAX(1)=1                                                     
       KMAX(2)=1                                                     
       KMAX(3)=1                                                     
      LL1=KXMAX+1                                                       
      LL2=KYMAX+1                                                       
      LL3=KZMAX+1                                                       
      TEST=1.D-05    
       NWAV=11111                                                   
       NWAV1=NWAV                                                        
      allocate ( kzz(3,nwav), absk(nwav) )
!                                                                       
      ILL1=-LL1+2                                                           
      ILL2=-LL2+2                                                           
      ILL3=-LL3+2    
      abmax=-999

!      write(6,*)'Debug'
!      do j=1,3
!        write(6,601)(BR2(J,K),K=1,3)
!      enddo
601   format(3f20.15)
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
!     Sort the radii
!     What is simplest is to write then to a file then use sort -n
!      open(unit=98,file='trap1')
!      do j=1,irad
!        write(98,981)radii(j),j
!981     format(F14.5,i10)
!      enddo
!      close(unit=98)
!      call system('sort -n trap1 > trap2 ')
!      open(unit=98,file='trap2')
!      do j=1,irad
!        read(98,*)t,icnt(j)
!      enddo
!     The above sorts both the radii, preserving irad order, but is not portable
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
!      open(unit=98,file='checker')
1000  rtest=radii(ilow)
      iup=iup+1
      if(iup.eq.irad)goto 1001
      if(radii(iup)-rtest.gt.test)then
!       radii are different, sort up to here
        iuse=iup-ilow
        if(iuse.gt.1)then
!                write(98,*)'Sorting from ',ilow,iup-1
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

!        do j=1,irad
!                i=icnt(j)
!98      format(4i5,f10.5)
!                write(98,98)i,ihkl1(1,i),ihkl1(2,i),ihkl1(3,i),radii(j)
!        enddo
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
!         AM(IND1)=0.0                                                   
!         DO 23 IND2=1,3                                                 
!  23     AM(IND1)=AM(IND1)+M(IND2)*BR2(IND1,IND2)                       
!         ABSM=ABSM+AM(IND1)**2                                          
         AMIND1=0.D0
         DO 23 IND2=1,3
   23           AMIND1=AMIND1+DBLE(M(IND2))*BR2(IND1,IND2)
         ABSM=ABSM+AMIND1*AMIND1
!         AHELP=AM(IND1)/PIA(IND1)                                       
         AHELP=AMIND1/PIA(IND1)
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
!
!     Use imcheck+1 as how far back to check at most
      DO 2 J=max(NK-imcheck-1,1),NK                                                       
      ABST=ABSM-ABSK(J)                                                 
      IF(ABS(ABST).LT.TEST)GOTO 5                                       
   2  CONTINUE                                                          
! .... NEW VECTOR Different                                           
201   goto 44

! .... NEW VECTOR EQU. OLD ONE ?
  5   CALL STERN                                                        
      DO 3 I=J,NK                                                       
       IF(KDELTA(KZZ(1,I))) goto 1
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
!                                                                       
      WRITE(6,1010) Nk                                              
!     reduce dimensions to actual value
         write(6,*) 'nwav1,kn',nwav1,nk
	call doreallocate(kzz, 3, nk)
	call doreallocate(absk, nk)
      allocate ( inst(nk) )
!     this is not needed in dstart
!      allocate ( tauk(nk*nsym) )
      kxmax=0
      kymax=0
      kzmax=0
      NWAVE=Nk                                                      
      NWAVE0=Nk                                                      
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
         CALL STERN                                                     
         INST(I)=NST                                                    
         WRITE(6,1020) I,KZZ(1,I),KZZ(2,I),KZZ(3,I),absk(i),nst             
         if(absk(i).gt.gmin) then
           nwave0=i-1
           gmin=1d10
         endif
         DO 32 JJ=1,NST   
!            write(6,1021)istm(1,jj),istm(2,jj),istm(3,jj),taup(jj)          
         if(iabs(istm(1,jj)).gt.kxmax) kxmax=iabs(istm(1,jj))
         if(iabs(istm(2,jj)).gt.kymax) kymax=iabs(istm(2,jj))
         if(iabs(istm(3,jj)).gt.kzmax) kzmax=iabs(istm(3,jj))
            IND=IND+1                                                   
!c            TAUK(IND)=TAUP(JJ)                                          
!c            KRCVAL( ISTM(1,JJ),ISTM(2,JJ),ISTM(3,JJ) )=IND              
 32      CONTINUE                                                       
 30   CONTINUE                                                          
 33   INDMAX=IND                                                        
      WRITE(6,1030) INDMAX                                              
      nwave1=nwave
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
 1020 FORMAT(7X,'KVEC(',I10,') = ',3I5,f12.6,i5)   
 1021 FORMAT(8X,'             ',3i5,2f10.5)
 1030 FORMAT(8X,'SIZE INCLUDING STAR MEMBERS = ',I10)                    
! 9000 FORMAT(I4,'TH VECTOR:',3I4)          
      END                                                               
