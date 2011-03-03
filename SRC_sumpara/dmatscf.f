SUBROUTINE dmatSCF(SCFFN,nat)
!
!     read scfdmat file and sums up 
!                                                                       
  IMPLICIT REAL*8 (A-H,O-Z)
  INCLUDE 'param.inc'
  !
  COMPLEX*16     :: xtrace
  CHARACTER*79   MARGN                                             
  !
  real*8, allocatable    :: orb(:,:,:),spi(:,:,:),pom(:,:,:,:)
  COMPLEX*16,ALLOCATABLE :: tra(:,:,:)
  real*8, allocatable    :: xop(:,:,:)
  dimension yy(-3:3)
  real*8,allocatable :: dmatr(:,:,:,:),dmati(:,:,:,:)
  
  CHARACTER*80 FNAME,SCFFN
  COMMON /IPROC/   IPROC
  !---------------------------------------------------------------------  
  allocate ( tra(NAT,3,0:3),orb(nat,4,0:3),spi(nat,4,0:3),pom(nat,3,2,0:3))
  allocate ( xop(NAT,0:3,4) )
  allocate ( dmatr(-3:3,-3:3,3,nat),dmati(-3:3,-3:3,3,nat))
  
  dmatr(-3:3,-3:3,1:3,1:nat)=0.d0
  dmati(-3:3,-3:3,1:3,1:nat)=0.d0
  tra(1:nat,1:3,0:3)=(0.d0,0.0d0)
  pom(1:nat,1:3,1:2,0:3)=0.d0
  spi(1:nat,1:4,0:3)=0.d0
  orb(1:nat,1:4,0:3)=0.d0
  xop=0.d0 
  !                                                                       
  !.....READ TAPE 8 = .scfdm UNTIL BOTTOM                                  
  !
  close(22)
  close(8)
  write(6,*) SCFfn
  OPEN(22,FILE=SCFfn,STATUS='unknown',FORM='formatted',err=999)
  
  iprocloop: DO ILOOP=1,IPROC
     katom=0
     WRITE(6,*)'SCFdmat: DOING PROCESSOR FILE: ',ILOOP
     CALL mknam(FNAME,SCFFN,ILOOP)
     write(6,*)'opening ',FNAME
     OPEN(8,FILE=FNAME,STATUS='old',FORM='formatted',err=999)
10   READ(8,700,END=20) MARGN
     write(6,*)MARGN
     if(margn(39:42).eq.'  L=') read(margn,'(42x,i2)') lorb
!print*,'L=',lorb
     IF(MARGN(2:4).EQ.'TRA') THEN
        if(MARGN(20:23).EQ.'UPUP') ispin=1
        if(MARGN(20:23).EQ.'UPDN') ispin=2
        if(MARGN(20:23).EQ.'DNDN') ispin=3
        write(6,*) ispin
        READ(MARGN,'(4x,i3,24x,2f10.5)')IATOM, xtrace
        tra(IATOM,ispin,lorb)=tra(IATOM,ispin,lorb)+Xtrace
        IF(ILOOP.EQ.IPROC)THEN
           if(ispin.eq.1)WRITE(22,719)IATOM,tra(IATOM,ispin,lorb)
           if(ispin.eq.2)WRITE(22,720)IATOM,tra(IATOM,ispin,lorb)
           if(ispin.eq.3)WRITE(22,721)IATOM,tra(IATOM,ispin,lorb)
        ENDIF
     ELSE IF(MARGN(2:4).EQ.'POM') THEN
        ispin=2
        if(MARGN(8:9)=='UP') ispin=1
        READ(MARGN,'(4x,i3,43x,3f9.5)')IATOM, X,y,z
        pom(IATOM,1,ispin,lorb)=pom(IATOM,1,ispin,lorb)+X
        pom(IATOM,2,ispin,lorb)=pom(IATOM,2,ispin,lorb)+y
        pom(IATOM,3,ispin,lorb)=pom(IATOM,3,ispin,lorb)+z
        IF(ILOOP.EQ.IPROC)THEN
           WRITE(22,723)IATOM,MARGN(8:9),(pom(IATOM,ityp,ispin,lorb),ityp=1,3)
        ENDIF
     ELSE IF(MARGN(2:4).EQ.'XOP') THEN
        READ(MARGN,'(4x,i3,i3,4f12.5)')IATOM,ll, X,y,z,zm
        xop(IATOM,ll,1)=xop(IATOM,ll,1)+X
        xop(IATOM,ll,2)=xop(IATOM,ll,2)+y
        xop(IATOM,ll,3)=xop(IATOM,ll,3)+z
        xop(IATOM,ll,4)=xop(IATOM,ll,4)+zm
        IF(ILOOP.EQ.IPROC)THEN
           WRITE(22,724)IATOM,ll,(xop(IATOM,ll,ityp),ityp=1,4)
        ENDIF
     else IF(MARGN(2:4).EQ.'SPI') THEN
        READ(MARGN,'(4x,i3,15x,3f10.5,16x,f9.5)')IATOM, X,y,z,zm
        spi(IATOM,1,lorb)=spi(IATOM,1,lorb)+X
        spi(IATOM,2,lorb)=spi(IATOM,2,lorb)+y
        spi(IATOM,3,lorb)=spi(IATOM,3,lorb)+z
        spi(IATOM,4,lorb)=spi(IATOM,4,lorb)+zm
        IF(ILOOP.EQ.IPROC)THEN
           WRITE(22,710)IATOM,(spi(IATOM,ityp,lorb),ityp=1,4)
        ENDIF
     ELSE IF(MARGN(2:4).EQ.'ORB') THEN
        READ(MARGN,'(4x,i3,18x,3f9.5,16x,f9.5)')IATOM,x,y,z,zm
        orb(IATOM,1,lorb)=orb(IATOM,1,lorb)+x
        orb(IATOM,2,lorb)=orb(IATOM,2,lorb)+y
        orb(IATOM,3,lorb)=orb(IATOM,3,lorb)+z
        orb(IATOM,4,lorb)=orb(IATOM,4,lorb)+zm
        IF(ILOOP.EQ.IPROC) WRITE(22,253)IATOM,(orb(IATOM,ityp,lorb),ityp=1,4)
     ELSE IF(MARGN(17:20).EQ.'UPUP')THEN
        IF(MARGN(29:32).EQ.'real') READ(MARGN,'(42x,i2)') l
        katom=katom+1
        do i=-l,l
           read(8,*) (yy(j),j=-l,l)
           do j=-l,l
              dmatr(i,j,1,katom)=dmatr(i,j,1,katom)+yy(j)
           enddo
        enddo
        read(8,*)
        do i=-l,l
           read(8,*) (yy(j),j=-l,l)
           do j=-l,l
              dmati(i,j,1,katom)=dmati(i,j,1,katom)+yy(j)
           enddo
        enddo
        IF(ILOOP.EQ.IPROC)THEN
           write(22,700) MARGN
           do i=-l,l
              write(22,726) (dmatr(i,j,1,katom),j=-l,l)
           enddo
           write(22,'(" Density matrix UPUP block, imag part ")')
           do i=-l,l
              write(22,726) (dmati(i,j,1,katom),j=-l,l)
           enddo
        ENDIF
     ELSE IF(MARGN(17:20).EQ.'UPDN')THEN
        IF(MARGN(29:32).EQ.'real') READ(MARGN,'(42x,i2)') l
        do i=-l,l
           read(8,*) (yy(j),j=-l,l)
           do j=-l,l
              dmatr(i,j,2,katom)=dmatr(i,j,2,katom)+yy(j)
           enddo
        enddo
        read(8,*)
        do i=-l,l
           read(8,*) (yy(j),j=-l,l)
           do j=-l,l
              dmati(i,j,2,katom)=dmati(i,j,2,katom)+yy(j)
           enddo
        enddo
        IF(ILOOP.EQ.IPROC)THEN
           write(22,700) MARGN
           do i=-l,l
              write(22,726) (dmatr(i,j,2,katom),j=-l,l)
           enddo
           write(22,'(" Density matrix UPDN block, imag part ")')
           do i=-l,l
              write(22,726) (dmati(i,j,2,katom),j=-l,l)
           enddo
        ENDIF
     ELSE IF(MARGN(17:20).EQ.'DNDN')THEN
        IF(MARGN(29:32).EQ.'real') READ(MARGN,'(42x,i2)') l
        do i=-l,l
           read(8,*) (yy(j),j=-l,l)
           do j=-l,l
              dmatr(i,j,3,katom)=dmatr(i,j,3,katom)+yy(j)
           enddo
        enddo
        read(8,*)
        do i=-l,l
           read(8,*) (yy(j),j=-l,l)
           do j=-l,l
              dmati(i,j,3,katom)=dmati(i,j,3,katom)+yy(j)
           enddo
        enddo
        IF(ILOOP.EQ.IPROC)THEN
           write(22,700) MARGN
           do i=-l,l
              write(22,726) (dmatr(i,j,3,katom),j=-l,l)
           enddo
           write(22,'(" Density matrix DNDN block, imag part ")')
           do i=-l,l
              write(22,726) (dmati(i,j,3,katom),j=-l,l)
           enddo
        ENDIF
     ELSE
        IF(ILOOP.EQ.IPROC)THEN
!           call strln(MARGN,LEN)
           len=len_trim(margn)
           IF(LEN.LT.1) LEN=1 
           WRITE(22,700)MARGN(1:LEN)
        ENDIF
     ENDIF
     GOTO 10
20     CLOSE(8)

  ENDDO iprocloop
  
999 continue
  !                                                                       
  RETURN
  
251 FORMAT(':ORB',i3.3,':  UP-ORBITAL MOMENT in global orthog. ', &
         'system=',3f9.5)               
252 FORMAT(':ORB',i3.3,':  DN-ORBITAL MOMENT in global orthog. ', &
         'system=',3f9.5)               
253 FORMAT(':ORB',i3.3,':  ORBITAL MOMENT:',3f9.5,' PROJECTION ON M', &
         F9.5)               
700 FORMAT(A79)        
710 FORMAT(':SPI',i3.3,':','  SPIN MOMENT:', 3F10.5,' PROJECTION ON M',f9.5)         
719 FORMAT(':TRA',i3.3,':','  TRACE of UPUP MATRIX=', &
         2F10.5)         
720 FORMAT(':TRA',i3.3,':','  TRACE of UPDN MATRIX=', &
         2F10.5)         
721 FORMAT(':TRA',i3.3,':','  TRACE of DNDN MATRIX=', &
         2F10.5)         
723 FORMAT(':POM',i3.3,a2,':Partial ORBITAL MOMENT in global orthog. system=',3f9.5)
724 FORMAT(':XOP',i3.3,i3,4f12.5)
725 FORMAT('       DENSITY MATRIX FOR ATOM ',i3,'    L=',i3)
726 FORMAT(7x,7f9.5)
END SUBROUTINE dmatSCF
