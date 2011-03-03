SUBROUTINE read_denmat(ifilenum,nat,ndm,jatom1,ll,alx,aly,alz,dmat)
  IMPLICIT NONE
  INTEGER, INTENT(IN)     :: nat,ifilenum
  INTEGER, INTENT(OUT)    :: ndm,jatom1(nat,4,3),ll(nat,4,3)
  REAL*8, INTENT(OUT)     :: alx(nat,4,3),aly(nat,4,3),alz(nat,4,3)
  COMPLEX*16, INTENT(OUT) :: dmat(nat,4,-3:3,-3:3,3)
  INTEGER                 :: inatm,inatmold,isp,ifile,iorb,jatom,m,l
  INTEGER                 :: ios,mp
  ndm=0
  DO isp=1,3
  inatmold=0
  jatom=0
     ifile=ifilenum+isp
     READ(ifile,*,IOSTAT=ios) inatm
     IF(ios/=0.OR.inatm==0) EXIT
     ndm=isp
     REWIND(ifile)
     DO
        READ(ifile,*,IOSTAT=ios) inatm
        IF(ios /= 0.OR.inatm==0) EXIT
        IF (inatm.EQ.inatmold) THEN
           iorb=iorb+1
        ELSE
           iorb=1
           inatmold=inatm
           jatom=jatom+1
        ENDIF
        jatom1(jatom,iorb,isp)=inatm  
        READ(ifile,*) ll(jatom,iorb,isp),alx(jatom,iorb,isp),aly(jatom,iorb,isp),alz(jatom,iorb,isp)   
        l=ll(jatom,iorb,isp)
        DO m=-l,l
           READ(ifile,201)(dmat(jatom,iorb,m,mp,isp),mp=-l,l)
        ENDDO
     ENDDO
  ENDDO
201 format(2(2e16.8,2x))
END SUBROUTINE read_denmat
