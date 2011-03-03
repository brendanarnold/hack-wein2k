SUBROUTINE write_denmat(ifilenum,nat,ndm,jatom1,ll,alx,aly,alz,dmat)
  IMPLICIT NONE
  INTEGER, INTENT(IN)     :: nat,ifilenum
  INTEGER, INTENT(IN)     :: ndm,jatom1(nat,4,3),ll(nat,4,3)
  REAL*8, INTENT(IN)      :: alx(nat,4,3),aly(nat,4,3),alz(nat,4,3)
  COMPLEX*16, INTENT(IN)  :: dmat(nat,4,-3:3,-3:3,3) 
  INTEGER                 :: idm,ifile,jatom,iorb,l,m,mp
  DO idm=1,ndm
     ifile=ifilenum+idm
     DO jatom=1,nat
        IF(jatom1(jatom,1,idm).EQ.0) EXIT
        DO iorb=1,4
           IF(jatom1(jatom,iorb,idm).EQ.0) EXIT      
           WRITE(ifile,100) jatom1(jatom,iorb,idm)
           WRITE(ifile,102) ll(jatom,iorb,idm),alx(jatom,iorb,idm),aly(jatom,iorb,idm),alz(jatom,iorb,idm)   
           l=ll(jatom,iorb,idm)
           DO m=-l,l
              WRITE(ifile,201)(dmat(jatom,iorb,m,mp,idm),mp=-l,l)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
100   format(i5,' atom density matrix')
102   format(i5,3f10.6, ' L, Lx,Ly,Lz in global orthogonal system')
201 format(2(2e16.8,2x))
END SUBROUTINE write_denmat
