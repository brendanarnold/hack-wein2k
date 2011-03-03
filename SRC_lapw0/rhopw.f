      subroutine rhopw(rhok,yka,ABSK,Rhopw0r,jrj,rhopw0)
      IMPLICIT REAL*8     (A-H,O-Z)
      COMPLEX*16 rhok,yka
      real*8 rhopw0(*), rhopw0r(*),bes(7)
      do i=1,jrj+100
        ARG=Rhopw0r(i)*ABSK
        CALL SPHBES(6,ARG,BES)
        rhopw0(i)= rhopw0(i) + rhOK*BES(1)*YKA
      enddo
      return
      end


