      subroutine normal(rhok,nkk,jspin,nat,mult,zz,qpw,qel,qmx,margn)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      CHARACTER*67       ERRMSG
      CHARACTER*79       MARGN
      complex*16 rhok(nwav,2)
      dimension qpw(2),qel(nato,2),zz(nato),mult(nato)
!
!.....renormalize charge density
!
      q=0.d0
      qpw1=0.d0
      do 20 ispin=1,jspin
         z=0.d0
         q=q+qpw(ispin)
         qpw1=qpw1+qpw(ispin)
            do 20 jatom=1,nat
            z=z+zz(jatom)*mult(jatom)
            q=q+qel(jatom,ispin)*mult(jatom)
 20   continue
!      if ((qmx.gt.1.0D-5).and.(ABS(z-q).gt.0.5d0))
!     *   goto 900
      write(6,2042)  margn,z,q,z/q
 2042 format(/,':NEC',a2,':',1x,'NUCLEAR AND ELECTRONIC CHARGE ',3f10.5)
      qfix=(qpw1+z-q)/qpw1
      do 30 ispin=1,jspin
      do 30 i=1,nkk
 30      rhok(i,ispin)=rhok(i,ispin)*qfix
      return
  900 CALL OUTERR('NORMAL','nuclear and electronic charge differs')
      WRITE (ERRMSG,9000) Z
      CALL OUTERR('NORMAL',ERRMSG)
      WRITE (ERRMSG,9010) Q
      CALL OUTERR('NORMAL',ERRMSG)
      STOP 'NORMAL - Error'
 9000 FORMAT('NUCLEAR CHARGE =',F10.5)
 9010 FORMAT('ELECTRONIC CHARGE =',F10.5)
      end
