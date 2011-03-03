      subroutine normal(rhok,nkk,jspin,nat,mult,zz,qpw,qel,qmx,margn &
      ,norm,supn)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      CHARACTER*67       ERRMSG
      CHARACTER*79       MARGN
      character*5 norm
      complex*16 rhok(nkk,2)
      dimension qpw(2),qel(nat,2),zz(nat),mult(nat)
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
      z=z+supn
!      if ((qmx.gt.1.0D-5).and.(ABS(z-q).gt.0.5d0))
!     *   goto 900
      write(6,2042)  margn,z,q,z/q
      write(21,2042)  margn,z,q,z/q
 2042 format(/,':NEC',a2,':',1x,'NUCLEAR AND ELECTRONIC CHARGE ',3f10.5)
      if(abs(supn).gt.1.d-6) write(21,*)  &
                            '      Background charge of',supn
      if(norm.eq.'NO   ') then
        write(21,2043) 
 2043   format('          Charge not renormalized due to INPUT-option')
        return
      endif
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
