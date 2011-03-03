      SUBROUTINE REDUC(v,iat,ipos,er)
      use sphe
      use atpos
      implicit none
      include 'param.inc'
!                                                                       
!     REDUCES VECTOR V TO EQUIVALENT SMALLEST VECTOR                    

      real*8 v,vhelp,vtest,vtmp
      real*8 r,rtest,r0

      integer iat,jatom
      integer i,j,ipos

      logical er,deb

      COMMON /DEBUG/ deb
      DIMENSION V(3),VHELP(3),VTEST(3), VTMP(3)                

      jatom=iabs(iatnr(iat))
      do  j=1,3                                                       
        vhelp(j)=v(j)                                                     
        vtmp(j)=v(j)-pos(j,iat)                                                     
      enddo
      R=1D9
      do I=1,NPOS                                                    
        rtest=0.d0
        do J=1,3                                                       
          vtest(j)=vtmp(j)-atp(j,i)
          rtest=rtest+vtest(j)*vtest(j)
        end do
        if (rtest.lt.r) then                                               
          r=rtest                                                           
          do j=1,3                                                       
            vhelp(j)=vtest(j)
            ipos=i
          end do
        end if                                                
      end do
      do j=1,3                       
        v(J)=vhelp(j)
      end do
      er=.false.
      if(r.gt.(rmt(jatom)*rmt(jatom))) er=.true.
      RETURN                                                            
      END                                                               

