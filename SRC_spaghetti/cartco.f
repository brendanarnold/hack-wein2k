      subroutine cartco (v,nv,lattice,a,b,c,alpha,beta,gamma)
!
      IMPLICIT REAL*8 (A-H,O-Z)
!      INCLUDE 'param.inc'
!
      character  lattice*4
!
      dimension  v(3,*),vt(3)
      dimension br2(3,3)
!-----------------------------------------------------------------------
!
      pi=acos(-1.d0)
      if (lattice(1:1).eq.'H')  then
        do 10 j=1,nv
!..wrong old version
!            v(1,j)=2.d0/sqrt(3.d0)*2.d0*pi/a *v(1,j)
!            v(2,j)=2.d0*pi/b *v(2,j) + v(1,j)/2.d0
!  correct is
            v(1,j)=2.d0/sqrt(3.d0)*2.d0*pi/a *v(1,j)+ &
                                 v(2,j)/sqrt(3.d0)*2.d0*pi/b
            v(2,j)=2.d0*pi/b *v(2,j) 
            v(3,j)=2.d0*pi/c *v(3,j)
 10      continue
      else if (lattice(1:1).eq.'R')  then
      DO 11 j=1,nv                                                 
        VT(1)=v(1,j)                                                        
        VT(2)=v(2,j)                                                        
        VT(3)=v(3,j)                                                        
        v(1,j)=vt(1)*1.d0/SQRT(3.d0)+vt(2)*1.d0/SQRT(3.d0) &
               -2.d0/SQRT(3.d0)*vt(3)             
        v(2,j)=-vt(1)+vt(2)             
        v(3,j)=vt(1)+vt(2)+vt(3)             
        v(1,j)=2.d0*pi/a *v(1,j)
	v(2,j)=2.d0*pi/b *v(2,j)
	v(3,j)=2.d0*pi/c *v(3,j)
 11   CONTINUE                                                          
      else
!     TRICLINIC CASE
      SINBC=SIN(alpha*PI/180.0D0)
      COSAB=COS(gamma*PI/180.0D0)
      COSAC=COS(beta*PI/180.0D0)
      COSBC=COS(alpha*PI/180.0D0)
      SINBC2=SINBC*SINBC
      COSAB2=COSAB*COSAB
      COSAC2=COSAC*COSAC
      COSBC2=COSBC*COSBC
      WURZEL=sqrt(SINBC2-COSAC2-COSAB2+2*COSBC*COSAC*COSAB)
      BR2(1,1)= SINBC/WURZEL*(2.d0*pi/a)
      BR2(1,2)= (-COSAB+COSBC*COSAC)/(SINBC*WURZEL)*(2.d0*pi/b) 
      BR2(1,3)= (COSBC*COSAB-COSAC)/(SINBC*WURZEL)*(2.d0*pi/c)
      BR2(2,1)= 0.0
      BR2(2,2)= (2.d0*pi/b)/SINBC
      BR2(2,3)= -(2.d0*pi/c)*COSBC/SINBC
      BR2(3,1)= 0.0
      BR2(3,2)= 0.0
      BR2(3,3)= (2.d0*pi/c)
        write(6,*) 'lattice-type with triclinic basis assumed'
        write(6,'(3f10.5)') br2
        do 20 j=1,nv
            v1=v(1,j)
            v2=v(2,j)
            v3=v(3,j)
!....corrected by P.Blaha, 14.11.00, according to harmon.f
!            v(1,j)=BR2(1,1)*v1+BR2(2,1)*v2+BR2(3,1)*v3
!            v(2,j)=BR2(1,2)*v1+BR2(2,2)*v2+BR2(3,2)*v3
!            v(3,j)=BR2(1,3)*v1+BR2(2,3)*v2+BR2(3,3)*v3
            v(1,j)=BR2(1,1)*v1+BR2(1,2)*v2+BR2(1,3)*v3
            v(2,j)=BR2(2,1)*v1+BR2(2,2)*v2+BR2(2,3)*v3
            v(3,j)=BR2(3,1)*v1+BR2(3,2)*v2+BR2(3,3)*v3
!     write(*,'(i4,3f10.5,5x,3f10.5)') j,v1,v2,v3,v(1,j),v(2,j),v(3,j)
!       THIS IS A LATTICE-TYPE WITH RECTANGULAR BASIS VECTORS
!        write(6,*) 'lattice-type with rectangular basis assumed'
!        do 20 j=1,nv
!            v(1,j)=2.d0*pi/a *v(1,j)
!            v(2,j)=2.d0*pi/b *v(2,j)
!            v(3,j)=2.d0*pi/c *v(3,j)
 20      continue
      endif
      return
      end
