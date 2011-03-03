!BOP
! !ROUTINE: Rotate
! !INTERFACE:
    subroutine Rotate (vector,rotmat,rotvec)
! !INPUT/OUTPUT PARAMETERS:
!   vector  :   vector before rotation
!   rotvec  :   vector after rotation
!   rotmat  :   3*3 rotation matrix
! !DESCRIPTION:
!     Rotate performs a rotation of the vector from the general
!     Carthesian coordination system into the local one of the
!     jatom-th sphere.
!     This subroutine is only required for nonsymmorphic cases.
! !REVISION HISTORY:
!     taken 30 1 98 from  SRC\_lapw2/rotate.f:                         
!     Updated November 2004 (Kevin Jorissen)
!EOP
      implicit none
!  IN/OUTPUT
      real*8 vector(3),rotvec(3),rotmat(3,3)                         
!  LOCALS
      integer jcoord,j
	  real*8  dotpro

      do jcoord=1,3                                                  
         dotpro=dble(0)                                                   
         do J=1,3                                                    
            dotpro=dotpro + vector(J)*rotmat(jcoord,J)                  
         enddo                                                       
         rotvec(jcoord)=dotpro                                          
      enddo
	                                                            
      return                                                            
      end                                                               

