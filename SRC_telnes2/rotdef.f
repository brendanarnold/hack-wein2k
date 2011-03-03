!BOP
! !ROUTINE: RotDef
! !INTERFACE:
    subroutine RotDef
! !USES:
    use input, only : natom
    use struct
	use rotation_matrices,only : rotij, tauij
	use crash
! !DESCRIPTION:
!   RotDef generates the rotation-matrices rotloc(3,3,jatom) for
!   nonsymmorphic structures.  These matrices perform rotations
!   from the general coordination system into local systems with
!   specified point symmetry.
!   The matrices rotij(3,3,index) transform the position of an
!   atom to its corresponding position of an equivalent atom.
!
!   This routine is taken from the SRC\_lapw2 package of wien2k\_04.10.
!   it is supposed to produce the same output as it would do in lapw2.
!   For use in the elnes package, i have altered the used modules.
!   Note that rotdef is here called without any arguments.
!   i restrict rotdef to work only on the atom kind of which we want the spectrum.
!
!   As rotdef is called by inilpw, we are at the correct position in the struct-file to start reading sym ops.

! !REVISION HISTORY:
!   Taken from SRC\_lapw2/rotdef.f of wien2k\_04.10.
!   Adapted November 2004 (Kevin jorissen)
!EOP

	  implicit none

!  lOCAl VARiABlES
      real*8 toler,one,x,x1,y1,z1,y,z,half,zero
	  integer j,j1,j2,i1,l,index,ncount,jatom,index1,m,i


	  toler=dble(0.0001)
	  one=dble(1)
      zero=dble(0)
	  half=dble(0.500000001d0)


!-----------------------------------------------------------------------
!.....read in symmetry operations and nonprimitive translations
      read(20,10) nsym                                                  
      call make_symmat(nsym)
      do j=1,nsym                                                  
        read(20,11) ( (iz(j1,j2,j),j1=1,3),tau(j2,j), j2=1,3 ) 
      enddo         

 10   format(i4)                                                        
 11   format(3(3i2,F10.5/))                                              


!.....prepare matrices :

      index=0                                                           
      ncount=0                                                          
      do jatom=1,nat                                                 

        if (jatom.eq.natom) then

         index1=index+1                                                 
         do m=1,mult(jatom)                                          
            index=index+1                                               
            do i=1,nsym                                              
            x=zero                                                      
            do j=1,3                                                 
              x=x+iz(j,1,i)*pos(j,index1)                               
            enddo
            x=x+tau(1,i) + one                                           
            x= mod (x,one)                                              
            y=zero                                                      
            do j=1,3                                                 
              y=y+iz(j,2,i)*pos(j,index1)                               
            enddo
            y=y+tau(2,i) + one                                           
            y= mod (y,one)                                              
            z=zero                                                      
            do j=1,3                                                 
            z=z+iz(j,3,i)*pos(j,index1)
            enddo
            z=z+tau(3,i) + one                                           
            z= mod (z,one)                                              
            x1=abs(x-pos(1,index))                                      
            y1=abs(y-pos(2,index))                                      
            z1=abs(z-pos(3,index))                                      
            if(x1.lt.toler.and.y1.lt.toler.and.z1.lt.toler) then        
              ncount=ncount+1                                             
              do j=1,3                                                 
                tauij(j,m)=tau(j,i)     ! m replaces index in old version                               
                do l=1,3                                                 
                  rotij(j,l,m)=iz(j,l,i) ! m replaces index in old version
                enddo
			  enddo                                
              goto 30                                                     
            end if                                                      
!....check positions for centered lattices
            if(lattic(1:1).eq.'B') then
              x1=mod(x1+half,one)
              y1=mod(y1+half,one)
              z1=mod(z1+half,one)
              if(x1.lt.toler.and.y1.lt.toler.and.z1.lt.toler) then        
                ncount=ncount+1                                             
                do j=1,3                                                 
                  tauij(j,m)=tau(j,i)  ! m replaces index in old version
                  do l=1,3                                                 
                    rotij(j,l,m)=iz(j,l,i)  ! m replaces index in old version        
                    goto 30                                                     
			      enddo
				enddo
              end if 
            endif                                    
            if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CXY') then
              x1=mod(x1+half,one)
              y1=mod(y1+half,one)
              if(x1.lt.toler.and.y1.lt.toler.and.z1.lt.toler) then        
              ncount=ncount+1                                             
              do j=1,3                                                 
                tauij(j,m)=tau(j,i)   ! m replaces index in old version
                do l=1,3                                                 
                  rotij(j,l,m)=iz(j,l,i)  ! m replaces index in old version     
                enddo
			  enddo
              goto 30                                                     
              end if                                     
              x1=mod(x1+0.5d0,one)
              y1=mod(y1+0.5d0,one)
            endif                 
            if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CXZ') then
              x1=mod(x1+half,one)
              z1=mod(z1+half,one)
              if(x1.lt.toler.and.y1.lt.toler.and.z1.lt.toler) then        
                ncount=ncount+1                                             
                do j=1,3                                                 
                  tauij(j,m)=tau(j,i)  ! m replaces index in old version                                   
                  do l=1,3                                                 
                    rotij(j,l,m)=iz(j,l,i)   ! m replaces index in old version                                
                  enddo
			    enddo
                goto 30                                                     
              end if                                                      
              x1=mod(x1+0.5d0,one)
              z1=mod(z1+0.5d0,one)
            endif                 
            if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CYZ') then
              y1=mod(y1+half,one)
              z1=mod(z1+half,one)
              if(x1.lt.toler.and.y1.lt.toler.and.z1.lt.toler) then        
              ncount=ncount+1                                             
              do j=1,3                                                 
                tauij(j,m)=tau(j,i)  ! m replaces index in old version                                   
                do l=1,3                                                 
                  rotij(j,l,m)=iz(j,l,i)  ! m replaces index in old version                               
                enddo
			  enddo
              goto 30                                                     
              end if                                                      
            end if
            enddo                                                    

!           Error: no symmetry operation found
            goto 900
               
   30    enddo               

       else
         index=index+mult(jatom)
       endif
                                         
      enddo                                                          
      if(ncount.NE.mult(natom)) goto 910    ! statement changed Kj
      RETURN                                                            

!        Error messages

  900 call outerr('rotdef','no symmetry operation found.')
      write(errmsg,9000) jatom, index          
      call outerr('rotdef',errmsg)
      write(errmsg,9010) (pos(i1,jatom),i1=1,3) 
      call outerr('rotdef',errmsg)
      write(errmsg,9020) (pos(i1,index),i1=1,3) 
      call outerr('rotdef',errmsg)
      stop 'rotdef - Error'
  910 call outerr('rotdef','rotij not defined for all atoms of basis.')
      write (errmsg,9030) ncount
      call outerr('rotdef',errmsg)
      stop 'rotdef - Error'
!
 9000 format ('for jatom, index',i2,i2)
 9010 format ('atomposition of jatom',3F12.7)
 9020 format ('atomposition of index',3F12.7)
 9030 format ('ncount=',i2)
      end                                                               
