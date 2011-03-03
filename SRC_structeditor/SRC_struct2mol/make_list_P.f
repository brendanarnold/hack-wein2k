subroutine make_list_P(res,add,zero,molf,translate)

   USE struct
   USE atom_list, only : br

   REAL*8        :: res,translate(3)
   INTEGER       :: i,j,k,ind, add
   REAL*8        :: movex,movey,movez,moveo
   LOGICAL       :: tmovex,tmovey,tmovez
   REAL*8        :: gamma0,cosg1
   REAL*8        :: pi,a2pi,a1,a2,a3,zero
   CHARACTER     :: molf*3

   pi=4.0*atan(1.0d0)
   a2pi=pi/180.0d0
   
   a1=alpha(1)*a2pi
   a2=alpha(2)*a2pi
   a3=alpha(3)*a2pi

   aa=aa*res*b2a
   bb=bb*res*b2a
   cc=cc*res*b2a

   cosg1=(cos(a3)-cos(a1)*cos(a2))/sin(a1)/sin(a2)
   gamma0=acos(cosg1)
   
   BR(1,1)=aa*1.0d0*sin(gamma0)*sin(a2)
   BR(1,2)=aa*1.0d0*cos(gamma0)*sin(a2)
   BR(1,3)=aa*1.0d0*cos(a2)    
   BR(2,1)=0.0d0              
   BR(2,2)=bb*1.0d0*sin(a1)
   BR(2,3)=bb*1.0d0*cos(a1)
   BR(3,1)=0.0d0              
   BR(3,2)=0.0d0              
   BR(3,3)=cc*1.0d0              

   write(6,*) 'lattice vectors'
   write(6,*) '          x          y          z'
   write(6,'(a,3f10.5)') ' a    ',(br(1,j),j=1,3)
   write(6,'(a,3f10.5)') ' b    ',(br(2,j),j=1,3)
   write(6,'(a,3f10.5)') ' c    ',(br(3,j),j=1,3)
   write(6,*)
   write(6,'(a,3f10.5)') 'aa, bb,cc =',aa,bb,cc
   write(6,'(a,3f10.5)') 'alphas    =',(alpha(i),i=1,3)
   write(6,*)
   
   ind=0
   DO i=1,nat
      DO j=1,mult(i)
         ind=ind+1
         pos(1:3,ind)=pos(1:3,ind)+translate(1:3)
         do k=1,3
            if (pos(k,ind).lt.0.0d0) pos(k,ind)= pos(k,ind)+1.0d0
            if (pos(k,ind).gt.1.0d0) pos(k,ind)= pos(k,ind)-1.0d0
         enddo
         call addatom2list(pos(1,ind),i,ind,0.0d0,0.0d0,0.0d0,molf)
         if (add.gt.0) then
            movex=0.0d0 
            movey=0.0d0 
            movez=0.0d0
            moveo=1.0d0!0.999d0 
            if (molf.eq.'pov') moveo=1.0d0
            if (pos(1,ind).lt.zero) movex=moveo
            if (pos(2,ind).lt.zero) movey=moveo
            if (pos(3,ind).lt.zero) movez=moveo
            if (abs(pos(1,ind)-1.0d0).lt.zero) movex=-moveo
            if (abs(pos(2,ind)-1.0d0).lt.zero) movey=-moveo
            if (abs(pos(3,ind)-1.0d0).lt.zero) movez=-moveo
            tmovex=.false.
            tmovey=.false.
            tmovez=.false. 
            if (pos(1,ind).lt.zero) tmovex=.true.
            if (pos(2,ind).lt.zero) tmovey=.true.
            if (pos(3,ind).lt.zero) tmovez=.true.
            if (abs(pos(1,ind)-1.0d0).lt.zero) tmovex=.true.
            if (abs(pos(2,ind)-1.0d0).lt.zero) tmovey=.true.
            if (abs(pos(3,ind)-1.0d0).lt.zero) tmovez=.true.
            
            if (tmovex)&
                 call addatom2list(pos(1,ind),i,ind,movex,0.0d0,0.0d0,molf)
            
            if (tmovey)&
                 call addatom2list(pos(1,ind),i,ind,0.0d0,movey,0.0d0,molf)
            
            if (tmovez)&
                 call addatom2list(pos(1,ind),i,ind,0.0d0,0.0d0,movez,molf)
            
            if (tmovex.and.tmovey)&
                 call addatom2list(pos(1,ind),i,ind,movex,movey,0.0d0,molf)
            
            if (tmovex.and.tmovez)&
                 call addatom2list(pos(1,ind),i,ind,movex,0.0d0,movez,molf)
            
            if (tmovey.and.tmovez)&
                 call addatom2list(pos(1,ind),i,ind,0.0d0,movey,movez,molf)
            
            if (tmovex.and.tmovey.and.tmovez)&
                 call addatom2list(pos(1,ind),i,ind,movex,movey,movez,molf)
            
         endif
      ENDDO
   ENDDO
   
 end subroutine make_list_P
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
 
