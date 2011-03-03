subroutine make_list_H(res,add,zero,molf,translate)

   USE struct
   USE atom_list, only : br

   REAL*8        :: res,translate(3)
   INTEGER       :: i,j,k,ind, add
   REAL*8        :: movex,movey,movez,moveo
   LOGICAL       :: tmovex,tmovey,tmovez
   REAL*8        :: gamma0,cosg1
   REAL*8        :: pi,a2pi,a1,a2,a3,zero,tr(3),pos_tmp(3)
   CHARACTER     :: molf*3

   pi=4.0*atan(1.0d0)
   a2pi=pi/180.0d0
   
   a1=alpha(1)*a2pi
   a2=alpha(2)*a2pi
   a3=alpha(3)*a2pi

   aa=aa*res*b2a
   bb=bb*res*b2a
   cc=cc*res*b2a

   BR(1,1)=aa*sqrt(3.0d0)/2.0d0
   BR(1,2)=-0.5d0*bb
   BR(1,3)=0.0d0
   BR(2,1)=0.0d0              
   BR(2,2)=bb
   BR(2,3)=0.0d0
   BR(3,1)=0.0d0              
   BR(3,2)=0.0d0              
   BR(3,3)=cc

   write(6,*) 'Br zone'
   write(6,*) '          x          y          z'
   write(6,'(a,3f10.5)') ' a    ',(br(1,j),j=1,3)
   write(6,'(a,3f10.5)') ' b    ',(br(2,j),j=1,3)
   write(6,'(a,3f10.5)') ' c    ',(br(3,j),j=1,3)
   write(6,*)
   write(6,'(a,3f10.5)') 'aa, bb,cc =',aa,bb,cc
   write(6,'(a,3f10.5)') 'alphas    =',(alpha(i),i=1,3)
   write(6,*)
   
   do is=1,1

      if (is.eq.1) then
         tr(1)=0.0d0
         tr(2)=0.0d0
         tr(3)=0.0d0
      endif

      ind=0
      DO i=1,nat
         DO j=1,mult(i)
            ind=ind+1
            pos_tmp(1:3)=pos(1:3,ind)+translate(1:3)+tr(1:3)
            do k=1,3
               if (pos_tmp(k).lt.0.0d0) pos_tmp(k)= pos_tmp(k)+1.0d0
               if (pos_tmp(k).gt.1.0d0) pos_tmp(k)= pos_tmp(k)-1.0d0
            enddo
            call addatom2list(pos_tmp,i,ind,0.0d0,0.0d0,0.0d0,molf)
            if (add.gt.0) then
               movex=0.0d0 
               movey=0.0d0 
               movez=0.0d0
               moveo=1.0d0!0.999d0 
               if (molf.eq.'pov') moveo=1.0d0
               if (pos_tmp(1).lt.zero) movex=moveo
               if (pos_tmp(2).lt.zero) movey=moveo
               if (pos_tmp(3).lt.zero) movez=moveo
               if (abs(pos_tmp(1)-1.0d0).lt.zero) movex=-moveo
               if (abs(pos_tmp(2)-1.0d0).lt.zero) movey=-moveo
               if (abs(pos_tmp(3)-1.0d0).lt.zero) movez=-moveo
               tmovex=.false.
               tmovey=.false.
               tmovez=.false. 
               if (pos_tmp(1).lt.zero) tmovex=.true.
               if (pos_tmp(2).lt.zero) tmovey=.true.
               if (pos_tmp(3).lt.zero) tmovez=.true.
               if (abs(pos_tmp(1)-1.0d0).lt.zero) tmovex=.true.
               if (abs(pos_tmp(2)-1.0d0).lt.zero) tmovey=.true.
               if (abs(pos_tmp(3)-1.0d0).lt.zero) tmovez=.true.
               
               if (tmovex)&
                    call addatom2list(pos_tmp,i,ind,movex,0.0d0,0.0d0,molf)
               
               if (tmovey)&
                    call addatom2list(pos_tmp,i,ind,0.0d0,movey,0.0d0,molf)
               
               if (tmovez)&
                    call addatom2list(pos_tmp,i,ind,0.0d0,0.0d0,movez,molf)
               
               if (tmovex.and.tmovey)&
                    call addatom2list(pos_tmp,i,ind,movex,movey,0.0d0,molf)
               
               if (tmovex.and.tmovez)&
                    call addatom2list(pos_tmp,i,ind,movex,0.0d0,movez,molf)
               
               if (tmovey.and.tmovez)&
                    call addatom2list(pos_tmp,i,ind,0.0d0,movey,movez,molf)
               
               if (tmovex.and.tmovey.and.tmovez)&
                    call addatom2list(pos_tmp,i,ind,movex,movey,movez,molf)
               
            endif
         ENDDO
      ENDDO
    
   enddo

 end subroutine make_list_H
