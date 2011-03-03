 subroutine pov(logncm,add)

   USE struct
   USE povprop
   USE atom_list  
   
   INTEGER       :: i,j,ind,k
   real*8        :: fi(ndif),teta(ndif),rnic,mm(nat),mmtmp
   logical       :: logncm,add,plotbox,plotca,plotbond
   real*8        :: mmm
   real*8        :: trox,troy,troz,fmagcol,rgba_tmp(3),ttr,ttr1,ttr2

   call init_povprop
   
     
   read(3,*)
   read(3,*) plotbox
   read(3,*) (boxrgb(i),i=1,3)
   read(3,*) iboxfin
   read(3,*) boxw
   read(3,*)
   read(3,*) plotca
   read(3,*) (cargb(i),i=1,3)
   read(3,*) icafin
   read(3,*) caw,cal
   read(3,*)
   read(3,*) (supercell(i),i=1,3)
   read(3,*)
   read(3,*) iatfin
   read(3,*) ncol
   do i=1,ncol
      read(3,*) rgb_name(i),(rgb_tmp(k,i),k=1,3)  
   enddo
   read(3,*)
   read(3,*) nvdwr
   do i=1,nvdwr
      read(3,*) vdwr_name(i),vdwr_tmp(i) 
   enddo
   read(3,*)
   read(3,*) ibondcol,(rgbbond(i),i=1,3)
   read(3,*) ibondfin
   read(3,*) lgradbond, bondw
   read(3,*) nbond
   do i=1,nbond
      read(3,*) bond_name(i,1), bond_name(i,2), bond_dis(i) 
   enddo
   read(3,*)
   read(3,*) imagfin
   read(3,*) imagcol, fmagcol
   read(3,*) (rgba_tmp(k),k=1,3)
   read(3,*) minmm,maxmm,mincm
   read(3,*) wmom

   if (logncm) then
      mmm=0.0d0 
      read(4,*)      
      read(4,*) (qspir(i),i=1,3)    
      k=0
      DO i=1,nat
         read(5,*) 
         read(5,*) rnic,rnic,mmtmp
            mm(i)=mmtmp
            mmm=max(abs(mm(i)),mmm)
         do j=1,mult(i)
            k=k+1
            read(4,*) fi(k),teta(k)
         enddo            
      enddo
      read(4,*)
      do i=1,nat
         if (abs(mm(i)).lt.mincm) then
            cyll(i)=0.0d0 
            cylw(i)=0.0d0
            coneh(i)=0.0d0
            conew(i)=0.0d0
         else
            cyll(i)=(mm(i)/mmm)*maxmm
            cylw(i)=wmom
            coneh(i)=wmom*2.7d0
            conew(i)=wmom*1.7d0
         endif
         if (cyll(i).lt.0.0d0) coneh(i)=-coneh(i)
      enddo     
   endif

   write(6,*) 'plot box',plotbox
   write(6,*) 'plot axes',plotca


!!$   write(2,'(a)') ' #include "finish.inc"'
!!$   write(2,'(a)') ' #include "colors.inc"'
!!$   write(2,'(a)') ' #include "metals.inc"'
!!$   write(2,'(a)') ' #include "glass.inc"'


   call generate_rgb
   call generate_vdwr
   call generate_bonds

   do i=1,nat
      rgba(1:3,i)=rgba_tmp(1:3)
   enddo
   if (imagcol.eq.0) rgba(1:3,1:nat)=fmagcol*rgb(1:3,1:nat)

   if (nbond.ne.0) plotbond=.true. 

   call declare_finish(plotbox,plotca,plotbond,logncm)
   call declare_atoms(logncm)
   
   trox=br(1,1)*supercell(1)+br(2,1)*supercell(2)+br(3,1)*supercell(3)
   troy=br(1,2)*supercell(1)+br(2,2)*supercell(2)+br(3,2)*supercell(3)
   troz=br(1,3)*supercell(1)+br(2,3)*supercell(2)+br(3,3)*supercell(3)
   trox=-trox/2.0d0 
   troy=-troy/2.0d0 
   troz=-troz/2.0d0 

   if (plotca) call declare_axes(trox,troy,troz)
   if (plotbox) call declare_box
   if (plotbond) call addbond

   write(2,*) 
   write(2,*) 
   write(2,'(a)') '// declare unit cell'
   write(2,*) 

   write(2,'(a,a,a)') '#declare UnitCell = merge {'
   if (plotbond) &
        write(2,'(a)') '        object { BOND_NET }'
   if (plotbox) &
        write(2,'(a)') '        object { BOX }'
   if (plotca) &
        write(2,'(a)') '        object { AXES}'

   do i=1,nlist

      do i1=1,supercell(1)      
         do i2=1,supercell(2)      
            do i3=1,supercell(3)      

               call addatom(i,teta,fi,logncm,add,i1,i2,i3)

            enddo
         enddo
      enddo

   enddo

   ttr=1.0d0*sqrt(trox**2+troy**2+troz**2)
   ttr1=1.3d0*ttr
   ttr2=-4.2d0*ttr  
      

   write(2,'(a)')     '                    }'

   write(2,*) 
   write(2,*) 
   write(2,*) 
   write(2,*)   
   write(2,'(a)') '// plot unit cell'
   write(2,*) 
   write(2,'(a,2(f10.4,a))')   ' camera {location <0,',ttr1,',',ttr2,'> look_at <0,0,0> angle 40}'
   write(2,*)
   write(2,'(a)')   ' light_source {<10,10,-20>, color rgb <1,1,1> '
   write(2,'(a)')   '               area_light <3,0,3> <-1,4,1> 5 5'
   write(2,'(a)')   '              }'
   write(2,*) 
   write(2,'(a)')   ' background { color rgb <1,1,1> }'
   write(2,*) 
   write(2,'(a)')   ' object {  UnitCell '
   write(2,'(a,3(f10.4,a))')&
                    '         translate <',trox,',',troy,',',troz,'>'
   write(2,'(a)')   '         rotate -90*x'
   write(2,'(a)')   '         rotate <0 30 0>'
   write(2,'(a)')   '        }'
   
 end subroutine pov
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine generate_rgb

   USE povprop
   USE struct

   integer  :: i,j,k
   character*5 :: rgbn,atn
   logical ok

   do i=1,ncol
      rgbn=rgb_name(i)
      do j=1,nat
         ok=.true. 
         atn=aname(j)
         do k=1,5
            if (rgbn(k:k).eq.'*') goto 1
            if (rgbn(k:k).ne.atn(k:k)) then
               ok=.false.
               goto 1
            endif 
         enddo
1        continue
         if (ok) rgb(1:3,j)=rgb_tmp(1:3,i)
      enddo
   enddo
   
 end subroutine generate_rgb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine generate_vdwr

   USE povprop
   USE struct

   integer  :: i,j,k
   character*5 :: rgbn,atn
   logical ok

   do i=1,nvdwr
      rgbn=vdwr_name(i)
      do j=1,nat
         ok=.true. 
         atn=aname(j)
         do k=1,5
            if (rgbn(k:k).eq.'*') goto 1
            if (rgbn(k:k).ne.atn(k:k)) then
               ok=.false.
               goto 1
            endif 
         enddo
1        continue
         if (ok) vdwr(j)=vdwr_tmp(i)
      enddo
   enddo
   
 end subroutine generate_vdwr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine generate_bonds

   USE atom_list 
   USE povprop
   USE struct

   integer  :: i,j1,j2,k1,k2,nb
   character*5 :: bn1,bn2,atn1,atn2
   logical ok1,ok2
   real*8  x1,x2,y1,y2,z1,z2,dd

   nb=0
   do i=1,nbond
      bn1=bond_name(i,1)
      bn2=bond_name(i,2)
      do j1=1,nlist
         atn1=aname(list_iat(j1))
         ok1=.true. 
         do k1=1,5
            if (bn1(k1:k1).eq.'*') goto 1
            if (bn1(k1:k1).ne.atn1(k1:k1)) then
               ok1=.false.
               goto 1
            endif
         enddo
1        continue
         if (ok1) then
            do j2=1,nlist
               if (j1.ne.j2) then
                  atn2=aname(list_iat(j2))
                  ok2=.true. 
                  do k2=1,5
                     if (bn2(k2:k2).eq.'*') goto 2
                     if (bn2(k2:k2).ne.atn2(k2:k2)) then
                        ok2=.false.
                        goto 2
                     endif
                  enddo
2                 continue
                  if (ok2) then
                     x1=list_rpos(1,j1)
                     x2=list_rpos(1,j2)
                     y1=list_rpos(2,j1)
                     y2=list_rpos(2,j2)
                     z1=list_rpos(3,j1)
                     z2=list_rpos(3,j2)                       
                     dd=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

                     if (dd.lt.bond_dis(i)) then
                        nb=nb+1
                        if (nb.gt.maxnb) then
                           write(6,*) 'nb.gt.maxnb',nb,maxnb
                           stop 'struct2mol - Error, nb.gt.maxnb'
                        endif
                        bond(nb,1)=j1
                        bond(nb,2)=j2
                        bondd(nb)=bond_dis(i)
                        write(6,'(a,2i5,a,2i5)') 'add bond',j1,j2,'      iat',list_iat(j1),list_iat(j2)
                     endif

                  endif
               endif
            enddo
         endif
      enddo
   enddo
   
   nbond=nb

   write(6,*) 'number of bonds',nb

 end subroutine generate_bonds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine declare_finish(plotbox,plotca,plotbond,logncm)

   USE povprop

   logical       :: logncm,plotbox,plotca,plotbond
   integer       :: iform1,iform2,iform3,iform4

   write(2,*) 
   write(2,*)    


   assign 3 to iform1
   if (iatfin.eq.1) assign 1 to iform1
   if (iatfin.eq.2) assign 2 to iform1
   if (iatfin.eq.3) assign 3 to iform1

   assign 3 to iform2
   if (imagfin.eq.1) assign 1 to iform2
   if (imagfin.eq.2) assign 2 to iform2
   if (imagfin.eq.3) assign 3 to iform2

   assign 3 to iform3
   if (iboxfin.eq.1) assign 1 to iform3
   if (iboxfin.eq.2) assign 2 to iform3
   if (iboxfin.eq.3) assign 3 to iform3

   assign 3 to iform4
   if (icafin.eq.1) assign 1 to iform4
   if (icafin.eq.2) assign 2 to iform4
   if (icafin.eq.3) assign 3 to iform4

   assign 3 to iform5
   if (ibondfin.eq.1) assign 1 to iform5
   if (ibondfin.eq.2) assign 2 to iform5
   if (ibondfin.eq.3) assign 3 to iform5

   write(2,iform1) ' #declare F_atoms = finish {','}'

   if (logncm) then
      write(2,iform2) ' #declare F_mag = finish {','}'
      write(2,*) 
      write(2,*)    
   endif

   if (plotbox) then
      write(2,iform3) ' #declare F_box = finish {','}'
      write(2,*) 
      write(2,*)    
   endif

   if (plotca) then
      write(2,iform4) ' #declare F_coordsys = finish {','}'
      write(2,*)    
      write(2,*)    
   endif   

   if (plotbond) then
      write(2,iform5) ' #declare F_bond = finish {','}'
      write(2,*)    
      write(2,*)    
   endif   

3  FORMAT(a,/,&
            '                              specular 0.6',/,&
            '                              roughness 0.002',/,&
            '                              ambient 0.2',/,&
            '                              diffuse 0.5',/,&
            '                              reflection {0.05, 1.0}',/,&
            '                              phong 0.9',/,&
            '                              phong_size 3',/,&
            '                              conserve_energy',&
          a)

2  FORMAT(a,/,&
            '                              specular 0.6',/,&
!            '                              roughness 0.002',/,&
            '                              ambient 0.4',/,&
            '                              diffuse 0.6',/,&
!            '                              reflection {0.0, 0.2}',/,&
            '                              phong 0.4',/,&
            '                              phong_size 3',/,&
            '                              conserve_energy',&
          a)

!!$2  FORMAT(a,/,&
!!$            '                              specular 0.1',/,&
!!$            '                              ambient 0.6',/,&
!!$            '                              diffuse 0.5',/,&
!!$            '                              phong 0.5',/,&
!!$            '                              phong_size 1',/,&
!!$            '                              conserve_energy',&
!!$          a,)

1  FORMAT(a,/,&
            '                              ambient 0.6',/,&
            '                              diffuse 0.5',/,&
          a)

 end subroutine declare_finish


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 subroutine declare_atoms(logncm)

   USE struct
   USE povprop

   logical       :: logncm
   integer i

   write(2,'(a)') '// declare atoms'
   write(2,*) 

   do i=1,nat      
      if (logncm.and.(abs(cyll(i)).gt.minmm)) then
         write(2,'(a,a,a)') '#declare ',aname(i)(1:3),' = merge {'
         write(2,'(a,f8.4)')&
              '                   sphere { <0,0,0>,',vdwr(i)
         write(2,'(a,3(f6.3,a))') &
              '                           pigment { rgb <',rgb(1,i),',',rgb(2,i),',',rgb(3,i),'> }'
         write(2,'(a)')&
              '                           finish { F_atoms } '
         write(2,'(a)')&
              '                   }'

         write(2,'(a,f8.4,a,f8.4,a)') &
              '                   cylinder { <0,0,0> <0,0,',cyll(i),'>',cylw(i),' open '
         write(2,'(a,3(f6.3,a))') &
              '                           pigment { rgb <',rgba(1,i),',',rgba(2,i),',',rgba(3,i),'> }'
         write(2,'(a)')&
              '                           finish { F_mag } '
         write(2,'(a)')&
              '                   }'
         write(2,'(a,f8.4,a,f8.4,a,f8.4,a)') & 
              '                   cone { <0,0,',cyll(i),'> ',conew(i),' <0,0,',coneh(i)+cyll(i),'> 0 '
         write(2,'(a,3(f6.3,a))') &
              '                           pigment { rgb <',rgba(1,i),',',rgba(2,i),',',rgba(3,i),'> }'
         write(2,'(a)')&
              '                           finish { F_mag } '
         write(2,'(a)')&
              '                   }'
         write(2,'(a)')&
              '               }'
      else
         write(2,'(a,a,a)') '#declare ',aname(i)(1:3),' = sphere {'
         write(2,'(a,f8.4)')&
              '                           <0,0,0>,',vdwr(i)
         write(2,'(a,3(f6.3,a))') &
              '                           pigment { rgb <',rgb(1,i),',',rgb(2,i),',',rgb(3,i),'> }'
         write(2,'(a)')&
              '                           finish { F_atoms } '
         write(2,'(a)')&
              '                   }'

      endif

   enddo

 end subroutine declare_atoms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine declare_box

   USE atom_list, only :br  
   USE povprop 
   integer j,i,k,k1,k2
   real*8     :: trbx,trby,trbz
   integer    :: is1,is2,is3

   write(2,'(a)') '#declare BOX = merge {'
   
   do is1=0,supercell(1)-1      
      do is2=0,supercell(2)-1      
         do is3=0,supercell(3)-1      
            
            trbx=br(1,1)*is1+br(2,1)*is2+br(3,1)*is3
            trby=br(1,2)*is1+br(2,2)*is2+br(3,2)*is3
            trbz=br(1,3)*is1+br(2,3)*is2+br(3,3)*is3


            do j=1,3

               if ((is1.gt.0).and.((j.eq.2).or.(j.eq.3))) go to 1
               if ((is2.gt.0).and.((j.eq.1).or.(j.eq.3))) go to 1
               if ((is3.gt.0).and.((j.eq.1).or.(j.eq.2))) go to 1
                  
               write(2,'(a,4(f10.4,a))') &
                    '          cylinder { <0,0,0> <',br(j,1)+trbx,',',br(j,2)+trby,',',br(j,3)+trbz,'>',boxw, ' '
               write(2,'(a)')&
                    '                     finish { F_box } '
               write(2,'(a,3(f6.3,a))') &
                    '                     pigment { rgb <',boxrgb(1),',',boxrgb(2),',',boxrgb(3),'> }'
               write(2,'(a)') &
                    '          }'
1              continue

               do i=1,2
                  if (j.eq.1) then
                     if (i.eq.1) then
                        if (is3.gt.0) go to 2
                        k=2
                     else
                        if (is2.gt.0) go to 2
                        k=3
                     endif
                  else if (j.eq.2) then
                     if (i.eq.1) then
                        if (is3.gt.0) go to 2
                        k=1
                     else 
                        if (is1.gt.0) go to 2
                        k=3
                     endif
                  else if (j.eq.3) then
                     if (i.eq.1) then
                        if (is2.gt.0) go to 2
                        k=1
                     else
                        if (is1.gt.0) go to 2
                        k=2
                     endif
                  endif
                  write(2,'(a,7(f10.4,a))') &
                       '          cylinder { <0,0,0> <',br(j,1),',',br(j,2),',',br(j,3),'>',boxw,&
                       '  translate <',br(k,1)+trbx,',',br(k,2)+trby,',',br(k,3)+trbz,'> '
                  write(2,'(a)')&
                       '                     finish { F_box } '
                  write(2,'(a,3(f6.3,a))') &
                       '                     pigment { rgb <',boxrgb(1),',',boxrgb(2),',',boxrgb(3),'> }'
                  write(2,'(a)') &
                       '          }'
2                 continue
               enddo
               
               if (j.eq.1) then
                  k1=2
                  k2=3
               else if (j.eq.2) then
                  k1=1
                  k2=3
               else if (j.eq.3) then
                  k1=1
                  k2=2
               endif
               write(2,'(a,7(f10.4,a))') &
                    '          cylinder { <0,0,0> <',br(j,1),',',br(j,2),',',br(j,3),'>',boxw,&
                    ' translate <',br(k1,1)+br(k2,1)+trbx,',',br(k1,2)+br(k2,2)+trby,',',br(k1,3)+br(k2,3)+trbz,'> '         
               write(2,'(a)')&
                    '                     finish { F_box } '
               write(2,'(a,3(f6.3,a))') &
                    '                     pigment { rgb <',boxrgb(1),',',boxrgb(2),',',boxrgb(3),'> }'
               write(2,'(a)') &
                    '          }'
               
            enddo
            
         enddo
      enddo
   enddo
   
   write(2,'(a)') &
        '          }'

 end subroutine declare_box

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine declare_axes(trox,troy,troz)

   USE povprop 
   integer i
   real*8        :: trox,troy,troz


   write(2,'(a)') '#declare AXIS = merge {'
   
   write(2,'(a,2(f10.4,a))') &
        '          cylinder { <0,0,0> <',cal,',0,0> ',caw,' '
   write(2,'(a)')&
        '                     finish { F_coordsys } '
   write(2,'(a,3(f6.3,a))')  &
        '                     pigment { rgb <',cargb(1),',',cargb(2),',',cargb(3),'> }'
   write(2,'(a)') &
        '          }'
   write(2,'(a,f10.4,a,f10.4,a,f10.4,a)') & 
        '          cone { <',cal,',0,0> ',caw*2.0d0,' <',cal+4.0d0*caw,',0,0> 0 '   
   write(2,'(a)')&
        '                 finish { F_coordsys } '
   write(2,'(a,3(f6.3,a))')  &
        '                 pigment { rgb <',cargb(1),',',cargb(2),',',cargb(3),'> }'
   write(2,'(a)') &
        '          }'
   write(2,'(a)') &
        '       }'

   write(2,'(a)') '#declare AXES = merge {'
   write(2,'(a)') &
        '                       object { AXIS }'
   write(2,'(a,f10.4,a)') &
        '                       text {ttf "timrom.ttf" "X" 0.3,0 rotate 90*x translate',cal+4.5d0*caw,'*x }'
   write(2,'(a)') &
        '                       object { AXIS rotate <0,0,90>}'
   write(2,'(a,f10.4,a)') &
        '                       text {ttf "timrom.ttf" "Y" 0.3,0 rotate 90*x translate',cal+4.5d0*caw,'*y }'
   write(2,'(a)')&
        '                       object { AXIS rotate <0,-90,0>}'
   write(2,'(a,f10.4,a)') &
        '                       text {ttf "timrom.ttf" "Z" 0.3,0 rotate 90*x translate',cal+4.5d0*caw,'*z }'
   write(2,'(a,3(f10.4,a))')   &
        '                       translate <',-trox,',',-troy,',',-troz,'>'
   write(2,'(a)') &
        '             }'

 end subroutine declare_axes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine addatom(i,teta,fi,logncm,add,i1,i2,i3)
   
   USE struct
   USE povprop 
   USE atom_list  

   integer  :: i,iat,ind,i1,i2,i3
   real*8   :: teta(ndif),fi(ndif)
   logical  :: logncm,add
   real*8   :: trx,try,trz,zero,pi,qq
   integer  :: movea,moveb,movec

   pi=2.0d0*asin(1.0d0)
   zero=1.0d-5

   if (add) then
      if ((i1.gt.1).and.(abs(list_pos(1,i)).lt.zero)) return
      if ((i2.gt.1).and.(abs(list_pos(2,i)).lt.zero)) return
      if ((i3.gt.1).and.(abs(list_pos(3,i)).lt.zero)) return
   endif

   iat = list_iat(i)
   ind = list_ind(i)

   trx=list_rpos(1,i)+(i1-1)*br(1,1)+(i2-1)*br(2,1)+(i3-1)*br(3,1)
   try=list_rpos(2,i)+(i1-1)*br(1,2)+(i2-1)*br(2,2)+(i3-1)*br(3,2)
   trz=list_rpos(3,i)+(i1-1)*br(1,3)+(i2-1)*br(2,3)+(i3-1)*br(3,3)

   movea=nint(list_pos(1,i)-pos(1,ind))
   moveb=nint(list_pos(2,i)-pos(2,ind))
   movec=nint(list_pos(3,i)-pos(3,ind))

   qq=2.0d0*pi*(qspir(1)*(movea+i1-1)+qspir(2)*(moveb+i2-1)+qspir(3)*(movec+i3-1))
   qq=qq*180.0d0/pi

   write(2,'(a,a)')       '        object {',aname(iat)(1:3)
   if (logncm) then
      write(2,'(a,f6.1,a)')  '                 rotate <0,',teta(ind),',0>'
      write(2,'(a,f6.1,a)')  '                 rotate <0,0,',fi(ind)+qq,'>'
   endif
   write(2,'(a,3(f10.4,a))') '                 translate <',trx,',',try,',',trz,'>'
   write(2,'(a)')         '               }'

   write(6,'(a,2i5)') ' attom added to pov',iat,ind

 end subroutine addatom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine addbond

   USE atom_list 
   USE povprop
   integer    :: ibond
   integer    :: i1,i2,ii1,ii2
   real*8     :: x1,x2,y1,y2,z1,z2   
   real*8     :: dd,dx,dy,dz,r1,r2,p1,p2
   real*8     :: tribx,triby,tribz
   real*8     :: trjbx,trjby,trjbz
   integer    :: is1,is2,is3

   write(2,*)
   write(2,'(a)') '#declare  BOND_NET = merge {'

   do is1=0,supercell(1)-1      
      do is2=0,supercell(2)-1      
         do is3=0,supercell(3)-1      
            
            tribx=br(1,1)*is1+br(2,1)*is2+br(3,1)*is3
            triby=br(1,2)*is1+br(2,2)*is2+br(3,2)*is3
            tribz=br(1,3)*is1+br(2,3)*is2+br(3,3)*is3


            do js1=0,supercell(1)-1      
               do js2=0,supercell(2)-1      
                  do js3=0,supercell(3)-1      
                     
                     trjbx=br(1,1)*js1+br(2,1)*js2+br(3,1)*js3
                     trjby=br(1,2)*js1+br(2,2)*js2+br(3,2)*js3
                     trjbz=br(1,3)*js1+br(2,3)*js2+br(3,3)*js3

                     if (abs(js1-is1).gt.1) cycle
                     if (abs(js2-is2).gt.1) cycle
                     if (abs(js3-is3).gt.1) cycle

            do ibond=1,nbond

               i1=bond(ibond,1)
               i2=bond(ibond,2)
               
               ii1=list_iat(i1)
               ii2=list_iat(i2)
               
               x1=list_rpos(1,i1)+tribx
               x2=list_rpos(1,i2)+trjbx
               y1=list_rpos(2,i1)+triby
               y2=list_rpos(2,i2)+trjby
               z1=list_rpos(3,i1)+tribz
               z2=list_rpos(3,i2)+trjbz  
               
               dx=x2-x1
               dy=y2-y1
               dz=z2-z1
               dd=sqrt(dx**2+dy**2+dz**2)
               
               if (dd.gt.bondd(ibond)) cycle

               r1=vdwr(ii1)
               r2=vdwr(ii2) 
               
               p1=dd-(dd/((r1/r2)+1.0d0))-lgradbond/2.0d0
               p2=p1+lgradbond/2.0d0
               
               p1=p1/dd
               p2=p2/dd
               
               write(2,*)
               write(2,'(2(a,2i3))') '// bond between atoms:',i1,i2, '   atoms',ii1,ii2 
               write(2,'(a,7(f10.4,a))') &
                    '                   cylinder { <0,0,0> <',dx,',',dy,',',dz,'>',bondw,' open '
               
               
               if (ibondcol.eq.0) then
                  write(2,'(a,3(f6.3,a))')  &
                       '                             pigment { rgb <',rgbbond(1),',',rgbbond(2),',',rgbbond(3),'> }'
               else
                  write(2,'(a,3(f10.4,a))') &
                       '                           pigment { gradient <',dx,',',dy,',',dz,'>'
                  
                  write(2,'(a)') &
                       '                                              color_map {'
                  write(2,'(a,3(f6.3,a))') &
                  '                                                   [ 0.0 rgb <',rgb(1,ii1),',',rgb(2,ii1),',',rgb(3,ii1),'> ]'
                  write(2,'(a,f8.4,a,3(f6.3,a))') &
                  '                                                   [',p1,' rgb <',rgb(1,ii1),',',rgb(2,ii1),',',rgb(3,ii1),'> ]'
                  write(2,'(a,f8.4,a,3(f6.3,a))') &
                  '                                                   [',p2,' rgb <',rgb(1,ii2),',',rgb(2,ii2),',',rgb(3,ii2),'> ]'
                  write(2,'(a,3(f6.3,a))') &
                  '                                                   [ 1.0 rgb <',rgb(1,ii2),',',rgb(2,ii2),',',rgb(3,ii2),'> ]'
                  write(2,'(a)')&
                  '                                               }'   
                  write(2,'(a,f10.4)') &
                       '                                      scale',dd
                  write(2,'(a)')&
                       '                            }'   
               endif
               write(2,'(a)')&
                    '                           finish { F_bond } '
               write(2,'(a,3(f10.4,a))')&
                    '                           translate <',x1+trbx,',',y1+trby,',',z1+trbz,'>'
               write(2,'(a)')&
                    '                    }'   
               
            enddo
            
            
         enddo
      enddo
   enddo   
   
         enddo
      enddo
   enddo   
   
   write(2,'(a)') '                     }'
 
 end subroutine addbond
