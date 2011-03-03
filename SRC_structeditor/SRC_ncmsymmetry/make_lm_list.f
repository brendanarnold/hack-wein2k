 subroutine make_lm_list

   USE struct
   use Ylm_rot
   use rotations, only : lmax, matom

   integer, allocatable    :: lm(:,:),lm_n(:,:,:),lmi_n(:)
   integer lmmax,isym
   real*8 :: clm_up,clm_dn,clm_up_i,clm_dn_i
   complex*16 :: clm_cr,clm_cr_i

   complex*16 :: f,fp,fac,sqrt2i,imag, sdm(2,2),srm(2,2)
   complex*16, allocatable ::  clm_y(:,:),clm_out(:,:)
   real*8     :: sqrt2,zero
   logical   ok,ok2
   logical, allocatable::   ok1(:)
   integer ii

   write(6,'(//,a)') '-----------------------------'
   write(6,'(a)')    '    Serching for lms         '
   write(6,'(a,/)')  '-----------------------------'

   zero=1.0d-6

   imag=(0.0d0,1.0d0)   
   sqrt2=sqrt(2.0d0)
   sqrt2i=sqrt2*imag

   ncom=(2*lmax+1)*lmax

   allocate (lm(2,ncom),lm_n(2,ncom,nat),lmi_n(nat),ok1(0:lmax))
   allocate ( clm_y(-lmax:lmax,1:3),clm_out(-lmax:lmax,1:3))

   lmi_n=0

   lmi=0
   do l=0,lmax
      lmi=lmi+1
      m=0
      if (lmi.gt.ncom) then
         write(6,*) 'lmi.gt.ncom:',ncom
         write(6,*) '            ',lmax,l
         stop 'error in make_lm_list: lmi.gt.ncom'  
      endif
      lm(1,lmi)=l
      lm(2,lmi)=m
      do m=1,l
         do ii=1,-1,-2            
            lmi=lmi+1
            if (lmi.gt.ncom) then
               write(6,*) 'lmi.gt.ncom:',ncom
               write(6,*) '            ',lmax,l
               stop 'error in make_lm_list: lmi.gt.ncom'  
            endif
            lm(1,lmi)=l
            if (ii.lt.0) lm(1,lmi)=-lm(1,lmi)
            lm(2,lmi)=m
         enddo
      enddo
   enddo

   lmmax=lmi

   iat=1
   do i=1,nat
      
     if (i.gt.1) iat=iat+mult(i-1)

     ok1(0:lmax)=.false.

      do lmi=1,lmmax
 
         l=abs(lm(1,lmi))
         m=lm(2,lmi)  
               
         if (m.eq.0) then
            fp=0.0d0
            f=1.0d0
         else if (mod(m,2).eq.0) then
            if (lm(1,lmi).gt.0) then
               fp=1.0d0
               f=1.0d0/sqrt2
            else
               fp=-1.0d0
               f=-imag/sqrt2               
            endif
         else
            if (lm(1,lmi).gt.0) then
               fp=-1.0d0
               f=-1.0d0/sqrt2
            else
               fp=1.0d0
               f=imag/sqrt2               
            endif
         endif
            
         clm_up=1.0d0
         clm_dn=-1.0d0
         clm_cr=(1.0d0,1.0d0)
         clm_out(-lmax:lmax,1:3)=(0.0d0,0.0d0)

         write(6,'(a,i3,a,i3,a,3i4)') 'atom=',i,'   lmi=',lmi,'   lm=',lm(1:2,lmi)
          
         do ipgo=1,npgopat(i)

!!$            write(*,*) ' ======== ipgo ======',ipgo

            clm_y(-lmax:lmax,1:3)=(0.0d0,0.0d0)
            do mm=-l,l
               fac=(Ylm_rot_mat(i,ipgo,l,mm,m)+Ylm_rot_mat(i,ipgo,l,mm,-m)*fp)*f
               clm_y(mm,1)=clm_y(mm,1)+fac*clm_up
               clm_y(mm,2)=clm_y(mm,2)+fac*clm_dn
               clm_y(mm,3)=clm_y(mm,3)+fac*clm_cr
            enddo
         
!!$            do mm=-l,l
!!$               write(*,*) '-- clm_y -------',l,mm
!!$               write(*,'(2f10.5)')  clm_y(mm,1)
!!$               write(*,'(2f10.5)')  clm_y(mm,2)
!!$               write(*,'(2f10.5)')  clm_y(mm,3)
!!$               write(*,*) '------------'
!!$            enddo

            do m1=-l,l

               mm=abs(m1)
               clm_up_i=0.0d0   
               clm_dn_i=0.0d0   
               clm_cr_i=0.0d0
               
               if (m1.eq.0) then
                  clm_up_i =clm_y(mm,1)
                  clm_dn_i =clm_y(mm,2)
                  clm_cr_i=clm_y(mm,3)        
               else if (mod(mm,2).eq.0) then                        
                  if (m1.gt.0) then
                     clm_up_i=(clm_y(mm,1)+clm_y(-mm,1))/(sqrt2)
                     clm_dn_i=(clm_y(mm,2)+clm_y(-mm,2))/(sqrt2)
                     clm_cr_i=(clm_y(mm,3)+clm_y(-mm,3))/(sqrt2)
                  else 
                     clm_up_i=(clm_y(mm,1)-clm_y(-mm,1))/(-sqrt2i)
                     clm_dn_i=(clm_y(mm,2)-clm_y(-mm,2))/(-sqrt2i)
                     clm_cr_i=(clm_y(mm,3)-clm_y(-mm,3))/(-sqrt2i)
                  endif
               else 
                  if (m1.gt.0) then
                     clm_up_i=(clm_y(mm,1)-clm_y(-mm,1))/(-sqrt2)
                     clm_dn_i=(clm_y(mm,2)-clm_y(-mm,2))/(-sqrt2)
                     clm_cr_i=(clm_y(mm,3)-clm_y(-mm,3))/(-sqrt2)
                  else
                     clm_up_i=(clm_y(mm,1)+clm_y(-mm,1))/(sqrt2i)
                     clm_dn_i=(clm_y(mm,2)+clm_y(-mm,2))/(sqrt2i)
                     clm_cr_i=(clm_y(mm,3)+clm_y(-mm,3))/(sqrt2i)
                  endif
               endif
               
               sdm(1,1)=clm_up_i
               sdm(2,2)=clm_dn_i
               sdm(1,2)=clm_cr_i
               sdm(2,1)=conjg(sdm(1,2))
               
               call apply_rotsym(sdm,i,ipgo)
               
               clm_up_i=sdm(1,1)
               clm_dn_i=sdm(2,2)
               clm_cr_i=sdm(1,2)
               
               clm_out(m1,1)=clm_out(m1,1)+clm_up_i
               clm_out(m1,2)=clm_out(m1,2)+clm_dn_i
               clm_out(m1,3)=clm_out(m1,3)+clm_cr_i

!!$            write(*,*) '-- clm_i --------',l,m,mm
!!$            write(*,'(2f10.5)')  clm_up_i
!!$            write(*,'(2f10.5)')  clm_dn_i
!!$            write(*,'(2f10.5)')  clm_cr_i
!!$            write(*,*) '------------'
               
            enddo
         enddo
         
         clm_out=clm_out/npgopat(i)         
         ok2=.true.         
         do mm=-l,l
            ok=(abs(clm_out(mm,1)).gt.zero)
            ok=ok.or.(abs(clm_out(mm,2)).gt.zero)
            ok=ok.or.(abs(clm_out(mm,3)).gt.zero)
            if (ok) then
               ok1(l)=.true.
               write(6,'(a,2i4)') '-- clm_out ---------------',l,mm
               write(6,'(2(2f10.5,3x))')  clm_out(mm,1),clm_out(mm,3)
               write(6,'(23x,2f10.5)')  clm_out(mm,2)
               write(6,'(a)') '--------------------------'
               if (lm(1,lmi).ge.0) then
                  ok2=ok2.and.(mm.eq.lm(2,lmi))
               else
                  ok2=ok2.and.(mm.eq.-lm(2,lmi))
               endif
            endif
         enddo

!!$         if (ok1) then
!!$            lmi_n(i)=lmi_n(i)+1
!!$            lm_n(1,lmi_n(i),i)=lm(1,lmi)
!!$            lm_n(2,lmi_n(i),i)=lm(2,lmi)
!!$         endif

!!$         if (ok2.and.ok1) then
!!$            lmi_n(i)=lmi_n(i)+1
!!$            lm_n(1,lmi_n(i),i)=lm(1,lmi)
!!$            lm_n(2,lmi_n(i),i)=lm(2,lmi)
!!$         else if (ok1) then
!!$            write(6,'(a,i4,)') 'lm not symmetric in point group:  lmi=',lmi,'   atom=',i
!!$            stop 'error in make lm_list:  lm not symmetric in point group'
!!$         endif
                     
      enddo

      do lmi=1,lmmax
 
         l=abs(lm(1,lmi))
         m=lm(2,lmi)  
      
         if (ok1(l)) then
            lmi_n(i)=lmi_n(i)+1
            lm_n(1,lmi_n(i),i)=lm(1,lmi)
            lm_n(2,lmi_n(i),i)=lm(2,lmi)
         endif

      enddo
   enddo

   write(6,'(//)')
   do i=1,nat
      write(6,'(a,i3,5x,100(2i3,2x):)') 'lm list for atom:',i,(lm_n(1:2,ii,i),ii=1,lmi_n(i)) 
   enddo

   do i=1,nat
      write(10,'(100i2:)') (lm_n(1:2,ii,i),ii=1,lmi_n(i)) 
   enddo
   write(10,'(a)') ' 14.          GMAX'
   write(10,'(a)') 'FILE        FILE/NOFILE  write recprlist'

   deallocate (lm,clm_y)

2032 FORMAT(a//)
1980 FORMAT(3X) 
2000 FORMAT(16X,I2//) 
2010 FORMAT(16X,I2,5X,I2/)
2021 FORMAT(3X,4E19.12)
2031 FORMAT(/)
2033 FORMAT(///)  
2060 FORMAT(/,13X,I6) 
2071 FORMAT(3X,3I5,2E19.12)  
1970 format(a)
78 FORMAT(1X,A20,A4,3F10.6,I5)    
77 FORMAT(1X,9(F9.7)) 
1990 FORMAT(3X,'ATOM NUMBER =',I2,5X,10A4)    
2001 FORMAT(3X,'NUMBER OF LM=',I2//) 
2011 FORMAT(3X,'CLM(R) FOR L=',I2,3X,'M=',I2/)
2076 FORMAT(3X,3I5,2E19.12,2e11.3)   
2061 FORMAT(1X,'NUMBER OF PW',I6) 
2051 FORMAT(3X,'TOTAL CHARGE DENSITY IN INTERSTITIAL')

 end subroutine make_lm_list

