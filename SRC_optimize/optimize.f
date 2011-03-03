      program structure

      character*79                rows
      CHARACTER*11      STATUS,FORM                                     
      CHARACTER*80      TITLE                                           
      CHARACTER*80      FNAME,master                                           
      integer     choice      
      character*80   orname
      character*1 backslash
      common /name/ orname
!----------------------------------------------
!     for new case edit the menue, 
!     and in subr.abc add dimensions and definitions for "vol" 
!     and new 'ichoi'-case formula
!----------------------------------------------
!     Menu
 1    write (*,*) '********************************************'
      write (*,*) '  GENERATES STRUCT-FILES AND optimize.job'
      write (*,*) 'PLEASE CHOOSE ONE OF THE FOLLOWING FEATURES:'
      write (*,*)
      write (*,*) '[1]  VARY VOLUME with CONSTANT RATIO A:B:C'
      write (*,*) '[2]  VARY C/A RATIO with CONSTANT VOLUME (tetr and hex lattices)'
      write (*,*) '[3]  VARY C/A RATIO with CONSTANT VOLUME and B/A (orthorh lattice)'
      write (*,*) '[4]  VARY B/A RATIO with CONSTANT VOLUME and C/A (orthorh lattice)'
      write (*,*) '[5]  VARY A and C (2D-case) (tetragonal or hexagonal lattice)'
      write (*,*) '[6]  VARY A, B and C (3D-case) (orthorhombic lattice)'
      write (*,*) '[7]  VARY A, B, C and Gamma (4D-case) (monoclinic lattice)'
      write (*,*) '[8]  VARY C/A RATIO and VOLUME (2D-case) (tetr and hex lattices)'
      write (*,*)
      write (*,*) '********************************************'
      write (*,*)
      read (*,*) choice
      if(choice.lt.1.or.choice.gt.8) goto 1
!
      iarg=iargc()
      if(iarg.ne.1) STOP 'Exactly one commandline argument must be given'
      call getarg(1,fname)
      OPEN(1,FILE=fname,STATUS='OLD',ERR=8000)
!      OPEN(1,FILE='//zeus/usr/lapw/def/nn.def',STATUS='OLD',ERR=8000)
 8003 READ(1,*,END=8001) IUNIT,FNAME,STATUS,FORM
      OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,iostat=ist, ERR=8002)
         if(iunit.eq.20) master=FNAME
      GOTO 8003
 8000 WRITE(*,*) ' ERROR IN OPENING OPTIMIZE.DEF !!!!'
      STOP 'optimize.DEF'
 8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
      WRITE(*,*) '    FILENAME: ',FNAME,'  STATUS: ',STATUS,'  FORM:',FORM,ist
      STOP 'OPEN FAILED'
 8001 CONTINUE
      do 8004 i=80,1,-1
 8004    if(fname(i:i).eq.'.') goto 8005
 8005 continue
      orname(1:i-1)=fname(1:i-1)
!
!...Test if case_initial.struct exists, use it or create it
!
      read(20,'(a)', END=10,ERR=10) status
      close (17)
      close (20)
      OPEN(17,FILE=master,STATUS='old',FORM=FORM,ERR=8002)
      goto 11
 10   continue
      close (20)
      OPEN(20,FILE=master,STATUS='unknown',FORM=FORM, ERR=8002)
 12   read (17,245,end=13,err=13) rows
      write(20,245) rows
      goto 12
 13   rewind (17)
 11   continue
!
      call abc(i-1,choice)
      write (*,*) 'Now run   optimize.job'
245   FORMAT (A79)
      end
        
        subroutine abc(ileng,ichoi)
!       ************** 
       common /name/ orname
   
       double precision   a,b,c,alpha,beta,gamma,va_ry(800),va_ry_v(800),va_ry_c(800),const       
       character*79                rows
       character*80                 orname
       character*7                 ext
       character*5                 vol(8)
       character*5                 zwname,zwname1 
       character*80                finame, f2name, joname(800)
       character*12                opname
       character*1 backslash
       integer                     i,ii,iii, numb  
      
       data ext /'.struct'/
       data vol /'_vol_','_coa_','_coa_','_boa_','_a+c_','_abc_','_mon_','_v+c_'/
       zwname1='xxxxx'
       if(ichoi.eq.8) then
          write(*,*) 'number of volumes (3-X): '
          read(*,*) numbv
          do i = 1,numbv
            write (*,*) 'PLEASE ENTER VALUE ',i, '(IN %)  '
            read (*,*) va_ry_v (i)
          enddo
          write(*,*) 'number of c/a ratios (3-X): '
          read(*,*) numbc
          do i = 1,numbc
            write (*,*) 'PLEASE ENTER VALUE ',i, '(IN %)  '
            read (*,*) va_ry_c (i)
          enddo
          numb=numbv*numbc
!          nodd=(numbc/2)*2-numbc+1
          numbc1=1
          numbv1=1
          do i=1,numb
           va_ry(i)=i
          enddo
       else if(ichoi.eq.7) then
 74       write(*,*) 'number of structures: 15, 81 (3x3x3x3), 256 (4x4x4x4)'
          read(*,*) numb
          if(numb.eq.15.or.numb.eq.81.or.numb.eq.256) goto 75
            write(*,*) 'You must enter a proper value. Do it again.'
            goto 74
 75         write(*,*) 'PLEASE enter a percentage change of a'
          read(*,*) va_ry1
          do i=1,numb
           va_ry(i)=i
          enddo
       else if(ichoi.eq.6) then
 76       write(*,*) 'number of structures: 10, 27 (3x3x3), 64 (4x4x4), 125 (5x5x5)'
          read(*,*) numb
          if(numb.eq.10.or.numb.eq.27.or.numb.eq.64.or.numb.eq.125) goto 77
            write(*,*) 'You must enter a proper value. Do it again.'
            goto 76
 77         write(*,*) 'PLEASE enter a percentage change of a'
          read(*,*) va_ry1
          do i=1,numb
           va_ry(i)=i
          enddo
       else if(ichoi.eq.5) then
 78       write(*,*) 'number of structures: 6, 9 (3x3), 16 (4x4), 25 (5x5), 36'
          read(*,*) numb
          inumb=sqrt(dble(numb))
          inumb=inumb*inumb
          if(numb.eq.6.or.numb.eq.inumb) goto 79
            write(*,*) 'You must enter a proper value. Do it again.'
            goto 78
 79         write(*,*) 'PLEASE enter a percentage change of a'
          read(*,*) va_ry1
          do i=1,numb
           va_ry(i)=i
          enddo
       else
       write (*,*) 'NUMBER OF STRUCTURE CHANGES ?'
       read (*,*) numb
      
       do i = 1,numb
         write (*,*) 'PLEASE ENTER VALUE ',i, '(IN %)  '
         read (*,*) va_ry (i)
       enddo
       endif

       do i = 1, numb
            write (zwname,'(F5.1)')va_ry(i)
            if(ichoi.eq.8) then
            write (zwname,'(F5.1)')va_ry_v(numbv1)
            write (zwname1,'(F5.1)')va_ry_c(numbc1)
!                if(numbc1.le.9) write (zwname,'(i2,"__",i1)')numbv1,numbc1
!                if(numbc1.gt.9) write (zwname,'(i2,"_",i2)')numbv1,numbc1
            endif
            finame = orname(1:ileng) //vol(ichoi) // zwname //zwname1
            do iii = 1,ileng+15
                if (finame(iii:iii).ne.' ') then
                  f2name (iii:iii) = finame(iii:iii)  
                else
                  f2name (iii:iii) = '_'
                endif  
            enddo
               ileng1=15
               if(ichoi.ne.8) ileng1=10
            joname(i) = f2name (1:ileng+ileng1)
            f2name= f2name(1:ileng+ileng1) // ext
            write (*,*) f2name 
            OPEN (21,file=f2name)
      
            do ii = 1,3
                  read (17,245) rows
                  write(21,245) rows
            enddo
        
            read (17,246) a, b, c, alpha, beta, gamma 
            if(ichoi.eq.1) then
              const = (1 + (va_ry (i)/100))**(1./3.)
              a = a*const
              b = b*const
              c = c*const            
            else if(ichoi.eq.2) then
              v = a * b * c
              valt = c/a
              a = (v/((1+va_ry(i)/100)*valt))**(1./3.)
              b = a
              c=((v/((1+va_ry(i)/100)*valt))**(1./3.))*((1+va_ry(i) &
             /100)*valt)           
            else if(ichoi.eq.3) then
              v = a * b * c
              valt = c/a
              balt = b/a
              a = (v/((1+va_ry(i)/100)*valt)/balt)**(1./3.)
              b = a*balt
              c=a*(1+va_ry(i)/100)*valt           
            else if(ichoi.eq.4) then
              v = a * b * c
              valt = c/a
              balt = b/a
              a = (v/((1+va_ry(i)/100)*balt)/valt)**(1./3.)
              b=a*(1+va_ry(i)/100)*balt           
              c = a*valt
            else if(ichoi.eq.5) then
              if(numb.eq.6) then
                if(i.lt.4) then
                a=a+a*va_ry1/100.d0*(i-2)
                elseif(i.eq.4) then
                c=c-c*va_ry1/100.d0
                elseif(i.eq.5) then
                c=c+c*va_ry1/100.d0
                elseif(i.eq.6) then
                a=a-a*va_ry1/100.d0
                c=c-c*va_ry1/100.d0
                endif
              else
                idiv=sqrt(dble(numb))
                a=a+a*va_ry1/100.d0*(mod(i-1,idiv)-idiv/2)
                c=c+c*va_ry1/100.d0*((i-1)/idiv-idiv/2)
              endif
              b=a
            else if(ichoi.eq.6) then
              if(numb.eq.10) then
                if(i.lt.4) then
                a=a+a*va_ry1/100.d0*(i-2)
                elseif(i.eq.4) then
                b=b-b*va_ry1/100.d0
                elseif(i.eq.5) then
                b=b+b*va_ry1/100.d0
                elseif(i.eq.6) then
                c=c-c*va_ry1/100.d0
                elseif(i.eq.7) then
                c=c+c*va_ry1/100.d0
                elseif(i.eq.8) then
                a=a-a*va_ry1/100.d0
                b=b-b*va_ry1/100.d0
                elseif(i.eq.9) then
                a=a-a*va_ry1/100.d0
                c=c-c*va_ry1/100.d0
                elseif(i.eq.10) then
                b=b-b*va_ry1/100.d0
                c=c-c*va_ry1/100.d0
                endif
              else
                idiv=(dble(numb+0.0001)**(1.d0/3.d0))
                a=a+a*va_ry1/100.d0*(mod(i-1,idiv)-idiv/2)
                b=b+b*va_ry1/100.d0*(mod((i-1)/idiv,idiv)-idiv/2)
                c=c+c*va_ry1/100.d0*((i-1)/idiv/idiv-idiv/2)
              endif
            else if(ichoi.eq.7) then
              if(numb.eq.15) then
                if(i.lt.4) then
                a=a+a*va_ry1/100.d0*(i-2)
                elseif(i.eq.4) then
                b=b-b*va_ry1/100.d0
                elseif(i.eq.5) then
                b=b+b*va_ry1/100.d0
                elseif(i.eq.6) then
                c=c-c*va_ry1/100.d0
                elseif(i.eq.7) then
                c=c+c*va_ry1/100.d0
                elseif(i.eq.8) then
                a=a-a*va_ry1/100.d0
                b=b-b*va_ry1/100.d0
                elseif(i.eq.9) then
                a=a-a*va_ry1/100.d0
                c=c-c*va_ry1/100.d0
                elseif(i.eq.10) then
                b=b-b*va_ry1/100.d0
                c=c-c*va_ry1/100.d0
                elseif(i.eq.11) then
                gamma=gamma-gamma*va_ry1/100.d0
                elseif(i.eq.12) then
                gamma=gamma+gamma*va_ry1/100.d0
                elseif(i.eq.13) then
                a=a-a*va_ry1/100.d0
                gamma=gamma-gamma*va_ry1/100.d0
                elseif(i.eq.14) then
                b=b-b*va_ry1/100.d0
                gamma=gamma-gamma*va_ry1/100.d0
                elseif(i.eq.15) then
                c=c-c*va_ry1/100.d0
                gamma=gamma-gamma*va_ry1/100.d0
                endif
              else
                idiv=(dble(numb+0.0001)**(1.d0/4.d0))
                a=a+a*va_ry1/100.d0*(mod(i-1,idiv)-idiv/2)
                b=b+b*va_ry1/100.d0*(mod((i-1)/idiv,idiv)-idiv/2)
                c=c+c*va_ry1/100.d0*(mod((i-1)/idiv/idiv,idiv)-idiv/2)
                gamma=gamma+gamma*va_ry1/100.d0*((i-1)/idiv/idiv/idiv-idiv/2)
              endif
            else if(ichoi.eq.8) then
              v = a * b * c
              v = v*(1+va_ry_v(numbv1)/100.d0)
              valt = c/a
              a = (v/((1+va_ry_c(numbc1)/100.d0)*valt))**(1.d0/3.d0)
              b = a
!              c=((v/((1+va_ry(i)/100.d0*numbc1)*valt))**(1.d0/3.d0))*((1+va_ry(i) &
 !            /100.d0*numbc1)*valt)

              c=v/a/a
print*,i,numbv1,numbc1,v,c/a
              numbc1=numbc1+1
              if(numbc1.gt.numbc) then
                numbc1=1
                numbv1=numbv1+1
              endif           
            endif
            write (*,246) a,b,c,gamma
            
            write (21,246)  a, b, c, alpha, beta, gamma 

            do ii=1,10000
                  read (17,245,END=244) rows
                  write (21,245) rows     
            end do
244         CONTINUE
            CLOSE (21)
            rewind (17)              
        enddo
                     
245     FORMAT (A79)
246     FORMAT (6F10.5)
!       creating job

        write (16,'(a)') '#!/bin/csh -f'
        write (16,*) '#   Modify this script according to your needs: '
        write (16,*) '#      Uncomment one of the lines ...'
        write (16,*) '#      Change run_lapw to runsp_lapw or use ' &
                               ,'different convergence criterium'
        write (16,*) '#      Change save_lapw -d XXX'
!        backslash=1H\
        backslash=achar(92)
        write (16,'("foreach i ( ",a1)') backslash
        do i=1,numb
        write (16,'(7x,a50,a1,a1)') joname(i)(1:ileng+ileng1+1),' ',backslash
        enddo
        write (16,*) ')'
        write (16,*) '    rm ', orname(1:ileng),ext,'         # NFS-bug'
        write (16,*) '    cp  $i', ext,' ', orname(1:ileng),ext
        write (16,*) ' '
        write (16,*) '# Please uncomment and adapt any of the lines below according to your needs '
        write (16,*) '#    cp  $i.clmsum ', orname(1:ileng),'.clmsum'
        write (16,*) '#    x dstart '
        write (16,*) '#    x dstart -c'
        write (16,*) '#    run_lapw -ec 0.0001 -in1new 3 -in1orig -renorm'
        write (16,*) '#    runsp_lapw -ec 0.0001'
        write (16,*) '#    min -I -j "run_lapw -I -fc 1.0 -i 40 "'
        write (16,*) '    run_lapw -ec 0.0001'
        write (16,*) '  '
        write (16,*) '    set stat = $status'
        write (16,*) '    if ($stat) then'
        write (16,*) '       echo "ERROR status in" $i'
        write (16,*) '       exit 1'
        write (16,*) '    endif'
        write (16,*) '    save_lapw  $i'
        write (16,*) '#    save_lapw  -f -d XXX $i'
        write (16,*) 'end'
        write (16,*) '   '
        
        close (16)
        return
        end
