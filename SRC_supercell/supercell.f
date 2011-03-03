      program supercell

!     ########################################      
!
!     Program generates supercell from a WIEN struct file.
!     (c) 2002 Alexander Riss (riss@m3plus.com)
!
!     ########################################


     
      IMPLICIT NONE

      character*60    fname
      character*70    fname2
      character*100   prob
      integer         numx,numy,numz,iatom,i,j,k,NATtar,nattar0
      real*8          vac(3),vacfac(3),xshift,yshift,zshift

      integer         maxnum
      parameter       (maxnum=500)
      
      character*80    TITLE
      character*4     LATTIC, lattar, IREL,CFORM, UNIT
      character*1     toplayer(3)
      real*8          AA,BB,CC,ALPHA(3),ROTLOC(3,3)
      real*8          POS(3,maxnum),oldpos(3,maxnum),bakpos(3,maxnum)
      real*8          RMT,R0,ZZ
      integer         NAT,IATNR(maxnum),MULT,ISPLIT,JRI
      character*10    ANAME
      integer         ATNR

      
      toplayer='N'      
!     ###
!     ### program input

      write(*,*) 'Program generates supercell from a WIEN struct file.'
      write(*,*)
      write (*,*) 'Filename of struct file: '
      read (*,17) fname
   17 format(A60)

      open (20,FILE=fname,STATUS='old',FORM='formatted',ERR=500)
      read(20,1000,ERR=505) TITLE
      read(20,1010,ERR=505) LATTIC,NAT,CFORM,IREL,UNIT
      cform='    '
      read(20,1020,ERR=505) AA,BB,CC,ALPHA(1),ALPHA(2),ALPHA(3)

      if (LATTIC.ne.'P' .and. LATTIC.ne.'B' .and.LATTIC.ne.'F' .and. & 
          LATTIC.ne.'H'.and. LATTIC.ne.'CXY'.and. LATTIC.ne.'CXZ') goto 510

      write (*,*)
      write (*,*) 'Number of cells in x direction: '
      read (*,*) numx
      write (*,*) 'Number of cells in y direction: '
      read (*,*) numy
      write (*,*) 'Number of cells in z direction: '
      read (*,*) numz

          write(*,*) 'Optional shift all atoms by the same amount (fractional coordinates). '
          write(*,*) 'Please enter x shift: '
          read(*,*) xshift 
          write(*,*) 'Please enter y shift: '
          read(*,*) yshift 
          write(*,*) 'Please enter z shift: '
          read(*,*) zshift 
      
!     ### read in target lattice type
      write (*,*)
      write (*,*) 'Current structure has lattice type ' // LATTIC
  101 if (LATTIC.eq.'H') then
        lattar=LATTIC
        if(numx.ne.numy) lattar='P'
      else if (LATTIC.eq.'CXY') then
        lattar='P'
      else if (LATTIC.eq.'CXZ') then
        lattar='P'
      else if((numx.eq.numy).and.(numx.eq.numz).and.(numx/2e0.eq.numx/2)) then
        write (*,*) 'Enter your target lattice type: (P,B,F)'
	      read (*,'(A4)') lattar
        if (lattar.eq.'') lattar=LATTIC
        if (lattar.eq.'p') lattar='P'
        if (lattar.eq.'f') lattar='F'
        if (lattar.eq.'b') lattar='B'
        if ((lattar.ne.'P').and.(lattar.ne.'B').and. (lattar.ne.'F')) goto 101
      else if ((numx.eq.numy).and.(numx.eq.numz)) then
        write (*,*) 'Enter your target lattice type: (P,B)'
        read (*,'(A4)') lattar
        if (lattar.eq.'') lattar=LATTIC
        if (lattar.eq.'p') lattar='P'
        if (lattar.eq.'b') lattar='B'
        if ((lattar.ne.'P').and.(lattar.ne.'B')) lattar='P'
!        if ((lattar.ne.'P').and.(lattar.ne.'B')) goto 101
      else
        lattar='P'
      end if
      write (*,*) 'Target lattice type will be ' // lattar


!     ### user may add vacuum  in x-, y- or z-direction
      if (lattar.eq.'P') then
        write(*,*)
 102    write (*,*) 'Add vacuum in x-direction for surface-slab [bohr]:'
        read (*,*,err=102) vac(1)
        if(vac(1).gt.0.00001d0) then
          write(*,*) 'Repeat atoms at x=0 at the top (y/n)'
          read(*,'(a)') toplayer(1)
          if(toplayer(1).eq.'y') toplayer(1)='Y'
        endif
 103    write (*,*) 'Add vacuum in y-direction for surface-slab [bohr]:'
        read (*,*,err=103) vac(2)
        if(vac(2).gt.0.00001d0) then
          write(*,*) 'Repeat atoms at y=0 at the top (y/n)'
          read(*,'(a)') toplayer(2)
          if(toplayer(2).eq.'y') toplayer(2)='Y'
        endif
 104    write (*,*) 'Add vacuum in z-direction for surface slab [bohr]:'
        read (*,*,err=104) vac(3)
        if(vac(3).gt.0.00001d0) then
          write(*,*) 'Repeat atoms at z=0 at the top (y/n)'
          read(*,'(a)') toplayer(3)
          if(toplayer(3).eq.'y') toplayer(3)='Y'
        endif
        vacfac(1)=AA*numx/(vac(1)+AA*numx)
        vacfac(2)=BB*numy/(vac(2)+BB*numy)
        vacfac(3)=CC*numz/(vac(3)+CC*numz)
      else if (lattar.eq.'H') then
        write(*,*)
 204    write (*,*) 'Add vacuum in z-direction for surface slab [bohr]:'
        read (*,*,err=204) vac(3)
        if(vac(3).gt.0.00001d0) then
          write(*,*) 'Repeat atoms at z=0 at the top (y/n)'
          read(*,'(a)') toplayer(3)
          if(toplayer(3).eq.'y') toplayer(3)='Y'
        endif
        vacfac(1)=1.d0
        vacfac(2)=1.d0
        vacfac(3)=CC*numz/(vac(3)+CC*numz)
      else
        vacfac(1)=1.d0
        vacfac(2)=1.d0
        vacfac(3)=1.d0
      end if

!     ### calculation of number of atoms in target file
      NATtar=NAT*numx*numy*numz
      if (LATTIC.eq.'P' .and. lattar.eq.'B') NATtar=NATtar/2     
      if (LATTIC.eq.'P' .and. lattar.eq.'F') NATtar=NATtar/4     
      if (LATTIC.eq.'B' .and. lattar.eq.'P') NATtar=NATtar*2     
      if (LATTIC.eq.'B' .and. lattar.eq.'F') NATtar=NATtar/2     
      if (LATTIC.eq.'F' .and. lattar.eq.'P') NATtar=NATtar*4     
      if (LATTIC.eq.'F' .and. lattar.eq.'B') NATtar=NATtar*2     
      if (LATTIC.eq.'CXY' .and. lattar.eq.'P') NATtar=NATtar*2     
      if (LATTIC.eq.'CXZ' .and. lattar.eq.'P') NATtar=NATtar*2     
      
      i=index(fname,'.')
      if (i.eq.0) i=1
      j=len(fname)
      fname2=fname(1:i-1) // '_super' // fname(i:j)
      
      AA=(AA*numx)+vac(1)
      BB=(BB*numy)+vac(2)
      CC=(CC*numz)+vac(3)
      write(21,1000,ERR=525) TITLE
      write(21,2010,ERR=525) lattar,NATtar,CFORM,IREL,UNIT
      write(21,2020,ERR=525) AA,BB,CC,ALPHA(1),ALPHA(2),ALPHA(3)

!    ### read atom information, call subroutine to write information every atom
      ATNR=0
      do i=1,NAT
        iatom=1
        read(20,1030,ERR=505) IATNR(iatom),(POS(k,iatom),k=1,3),MULT,ISPLIT
        do j=1,MULT-1
          iatom=iatom+1
          if (iatom.gt.maxnum) goto 530
          read(20,1031,ERR=505) IATNR(iatom),(POS(k,iatom),k=1,3)
        end do
        read(20,1050,ERR=505) ANAME,JRI,R0,RMT,ZZ 	
        aname(3:10)='        '
        read(20,1051,ERR=505) ((ROTLOC(k,j),k=1,3),j=1,3)
        call startcell(fname2,maxnum,numx,numy,numz,iatom,NAT, &
          ROTLOC,ISPLIT,MULT,IATNR,POS,oldpos,bakpos,vacfac, &
          ANAME,JRI,R0,RMT,ZZ,ATNR,LATTIC,lattar,toplayer)
      end do

!     correct for actual number of atoms, reread new struct file
      rewind 21
      read(21,1000,ERR=525) TITLE
      read(21,1010,ERR=525) lattar,NATtar,CFORM,IREL,UNIT
      read(21,2020,ERR=525) AA,BB,CC,ALPHA(1),ALPHA(2),ALPHA(3)
      do i=1,1000000
         iatom=1
         read(21,1030,end=600) IATNR(iatom),(POS(k,iatom),k=1,3),MULT,ISPLIT
        do j=1,MULT-1
          iatom=iatom+1
          if (iatom.gt.maxnum) goto 530
          read(21,1031,ERR=505) IATNR(iatom),(POS(k,iatom),k=1,3)
        end do
        read(21,1050,ERR=505) ANAME,JRI,R0,RMT,ZZ 	
        read(21,1051,ERR=505) ((ROTLOC(k,j),k=1,3),j=1,3)
      enddo
 600  nattar=i-1
      rewind 21
      open (22,FILE=fname2,STATUS='unknown',FORM='formatted',ERR=520)
      read(21,1000,ERR=525) TITLE
      read(21,1010,ERR=525) lattar,NATtar0,CFORM,IREL,UNIT
      read(21,2020,ERR=525) AA,BB,CC,ALPHA(1),ALPHA(2),ALPHA(3)
      write(22,1000) TITLE
      write(22,2010) lattar,NATtar,CFORM,IREL,UNIT
      write(22,2020) AA,BB,CC,ALPHA(1),ALPHA(2),ALPHA(3)
      do i=1,nattar
         iatom=1
         read(21,1030,end=600) IATNR(iatom),(POS(k,iatom),k=1,3),MULT,ISPLIT
         write(22,2030) IATNR(iatom),mod(POS(1,iatom)+xshift+1.d0,1.d0),mod(POS(2,iatom)+yshift+1.d0,1.d0),&
            mod(POS(3,iatom)+zshift+1.d0,1.d0),MULT,ISPLIT
        do j=1,MULT-1
          iatom=iatom+1
          if (iatom.gt.maxnum) goto 530
          read(21,1031,ERR=505) IATNR(iatom),(POS(k,iatom),k=1,3)
          write(22,2031) IATNR(iatom),mod(POS(1,iatom)+xshift+1.d0,1.d0),mod(POS(2,iatom)+yshift+1.d0,1.d0),&
               mod(POS(3,iatom)+zshift+1.d0,1.d0)
        end do
        read(21,1050,ERR=505) ANAME,JRI,R0,RMT,ZZ 	
        read(21,1051,ERR=505) ((ROTLOC(k,j),k=1,3),j=1,3)
        write(22,2050) ANAME,JRI,R0,RMT,ZZ 	
        write(22,2051) ((ROTLOC(k,j),k=1,3),j=1,3)
      enddo
      
!     ### set symmetry operations to zero
      write(22,2060,ERR=525) 0

      write(*,*)
      write(*,*) 'Supercell generated sucessfully.'
      write(*,*) 'Stored in struct file: ' // fname2
      write(*,*) 'You may need to replace an atom by an impurity or', &
                 ' distort the positions, ....'
      
      stop


!     ### errors leading to termination of program
  500 prob='Error accessing file: ' // fname
      call faterr(prob)
      stop
  505 prob='Error reading file: ' // fname
      call faterr(prob)
      stop
  510 prob='Unknown lattice type: ' // LATTIC
      call faterr(prob)
      stop
  520 prob='Error creating file: ' // fname2
      call faterr(prob)
      stop
  525 prob='Error writing to file: ' // fname2
      call faterr(prob)
      stop
  530 prob='Number of Atoms exceeds array memory.'
      call faterr(prob)
      stop


 1000 FORMAT(A80)
 1010 FORMAT(A4,23X,I3,1X,A4,/,13X,A4,6X,A4)
 2010 FORMAT(A4,'LATTICE,NONEQUIV. ATOMS',I3,1X,A4,/,'MODE OF CALC=',A4,' unit=',A4)
 1020 FORMAT(6F10.6)
 2020 FORMAT(6F10.6)
 1030 FORMAT(4X,I4,4X,F10.8,3X,F10.8,3X,F10.8,/,15X,I2,17X,I2)
 1031 FORMAT(4X,I4,4X,F10.8,3X,F10.8,3X,F10.8)
 1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)
 1051 FORMAT(20X,3F10.8)
 2060 FORMAT(I4,6X,'NUMBER OF SYMMETRY OPERATIONS')
 2030 FORMAT('ATOM',I4,': X=',F10.8,' Y=',F10.8,' Z=',F10.8,/,10X,'MULT=',I2,10X,'ISPLIT=',I2)
 2031 FORMAT('ATOM',I4,': X=',F10.8,' Y=',F10.8,' Z=',F10.8)
 2050 FORMAT(A10,1X,'NPT=',I5,2X,'R0=',F10.8,1X,'RMT=',F10.4,3X,'Z:',F5.1)
 2051 FORMAT('LOCAL ROT MATRIX:',3X,3F10.7,/,20X,3F10.7,/,20X,3F10.7)
 
      end

      subroutine faterr(prob)
        character*80 prob
        write(*,*)
        write(*,*) 'Fatal Error occured:'
        write(*,*) prob
        write(*,*) 'Program terminated.'
        write(*,*)
        return
      end

      
      
!     ###      
!     ### converts integer to a string
      character*4 function intstr(i)
        integer i
        if (i.lt.10) then
          write(intstr,'(I1)') i
        else if (i.lt.100) then
          write(intstr,'(I2)') i
        else if (i.lt.1000) then
          write(intstr,'(i3)') i
        else
          write(intstr,'(i4)') i
        end if
      end
      

      
!     ###
!     ### transform body-centered or face-centered cells to primitive cells,
!     ### then call supercell for multiplication and reducing to target lattice type as necessary
      subroutine startcell(fname2,maxnum,numx,numy,numz,iatom,  &
         NAT,ROTLOC,ISPLIT,MULT,IATNR,POS,oldpos,bakpos,vacfac, &
         ANAME,JRI,R0,RMT,ZZ,ATNR,LATTIC,lattar,toplayer)
        integer       maxnum,IATNR(maxnum),MULT,ISPLIT,JRI,NAT,ATNR
        integer       i,j,iatom,nlatA,nlatB
        integer       numx,numy,numz,iname
        real*8        POS(3,maxnum),oldpos(3,maxnum),bakpos(3,maxnum)
        real*8        RMT,R0,ZZ,ROTLOC(3,3)
        real*8        vacfac(3),latAdd(3,11)
        character*10  ANAME
        character*70  fname2
        character*4   LATTIC,lattar
        character*1   toplayer(3)

!       ### get from body-cenetred or face-centered lattice type to primitive lattice type by generating appropriate atoms 
        
!       ## primitive: atom1=atom 
        latAdd(1,1)=0.d0
        latAdd(2,1)=0.d0
        latAdd(3,1)=0.d0

!       ## body-centered: atom1=atom, atom2=atom+(1/2,1/2,1/2)        
        latAdd(1,2)=0.d0
        latAdd(2,2)=0.d0
        latAdd(3,2)=0.d0
        latAdd(1,3)=0.5d0
        latAdd(2,3)=0.5d0
        latAdd(3,3)=0.5d0
        
!       ## face-centered: atom1=atom, atom2=atom+(1/2,1/2,0), atom3=atom+(1/2,0,1/2), atom4=atom+(0,1/2,1/2)
        latAdd(1,4)=0.d0
        latAdd(2,4)=0.d0
        latAdd(3,4)=0.d0
        latAdd(1,5)=0.5d0
        latAdd(2,5)=0.5d0
        latAdd(3,5)=0.d0
        latAdd(1,6)=0.5d0
        latAdd(2,6)=0.d0
        latAdd(3,6)=0.5d0
        latAdd(1,7)=0.d0
        latAdd(2,7)=0.5d0
        latAdd(3,7)=0.5d0
!       ## cxy centered: atom1=atom, atom2=atom+(1/2,1/2,0)
        latAdd(1,8)=0.d0
        latAdd(2,8)=0.d0
        latAdd(3,8)=0.d0
        latAdd(1,9)=0.5d0
        latAdd(2,9)=0.5d0
        latAdd(3,9)=0.d0
!       ## cxz centered: atom1=atom, atom2=atom+(1/2,0,1/2)
        latAdd(1,10)=0.d0
        latAdd(2,10)=0.d0
        latAdd(3,10)=0.d0
        latAdd(1,11)=0.5d0
        latAdd(2,11)=0.d0
        latAdd(3,11)=0.5d0

        if (LATTIC.eq.'B') then
          nlatA=2
          nlatB=3
        else if (LATTIC.eq.'F') then
          nlatA=4
          nlatB=7
        else if (LATTIC.eq.'CXY') then
          nlatA=8
          nlatB=9
        else if (LATTIC.eq.'CXZ') then
          nlatA=10
          nlatB=11
        else
          nlatA=1
          nlatB=1
        end if

        
!       ### save old coordinates
        do i=1,MULT
          do j=1,3
            bakpos(j,i)=POS(j,i)
          end do
        end do
        
        iname=0 

!       ### add positions and call makecell
        do i=nlatA,nlatB
          do iatom1=1,MULT
            do j=1,3
              POS(j,iatom1)=bakpos(j,iatom1)+latAdd(j,i)
              if (POS(j,iatom1).ge.1) POS(j,iatom1)=POS(j,iatom1)-1 
            end do
          end do
          call makecell(fname2,maxnum,numx,numy,numz,iatom,NAT, &
            ROTLOC,ISPLIT,MULT,IATNR,POS,oldpos,vacfac,         &
            ANAME,JRI,R0,RMT,ZZ,ATNR,LATTIC,lattar,iname,toplayer)
        end do
        
        return
      end
      

      
!     ###
!     ### make supercell
      subroutine makecell(fname2,maxnum,numx,numy,numz,iatom, &
         NAT,ROTLOC,ISPLIT,MULT,IATNR,POS,oldpos,vacfac,      &
         ANAME,JRI,R0,RMT,ZZ,ATNR,LATTIC,lattar,iname,toplayer)
        integer       maxnum,IATNR(maxnum),MULT,ISPLIT,JRI,NAT,ATNR
        integer       i,j,ix,iy,iz,iatom
        integer       laname,iname,curmul
        integer       numx,numy,numz
        real*8        POS(3,maxnum),oldpos(3,maxnum)
        real*8        posx,posy,posz,vacfac(3)
        real*8        RMT,R0,ZZ,ROTLOC(3,3)
        character*10  ANAME,aname2
        character*100 prob
        character*70  fname2
        character*4   intstr
        character*4   LATTIC,lattar
        character*1   toplayer(3)
        logical       cut           

!       ### strip whitspaces from ANAME
        laname=1
        do i=1,len(ANAME)
          if (ANAME(i:i).ne.' ') laname=i
        end do
        if (laname.gt.8) laname=8
!       ### names should not be changed to third position (otherwise they would not be interpeted correctly)	
	      if (laname.lt.3) laname=3

!       ### save old coordinates
        do i=1,iatom
          do j=1,3
            oldpos(j,i)=POS(j,i)
          end do
        end do
        
        iz1=numz-1
        iy1=numy-1
        ix1=numx-1
        if(toplayer(1).eq.'Y') ix1=numx
        if(toplayer(2).eq.'Y') iy1=numy
        if(toplayer(3).eq.'Y') iz1=numz
!       ### calculate new coordinates now
        do iz=0,numz
          do iy=0,numy
            do ix=0,numx
            
!             ### check if some atoms should be "reduced", when transforming to a face-centered or body-centered type lattice
              iatom=1
              iatom1=1
              curmul=MULT
              do i=1,MULT
                cut=.false.
        if((iz.eq.numz).and.(toplayer(3).ne.'Y'.or.oldpos(3,iatom).ne.0.d0)) cut=.true.
        if((iy.eq.numy).and.(toplayer(2).ne.'Y'.or.oldpos(2,iatom).ne.0.d0)) cut=.true.
        if((ix.eq.numx).and.(toplayer(1).ne.'Y'.or.oldpos(1,iatom).ne.0.d0)) cut=.true.
                posx=0.d0+oldpos(1,iatom)/numx+1.d0*ix/numx
                posy=0.d0+oldpos(2,iatom)/numy+1.d0*iy/numy
                posz=0.d0+oldpos(3,iatom)/numz+1.d0*iz/numz
                iatom=iatom+1
                if (lattar.eq.'B') then
                  if (posz.ge.0.5d0) cut=.true.
                else if (lattar.eq.'F') then
                  if (posz.ge.0.5d0 .or. posy.ge.0.5d0) cut=.true.
                end if
                if (cut) then
                  curmul=curmul-1
                  goto 201
                end if
                POS(1,iatom1)=posx*vacfac(1)
                POS(2,iatom1)=posy*vacfac(2)
                POS(3,iatom1)=posz*vacfac(3)
               iatom1=iatom1+1
  201         end do
              
              if (curmul.lt.1) goto 202
              
!             ### get name of atom              
              ATNR=ATNR+1
              iname=iname+1
              if (iname.gt.100) then
                if (laname.gt.6) laname=6
                    end if
              if (iname.gt.10) then
                if (laname.gt.7) laname=7
              end if
              aname2= ANAME(1:laname) // '_' // intstr(iname)

!             ## write to file
              do i=1,curmul
                if (i.eq.1) then
                  write(21,2030,ERR=625) ATNR,(POS(k,i),k=1,3),curmul,ISPLIT
                else
                  write(21,2031,ERR=625)ATNR,(POS(k,i),k=1,3)
                end if
              end do
              write(21,2050,ERR=625) aname,JRI,R0,RMT,ZZ 	
              write(21,2051,ERR=625) ((ROTLOC(k,j),k=1,3),j=1,3)

  202       end do
          end do
        end do

        return
        
!       ### has to be defined here again

 2030   FORMAT('ATOM',I4,': X=',F10.8,' Y=',F10.8,' Z=',F10.8,/,10X, &
               'MULT=',I2,10X,'ISPLIT=',I2)
 2031   FORMAT('ATOM',I4,': X=',F10.8,' Y=',F10.8,' Z=',F10.8)
 2050   FORMAT(A10,1X,'NPT=',I5,2X,'R0=',F10.8,1X,'RMT=',F10.4,3X,'Z:',F5.1)
 2051   FORMAT('LOCAL ROT MATRIX:',3X,3F10.7,/,20X,3F10.7,/,20X,3F10.7)
 
  625   prob='Error writing to file: ' // fname2
        call faterr(prob)
        stop

      end


