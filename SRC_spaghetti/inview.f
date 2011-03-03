      subroutine inview (ymin,ymax,xcm,ycm,xoffs,yoffs,height, &
                      ef,lyflag,nb_min,nb_max, &
                      jatom,jatom_list,jtype,sizec,ticks,nticks)
!     **************************************************************
!
      use irr_param
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension jatom_list(*)
      character*40     sfont(7)
      common /prtft/   sfont
!-----------------------------------------------------------------------
!
!.....read paper and text configuration
      write(6,*) 'reading case.insp' 
      read(5,6010) tjunk
      if(tjunk(1:3).ne.'###' )  &
          stop 'case.insp: 1st line should start with ### '
      read(5,*) xoffs,yoffs
      write(6,*) xoffs,yoffs
      read(5,*) xcm,ycm
      read(5,*) ticks,nticks
!c      read(5,6020) ylabel
!c      read(5,6020) xlabel
      read(5,*) height, ifont
      read(5,*) dlwid, iprto(1),iprto(2) 
!
!.....read data configuration 
      read(5,6010) tjunk
      if(tjunk(1:3).ne.'###' )  &
          stop 'case.insp: 9th line should start with ### '
      read(5,*) ymin, ymax,lyflag
      read(5,*) iprtf, ef
      read(5,*) nb_min, nb_max
!
      jatom=0
 1    read(5,*,err=999,end=999) jatom_list(jatom+1), jtype, sizec
!      if(iprto(1).eq.0) jatom_list(jatom+1)=0
      jatom=jatom+1
      goto 1
 999  continue
!
! 
!.....checking input
      if(iprto(1).gt.3) stop 'error in case.insp: line  switch > 3'
      if(iprto(2).gt.5) stop 'error in case.insp: color switch > 5'
      if(iprtf   .gt.3) stop 'error in case.insp: fermi switch > 3'
!
!      lyflag=0
!      do i=1,26
!        if(ylabel(i:i+3).eq.'(eV)') lyflag=2
!        if(ylabel(i:i+3).eq.'[eV]') lyflag=2
!        if(ylabel(i:i+3).eq.' eV ') lyflag=2
!        if(ylabel(i:i+3).eq.'(Ry)') lyflag=1
!        if(ylabel(i:i+3).eq.'[Ry]') lyflag=1
!        if(ylabel(i:i+3).eq.' Ry ') lyflag=1
!      end do 
      if(lyflag.eq.1 ) then
        write(6,*) 'Energy in Ry'
      elseif(lyflag.eq.2 ) then
        write(6,*) 'Energy in eV'
      else
        stop 'error in case.insp, ylabel shall contain eV or Ry'
      endif 
!
 6010 format(a3)
 6020 format(a30)
      return
      end
