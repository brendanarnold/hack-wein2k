PROGRAM struct2mol

   USE struct
   USE atom_list

   CHARACTER*80  :: DEFFN, ERRFN
   CHARACTER*80  :: FNAME,VECFN  
   CHARACTER*11  :: STATUS,FORM
   CHARACTER*67  ::   ERRMSG   
   REAL*8          :: res, translate(3),zeroadd
   INTEGER       :: i,j,k,jj,ind, iunit,irecl, info, add
   CHARACTER     :: atom*10,atom2*10,molf*3
   LOGICAL       :: logncm,inexist,cfexist,addl,logfor


   CALL GTFNAM(DEFFN,ERRFN,IPROC)
   CALL ERRFLG(ERRFN,'Error in struct2mol')

   molf='   '
   logncm=.false.   
   logfor=.false.   
   cfexist=.false.

   OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)
10 CONTINUE
   READ (1,*,END=20,ERR=960) IUNIT,FNAME,STATUS,FORM,IRECL
   if (iunit.eq.5) then
       inquire (FILE=FNAME,EXIST=cfexist)
       if (.not.cfexist) then
          open(unit=5,file=FNAME)
          close(5)
       endif
    endif
   if (iunit.eq.3) then
      inquire (FILE=FNAME,EXIST=inexist)
      if (.not.inexist) then
         open(unit=3,file=FNAME)
         close(3)
      endif
   endif
   OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
   if (iunit.eq.2) then
      DO i=1,80
         IF (fname(i:i).eq.'.') THEN
            molf=fname(i+1:i+3)
            GOTO 10
         ENDIF
      ENDDO
   endif
   if (iunit.eq.4) logncm=.true.
   if (iunit.eq.8) logfor=.true.
   GOTO 10
20 CONTINUE
   CLOSE (1)

   if (logncm.and.logfor) then
      write(0,*) 'only one possible (logncm,logfor)'
      stop
   endif

   write(6,'(a,a,a,/)') 'format=>',molf,'<'
  
   CALL init_struct

   call testunit(3,inexist)
   if (logncm) call testunit(5,cfexist)
   
   if (.not.inexist) call create_in
   if (logncm.and.(.not.cfexist)) call create_cf

   do i=1,nat
      atom=aname(i)
      jj=0 
      atom2='            ' 
      do j=1,7
         if (atom(j:j).ne.' ') then 
            jj=jj+1
            atom2(jj:jj)=atom(j:j)
         endif
      enddo
      aname(i)=atom2
   enddo


   iunit=3
   read(3,*,err=960) res
   read(3,*,err=960) add,zeroadd
   read(3,*,err=960) (translate(i),i=1,3)

   addl=.true.
   if (add.eq.0) addl=.false.
   
   call init_atom_list

   write(6,'(/,a,/)') 'lattice type >',LATTIC 
   IF(LATTIC(1:1).EQ.'S'.OR.LATTIC(1:1).EQ.'P') THEN
      call make_list_P(res,add,zeroadd,molf,translate)
   else if (LATTIC(1:1).EQ.'F') then
      call make_list_F(res,add,zeroadd,molf,translate)
   else if (LATTIC(1:1).EQ.'B') then
      call make_list_B(res,add,zeroadd,molf,translate)
   else if (LATTIC(1:1).EQ.'H') then
      call make_list_H(res,add,zeroadd,molf,translate)
   else if (LATTIC(1:1).EQ.'R') then
      call make_list_R(res,add,zeroadd,molf,translate)
   else if (LATTIC(1:3).EQ.'CXY') then
      call make_list_CXY(res,add,zeroadd,molf,translate)
   else if (LATTIC(1:3).EQ.'CYZ') then
      call make_list_CYZ(res,add,zeroadd,molf,translate)
   else if (LATTIC(1:3).EQ.'CXZ') then
      call make_list_CXZ(res,add,zeroadd,molf,translate)
   else
      goto 980
   endif

   IF (molf.eq.'pdb') THEN
      
      call pdb

   ELSE IF (molf.eq.'xtl') THEN

      call xtl
 
   ELSE IF (molf.eq.'xyz') THEN
      
      call xyz
       
   ELSE IF (molf.eq.'dx') THEN
      
      call dataexp(logncm,addl,logfor)
       
   ELSE IF (molf.eq.'pov') THEN
      
      call pov(logncm,addl)
       
  ELSE
     GOTO 970
  ENDIF

  CLOSE(2)

  CALL ERRCLR(ERRFN)

  STOP 'struct2mol END'

!
!        error handling
!
910 INFO = 1
  WRITE (ERRMSG,9000) FNAME
  CALL OUTERR('struct2mol',ERRMSG)
  GOTO 999
920 INFO = 2
  WRITE (ERRMSG,9010) IUNIT
  CALL OUTERR('struct2mol',ERRMSG)
  WRITE (ERRMSG,9020) FNAME
  CALL OUTERR('struct2mol',ERRMSG)
  WRITE (ERRMSG,9030) STATUS, FORM
  CALL OUTERR('struct2mol',ERRMSG)
  GOTO 999
960 INFO = 7
  WRITE (ERRMSG,9040) iunit
  CALL OUTERR('struct2mol',ERRMSG)
  GOTO 999
970 INFO=8
  WRITE (ERRMSG,9050) molf
  CALL OUTERR('struct2mol',ERRMSG)
  GOTO 999
980 INFO=9
  WRITE (ERRMSG,9060) lattic
  CALL OUTERR('struct2mol',ERRMSG)
  GOTO 999
999 STOP 'struct2mol - Error'

9000 FORMAT('can''t open definition file ',A40)
9010 FORMAT('can''t open unit: ',I2)
9020 FORMAT('       filename: ',A50)
9030 FORMAT('         status: ',A,'  form: ',A)
9040 FORMAT('Error reading unit: ',i3)
9050 FORMAT('Unknown format: ',A3)
9060 FORMAT('Unknown lattice: ',A3)

  END

   subroutine testunit(iunit,test)
     
     logical test
     integer iunit

     if (test) then
        read(iunit,*,end=1) rnic
        goto 2
1       continue
        test=.false.
     endif
2    continue
     rewind(iunit)

   end subroutine testunit




