  program readwrite

   use struct
   use nstruct

   character*50 filename,mode,strname,dom
   logical domult

     call getarg(1,filename)
     call getarg(2,mode)
     call getarg(3,strname)

     if (trim(mode).eq.'read') then

        open(unit=20,file=trim(filename),status='old')

        call init_struct
        call latgen_struct

        close(20)

        call convert

        call write_octave(strname,domult)

     else if  (trim(mode).eq.'read_sym') then

         open(unit=20,file=trim(filename),status='old')
         call init_struct
         call write_octave_sym(strname)

     else

        call getarg(4,dom)

        if (trim(dom).eq.'1') then
           domult=.true.
        else
           domult=.false.
        endif   

        call read_octave(strname,domult)

        open(unit=20,file=trim(filename),status='unknown')

        call write_nstruct

        close(20)

     end if

   end program readwrite
