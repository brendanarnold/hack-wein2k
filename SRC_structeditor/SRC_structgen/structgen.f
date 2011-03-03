 program structgen

   use struct, only : write_struct, init_struct, tobohr, sgnum

   character*80 filein, line
   integer i,n,nato
   logical usespg

   nato=10000
   usespg=.true.

   call init_struct(nato)

   call getarg(1,filein)

   n=len_trim(filein)
   if (filein(n-3:n).eq.'.cif') then
      call scan_cif(filein,usespg,nato)
      filein=filein(1:n-4)
   elseif (filein(n-6:n).eq.'.octave') then
      call scan_octave(filein,nato)
   else
      call scan_in(filein,nato)
   endif

   if (usespg) then
      call spacegroup
      close(1)
      call getlattype(.false.)
   endif

   if (sgnum.eq.0) call getsgnum
   call getzz

   open(unit=20,file=trim(filein)//'.struct',status='unknown')

   call write_struct

   close(20)
   
 end program structgen
