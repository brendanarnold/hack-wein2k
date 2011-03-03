program ncmsymmetry

   use struct
   use rotations
   use Ylm_rot, only : init_Ylm_rot

   character*80  :: deffn, errfn
   character*80  :: fname,vecfn  
   character*11  :: status,form
   character*67  :: errmsg   

   ideb=1

   call gtfnam(deffn,errfn,iproc)
   call errflg(errfn,'error in ncmsymetry')

   open (1,file=deffn,status='old',err=910)
10 continue
   read (1,*,end=20,err=960) iunit,fname,status,form,irecl
   open (iunit,file=fname,status=status,form=form,err=920)
   goto 10
20 continue
   close (1)

   call init_struct
   call init_rotations

   call latgen_struct

   if (nsym.eq.0) then
      call symgen
   else         
      call order_symoper
   endif

   call test_symetry
   call init_Ylm_rot

   call make_point_groups

   if (.not.nmagmod) then
      if (fixnmag) then 
         call fix_nonmagatom_1
         call fix_nonmagatom_2
      endif
      call fix_rotloc
   endif

   call write_struct(21)
   call write_inncm

   call make_struct_klist
   call write_struct(22)

   
   call def_Ylm_rot_mat
   call def_spin_rot_mat

   call make_lm_list

   call errclr(errfn)

   stop 'ncmsymetry end'
!
!        error handling
!
910 info = 1
  write (errmsg,9000) fname
  call outerr('ncmsymetry',errmsg)
  goto 999
920 info = 2
  write (errmsg,9010) iunit
  call outerr('ncmsymetry',errmsg)
  write (errmsg,9020) fname
  call outerr('ncmsymetry',errmsg)
  write (errmsg,9030) status, form
  call outerr('ncmsymetry',errmsg)
  goto 999
960 info = 7
  write (errmsg,9040) fname
  call outerr('ncmsymetry',errmsg)
  goto 999

999 stop 'ncmsymetry - error'

9000 format('can''t open definition file ',a40)
9010 format('can''t open unit: ',i2)
9020 format('       filename: ',a50)
9030 format('         status: ',a,'  form: ',a)
9040 format('error reading file: ',a47)
9050 format('unknown format: ',a3)
9060 format('unknown lattice: ',a3)

end program ncmsymmetry
