 subroutine create_in

   use struct, only: nat

   integer i
   integer     nspec
   character*6, allocatable:: specname(:)
   real*8, allocatable::      specrad(:)
   real*8, allocatable::      speccr(:)
   real*8, allocatable::      specrgb(:,:)
   integer     nbo
   integer, allocatable::     indb(:,:) 
   real*8, allocatable::      bonddis(:)  
  
  allocate(specname(nat),specrad(nat),speccr(nat),specrgb(3,nat),&
       indb(2,nat*nat),bonddis(nat*nat))

  call getelem(nspec,specname,specrad,specrgb,speccr,nbo,indb,bonddis)

  write(3,'(a)') &
       ' 1.0                              # rescal ',&
       ' 1   0.1                             # add periodic',&
       ' 0.0 0.0 0.0                      #translation',&
       '# plot box (for povray)' ,&
       ' 1                                # 1/0 yes/no  ',&
       ' 0.8 0.8 0.8                      # rgb',&
       ' 3                                # finish',&
       ' 0.1                              # width  ',&
       '# plot xyz axes (for povray)',&
       ' 0                                # 1/0 yes/no',&
       ' 0 1 0                            # rgb',&
       ' 1                                # finish',&
       ' 0.1 7.0                          # width, length',&
       '# define supercell',&
       ' 1 1 1',&
       '# define finish and colors of atoms',&
       ' 3                                # finish'
  write(3,'(i5,a)')&
        nspec,'                             #  number of defs'
  do i=1,nspec
     write(3,'(a,3f7.3,a)')&
        specname(i),(specrgb(j,i),j=1,3),'     # atom name, rgb'
  enddo
  write(3,'(a)') &
       '# define radius of atoms'
  write(3,'(i5,a)')&
        nspec,'                             #  number of defs'
  do i=1,nspec
     write(3,'(a,f7.3,a)')&
       specname(i),specrad(i),'               # atom name, vdwr'
  enddo
  write(3,'(a)')&
       '# define bonds',&
       '  0 1 1 0.1                       # 1/0 - use atom colors, use this rgb',&
       '  3                               # finish',&
       '  0.1 0.15                         # lenght of gradient for atom colors, bond width'
  write(3,'(i5,a)')&
        nbo,'                             #  number of defs'
  do i=1,nbo
     write(3,'(2a,f7.3,a)')&
       specname(indb(1,i)),specname(indb(2,i)),bonddis(i),'             # atoms names, distans'
  enddo
  write(3,'(a)')&
       '# define color ... of momenta',&
       ' 3                                # finish',&
       ' 0 0.5                            # 0 - color as atom, scaled by factor, 1 defined bellow',&
       ' 0 0 0                            # rgb, valid only for 1    ',&
       ' 0.1 2.0 0.15                     # min on pict, max on pict, min consid ',&
       ' 0.2                              # width of cylinder'

       rewind(3)

 end subroutine create_in
