      subroutine info(nume,ikf,nef,ief)
!     last changes 10.07.00 ba
!                  01.11.00 ub (more flexible band selection)
!                  13.11.00 ub (test printing removed)
!
!     subroutine provides  
!      the k-points and bands indices selected
!
!     output:
!     ikf(k)   -- true if the k-th k-point is selected
!     nef(k)   -- the number of selected bands at the k-th k-point
!     ief(n,k) -- true if the n-th band at the k-th k-point is selected
!
      implicit none
      INCLUDE 'param.inc'
      integer     nume,i, ii, j, jj, k, kmaxf, ik0, ik1, ie0, ie1, j0, j1, &
                  ierr, k1, ik(nkpt), nmaxf,  nef(nkpt)
      integer,allocatable :: ie(:)
      logical     ikf(nkpt), ief(nume,nkpt), global

      allocate (ie(nume))

!     write header
!      write(6,1000)

!     find number of inequiv. atoms (from case.struct)
!      read(20,'(/A4,24X,I2)') lattic, nat
!      write(6,1010) lattic, nat

!     initialization of the selection masks
      do 10 j = 1, nkpt
         ikf(j) = .false.
         do 15 i = 1, nume
            ief(i,j) = .false.
   15    continue
   10 continue

!     which energy bands at which k-points to be filtered ?
      read(5,*,iostat=ierr) kmaxf, (ik(j),j=1,MIN(kmaxf,nkpt))
      if (ierr.lt.0) stop 'ERROR end-of-file encountered in unit 5'
      if (ierr.eq.0) then ! global band selection mode
         global = .true.
         write(6,1020) (ik(j),j=1,MIN(kmaxf,nkpt))
         k1 = 1
      else ! individual band selection mode
         global = .false.
         write(6,1025)
         k1 = kmaxf
      endif
      if (kmaxf.le.0.or.kmaxf.gt.nkpt)  &
         stop 'ERROR # k-points in case.inf'
!     loop over all band selection blocks
      do 20 k = 1, k1 ! (only one block in global mode)
         if (global) then
            read(5,*) nmaxf, (ie(i),i=1,MIN(nmaxf,nume))
            write(6,1030) (ie(i),i=1,MIN(nmaxf,nume))
            j0 = 1
            j1 = kmaxf
         else
            read(5,*) ik(k), nmaxf, (ie(i),i=1,MIN(nmaxf,nume))
            write(6,1035) ik(k), (ie(i),i=1,MIN(nmaxf,nume))
            if (ik(k).lt.0) &
               stop 'ERROR individual k-point in case.inf'
            j0 = k
            j1 = k
         end if
         if (nmaxf.le.0.or.nmaxf.gt.nume)  &
            stop 'ERROR # bands in case.inf'
!        loop over all k-point items within a band selection block
         ik0 = -1
         do 30 j = j0, j1 ! (only the current k-point in individual mode)
            if (ik(j).eq.0.or.ABS(ik(j)).gt.nkpt)  &
               stop 'ERROR k-point in case.inf'
            if (ik(j).lt.0) then 
               if (ik0.lt.0 .or. -ik(j).le.ik0) &
                  stop 'ERROR k-point range in case.inf'
               ik0 = ik0 + 1
               ik1 = -ik(j)
            else
               ik0 = ik(j)
               ik1 = ik(j)
            end if
!           loop over all k-points within a k-point range
            do 35 jj = ik0, ik1
               ikf(jj) = .true.
!              loop over all band index items of a band selection block
               ie0 = -1
               do 40 i = 1, nmaxf
                  if (ie(i).eq.0.or.ABS(ie(i)).gt.nume)  &
                     stop 'ERROR band index in case.inf'
                  if (ie(i).lt.0) then
                     if (ie0.lt.0 .or. -ie(i).le.ie0) &
                        stop 'ERROR band index range in case.inf'
                     ie0 = ie0 + 1
                     ie1 = -ie(i)
                  else
                     ie0 = ie(i)
                     ie1 = ie(i)
                  end if
!                 loop over all bands within a band index range
                  do 45 ii = ie0, ie1
                    ief(ii,jj) = .true.
   45             continue
                  ie0 = ie(i)
   40          continue
   35       continue
            ik0 = ik(j) 
   30    continue
   20 continue

!     post-processing of the selection masks
      kmaxf=0
      do 50 j = 1, nkpt
         if (ikf(j)) then
            nmaxf = 0         
            do 55 i = 1, nume
               if (ief(i,j)) nmaxf = nmaxf + 1
   55       continue
            nef(j) = nmaxf
            ikf(j) = nmaxf .gt. 0
         end if
         kmaxf = kmaxf + nef(j)
   50 continue
      write(6,1040) kmaxf

      close(5)
      close(20)

! 1000 format(' filter "case.vector"'/ &
!             ' ~~~~~~~~~~~~~~~~~~~~')
! 1010 format(' lattice = ',A4,'   nat =',I3/)
 1020 format(' global band selection'/ &
             ' special k-points:',15I4:/(18X,15I4))
 1025 format(' individual band selection')
 1030 format(' special bands   :',15I4:/(18X,15I4))
 1035 format(' k-point',I4,' -- selected bands:',12I4:/(31X,12I4)) 
 1040 format(' total number of selected states:',I7)

      return
      end
