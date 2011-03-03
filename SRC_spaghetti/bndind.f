      subroutine bndind(jk,nevl)
      use irr_param
      IMPLICIT REAL*8 (A-H,O-Z)
!
      logical        flg2
!      dimension      skir(3),eneik(NEVL,NKP)
!      dimension      ngrp(4,NEVL,NKP),ngde(4,NEVL,NKP)
!      dimension      ndeg(NEVL,NKP)
!      dimension      neik(NKP),ikg(NKP),iste(NEVL)
      integer,allocatable::     mgrp(:,:,:),mgde(:,:,:) !16,NEVL,2)
      integer,allocatable::      isqe(:) ! NEVL)
!                     ,iscol2(NEVL)
!      character*3    grpnam(NKP) 
!      common /bndid/ ngrp,ngde,ndeg,neik,iste,grpnam
!      common /eneir/ skir,eneik
!   
!-------------------------------------------------------------------
!
      allocate (     mgrp(16,NEVL,2),mgde(16,NEVL,2))
      allocate (     isqe(NEVL))
!
      iste = 0
      mgrp = 0
!  
!.....if k_i and k_i+1 belong to the same point group
      if(grpnam(jk).eq.grpnam(jk+1)) then
        do 42 je=1,neik(jk)
        do 42 ig=1,4
 42        mgrp(ig,je,1)=ngrp(ig,je,jk)
        do 43 je2=1,neik(jk+1)
        do 43 ig=1,4
 43        mgrp(ig,je2,2)=ngrp(ig,je2,jk+1)
      else
!.....if k_i and k_i+1 do not belong to the same point group
        call pgrpnr(grpnam(jk),  igrp1)
        call pgrpnr(grpnam(jk+1),igrp2)
!
        if(igrp1.lt.igrp2) then
          do 46 je=1,neik(jk)
          do 46 ig=1,4
 46          mgrp(ig,je,1)=ngrp(ig,je,jk)
        else
          do 47 je2=1,neik(jk+1)
          do 47 ig=1,4
 47          mgrp(ig,je2,2)=ngrp(ig,je2,jk+1)
        endif
!.......one of the point group is E
        if(igrp1.eq.1.or.igrp2.eq.1) then
          if(igrp1.eq.1) then
            write(6,*) 'k-points:',jk+1,jk, &
             ' Point group rel.: ',grpnam(jk+1)(1:3),' => ',grpnam(jk)(1:3)
            do 52 je=1,neik(jk+1)
              iste(je) = je
 52           iscol2(je) = ngrp(1,je,jk+1)-1
          else 
            write(*,*)'here',jk,grpnam(jk)(1:3)
            write(6,*) 'k-points:',jk,jk+1, &
             ' Point group rel.: ',grpnam(jk)(1:3),' => ',grpnam(jk+1)(1:3)
            do 56 je2=1,neik(jk)
              iste(je2) = je2
 56           iscol2(je2) = ngrp(1,je2,jk)-1
          endif
!
!.......set general paramenters
        else
          call seppt(jk,igrp1,igrp2,mgrp,nevl)
        endif
      endif
!
!.....find connecting lines and  
      do 50 je2=1,neik(jk+1)
 50      isqe(je2) = 0
!
      do 100 je=1, neik(jk)
        flg2=.true.
        if(mgrp(1,je,1).eq.0) flg2=.false.
        je2=1
        do while(flg2.and.je2.le.neik(jk+1))
          if(isqe(je2).eq.0.and.mgrp(1,je2,2).ne.0) then
!
            do 110 ig=1,16
            do 110 ig2=1,16
              if(mgrp(ig,je,1).eq.mgrp(ig2,je2,2).and. &
                 mgrp(ig,je,1).ne.0 ) then
                    iste(je)   = je2
                    isqe(je2)  = 1
                    iscol2(je) = mgrp(ig, je, 1)-1
                    goto 112
              endif
 110        continue
          endif
          je2=je2+1
        enddo
 112    continue
 100  continue
!
 900  continue
      deallocate (     mgrp,mgde)
      deallocate (     isqe)

      return
      end 
