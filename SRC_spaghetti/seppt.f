      subroutine seppt(jk,igrp1,igrp2,mgrp,nevl) 
      use irr_param
      IMPLICIT REAL*8 (A-H,O-Z)

      character*1    dir1, dir2, invt
      dimension      mgrp(16,NEVL,2), ncm(-14:14,4)
      logical        flg1, flg2
      logical        C2(3,3),  C3(2,2), C4(2,2), C6(2,2)
      logical        IC2(3,3),IC3(2,2),IC4(2,2),IC6(2,2)
      common /cflg/  C2, C3, C4, C6, IC2, IC3, IC4, IC6
!                                                              
!        dimension      ngrp(4,NEVL,NKP),ngde(4,NEVL,NKP)
!        dimension      ndeg(NEVL,NKP)
!        dimension      neik(NKP),ikg(NKP),iste(NEVL)
!        character*3    grpnam(NKP) 
!        common /bndid/ ngrp,ngde,ndeg,neik,iste,grpnam
!
!        C2(1,1)     -   C2  -> C2
!          (1,2)     -   C2  -> C2`
!          (1,3)     -   C2  -> C2"
!          (2,1)     -   C2` -> C2    etc
!--------------------------------------------------------------
! 
      C2 =.false.
      C3 =.false.
      C4 =.false.
      C6 =.false.
      IC2=.false.
      IC3=.false.
      IC4=.false.
      IC6=.false.
!
!.....use compatibility relations to divide the IR of  
!     point group for k=jk1,jk2 into IRs of the smaller point group.     
!
      if(igrp1.lt.igrp2) then
        jec = 2  
        jk1 = jk+1 
        jk2 = jk
        i1  = igrp2   
        i2  = igrp1   
      else
        jec = 1
        jk1 = jk
        jk2 = jk+1
        i1  = igrp1   
        i2  = igrp2 
      endif
!
!.....find permutation of classes.
      write(6,*) 'k-points:',jk1,jk2, &
             ' Point group rel.: ',grpnam(jk1)(1:3),' => ',grpnam(jk2)(1:3)
      do 70 itab1=1,ntab4(jk1)
      if( cnam(itab1,jk1)(3:3).eq.'C') then
      do 71 itab2=1,ntab4(jk2)
        if( cnam(itab1,jk1)(3:4).eq.cnam(itab2,jk2)(3:4)  ) then 
        flg1=.true.
!
          do 75 ige1=1,nlkg(itab1,jk1)
          do 75 ige2=1,nlkg(itab2,jk2)
          if( lkg(itab1,ige1,jk1).eq.lkg(itab2,ige2,jk2) ) then
            read( cnam(itab1,jk1)(4:4), 500) icls
            read( cnam(itab1,jk1)(5:5), 501) dir1
            read( cnam(itab2,jk2)(5:5), 501) dir2
            read( cnam(itab2,jk2)(2:2), 501) invt                  

            if(invt.ne.'I') then
              if(icls.eq.2.and.dir1.eq.' '.and.dir2.eq.' ')  C2(1,1)=.true.
              if(icls.eq.2.and.dir1.eq.' '.and.dir2.eq.'`')  C2(1,2)=.true.
              if(icls.eq.2.and.dir1.eq.' '.and.dir2.eq.'"')  C2(1,3)=.true.
              if(icls.eq.2.and.dir1.eq.'`'.and.dir2.eq.' ')  C2(2,1)=.true.
              if(icls.eq.2.and.dir1.eq.'`'.and.dir2.eq.'`')  C2(2,2)=.true.
              if(icls.eq.2.and.dir1.eq.'`'.and.dir2.eq.'"')  C2(2,3)=.true.
              if(icls.eq.2.and.dir1.eq.'"'.and.dir2.eq.' ')  C2(3,1)=.true.
              if(icls.eq.2.and.dir1.eq.'"'.and.dir2.eq.'`')  C2(3,2)=.true.
              if(icls.eq.2.and.dir1.eq.'"'.and.dir2.eq.'"')  C2(3,3)=.true.
!
              if(icls.eq.3.and.dir1.eq.'+'.and.dir2.eq.'+')  C3(1,1)=.true.
              if(icls.eq.3.and.dir1.eq.'+'.and.dir2.eq.'-')  C3(1,2)=.true.
              if(icls.eq.3.and.dir1.eq.'-'.and.dir2.eq.'+')  C3(2,1)=.true.
              if(icls.eq.3.and.dir1.eq.'-'.and.dir2.eq.'-')  C3(2,2)=.true.
!
              if(icls.eq.4.and.dir1.eq.'+'.and.dir2.eq.'+')  C4(1,1)=.true.
              if(icls.eq.4.and.dir1.eq.'+'.and.dir2.eq.'-')  C4(1,2)=.true.
              if(icls.eq.4.and.dir1.eq.'-'.and.dir2.eq.'+')  C4(2,1)=.true.
              if(icls.eq.4.and.dir1.eq.'-'.and.dir2.eq.'-')  C4(2,2)=.true.
!
              if(icls.eq.6.and.dir1.eq.'+'.and.dir2.eq.'+')  C6(1,1)=.true.
              if(icls.eq.6.and.dir1.eq.'+'.and.dir2.eq.'-')  C6(1,2)=.true.
              if(icls.eq.6.and.dir1.eq.'-'.and.dir2.eq.'+')  C6(2,1)=.true.
              if(icls.eq.6.and.dir1.eq.'-'.and.dir2.eq.'-')  C6(2,2)=.true.
            else
              if(icls.eq.2.and.dir1.eq.' '.and.dir2.eq.' ') IC2(1,1)=.true.
              if(icls.eq.2.and.dir1.eq.' '.and.dir2.eq.'`') IC2(1,2)=.true.
              if(icls.eq.2.and.dir1.eq.' '.and.dir2.eq.'"') IC2(1,3)=.true.
              if(icls.eq.2.and.dir1.eq.'`'.and.dir2.eq.' ') IC2(2,1)=.true.
              if(icls.eq.2.and.dir1.eq.'`'.and.dir2.eq.'`') IC2(2,2)=.true.
              if(icls.eq.2.and.dir1.eq.'`'.and.dir2.eq.'"') IC2(2,3)=.true.
              if(icls.eq.2.and.dir1.eq.'"'.and.dir2.eq.' ') IC2(3,1)=.true.
              if(icls.eq.2.and.dir1.eq.'"'.and.dir2.eq.'`') IC2(3,2)=.true.
              if(icls.eq.2.and.dir1.eq.'"'.and.dir2.eq.'"') IC2(3,3)=.true.
!
              if(icls.eq.3.and.dir1.eq.'+'.and.dir2.eq.'+') IC3(1,1)=.true.
              if(icls.eq.3.and.dir1.eq.'+'.and.dir2.eq.'-') IC3(1,2)=.true.
              if(icls.eq.3.and.dir1.eq.'-'.and.dir2.eq.'+') IC3(2,1)=.true.
              if(icls.eq.3.and.dir1.eq.'-'.and.dir2.eq.'-') IC3(2,2)=.true.
!
              if(icls.eq.4.and.dir1.eq.'+'.and.dir2.eq.'+') IC4(1,1)=.true.
              if(icls.eq.4.and.dir1.eq.'+'.and.dir2.eq.'-') IC4(1,2)=.true.
              if(icls.eq.4.and.dir1.eq.'-'.and.dir2.eq.'+') IC4(2,1)=.true.
              if(icls.eq.4.and.dir1.eq.'-'.and.dir2.eq.'-') IC4(2,2)=.true.
!
              if(icls.eq.6.and.dir1.eq.'+'.and.dir2.eq.'+') IC6(1,1)=.true.
              if(icls.eq.6.and.dir1.eq.'+'.and.dir2.eq.'-') IC6(1,2)=.true.
              if(icls.eq.6.and.dir1.eq.'-'.and.dir2.eq.'+') IC6(2,1)=.true.
              if(icls.eq.6.and.dir1.eq.'-'.and.dir2.eq.'-') IC6(2,2)=.true.
            endif
            if(flg1) then
              if(invt.ne.'I') then
                write(6,510) icls,dir1,icls,dir2
              else
                write(6,511) icls,dir1,icls,dir2 
              endif
              if( mod(itab2,4).eq.0 ) write(6,*)
            endif
            flg1=.false.
          endif
 75       continue  
!
        endif
 71   continue
      endif
 70   continue
      write(6,*)
!
!.....get the compatibility relation
      call comprel( i1, i2, ncm )
!
!.....find the correct comp.rel. table.
      do 100 ie=1,neik(jk1)
      do 101 igg=1,16
 101     mgrp(igg,ie,jec) = 0
      igg=0
      do 110 ig=1,4
        npp = ngrp(ig,ie,jk1)
!
        if(npp.ne.0) then
          do 120 ig2=1,4
            npp2=ncm(npp,ig2)
!            write(6,*) 'clape ',npp,npp2
!
            if(npp2.ne.0) then
              igg=igg+1
              mgrp(igg,ie,jec) = npp2             
            endif
 120      continue
        endif
 110  continue
 100  continue
!
      return
 500  format(i1)
 501  format(a1)
 510  format(3x, 'C',i1,a1,' ->  C',i1,a1,$)
 511  format(3x,'IC',i1,a1,' -> IC',i1,a1,$)
      end
