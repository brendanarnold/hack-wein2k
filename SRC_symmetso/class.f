      subroutine class(jatom,isym,s,lattic,spos,mult,isplit, &
                      name,rotold,iatnr)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      logical s(48)                    
      character*4 lattic
      character*79 name
      dimension lm(2,49),zvec(3),yvec(3),rotloc(3,3),rotold(3,3)
      dimension spos(3,*)
      dimension iatnr(*)
      label=0
      do 1 i=1,48
      lm(1,i)=0
 1    lm(2,i)=0
      lmmax=0
 900  format('lm:',121(2i2,1x))
!
      if(lattic(1:1).eq.'H'.or.lattic(1:1).eq.'R') then
         if(s(2).and.s(24).and.isym.eq.24) then
            write(6,*) ' pointgroup is 6/mmm (neg. iatnr!!)'
            write(6,*) ' axes should be: 6 || z, m n z, m n y'
            ixy=2
            call dirdeh(zvec,s,'  ',label)
 810        call dirdeh(yvec,s,'my',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            if(label.ne.0) goto 810
            isplit=4
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(2,0,6,0,1,lm,lmmax)
         elseif(s(2).and.s(24).and.isym.eq.12) then
            write(6,*) ' pointgroup is 6/m (neg. iatnr!!)'
            write(6,*) ' axes should be: 6 || z, m n z'
               ixy=0
            call dirdeh(zvec,s,'  ',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=4
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(2,0,6,0,1,lm,lmmax)
            call kurki(2,0,6,0,-2,lm,lmmax)
         elseif(s(24).and.isym.eq.12) then
            write(6,*) ' pointgroup is -6m2 (neg. iatnr!!)'
            write(6,*) ' axes should be: -6 || z, m n y'
            ixy=2
            call dirdeh(zvec,s,'  ',label)
 811        call dirdeh(yvec,s,'my',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            if(label.ne.0) goto 811
            isplit=4
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(2,0,6,0,1,lm,lmmax)
            call kurki(2,1,6,3,1,lm,lmmax)
         elseif(s(4).and.s(10).and.isym.eq.12) then
            write(6,*) ' pointgroup is 622 (neg. iatnr!!)'
            write(6,*) ' axes should be: 6 || z, m n y'
            ixy=2
            call dirdeh(zvec,s,'  ',label)
 812        call dirdeh(yvec,s,'my',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            if(label.ne.0) goto 812
            isplit=4
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(2,0,6,0,1,lm,lmmax)
            call kurki(2,1,6,0,-1,lm,lmmax)
         elseif(s(18).and.s(10).and.isym.eq.12) then
            write(6,*) ' pointgroup is 6mm (neg. iatnr!!)'
            write(6,*) ' axes should be: 6 || z, m n y'
            ixy=2
            call dirdeh(zvec,s,'  ',label)
 813        call dirdeh(yvec,s,'my',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            if(label.ne.0) goto 813
            isplit=4
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(1,0,6,0,1,lm,lmmax)
         elseif(s(2).and.isym.eq.12) then
            write(6,*) ' pointgroup is -3m (neg. iatnr!!)'
            write(6,*) ' axes should be: -3 || z, m n y'
            ixy=2
            call dirdeh(zvec,s,'  ',label)
 814        call dirdeh(yvec,s,'my',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            if(label.ne.0) goto 814
            isplit=4
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(2,0,3,0,1,lm,lmmax)
         elseif(s(2).and.isym.eq.6) then
            write(6,*) ' pointgroup is -3 (neg. iatnr!!)'
            write(6,*) ' axes should be: -3 || z'
            ixy=0
            call dirdeh(zvec,s,'  ',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=4
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(2,0,3,0,1,lm,lmmax)
            call kurki(2,0,3,0,-2,lm,lmmax)
         elseif(s(24).and.isym.eq.6) then
            write(6,*) ' pointgroup is -6 (neg. iatnr!!)'
            write(6,*) ' axes should be: -6 || z'
               ixy=0
            call dirdeh(zvec,s,'  ',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=4
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(2,0,6,0,1,lm,lmmax)
            call kurki(2,0,6,0,-2,lm,lmmax)
            call kurki(2,1,6,3,1,lm,lmmax)
            call kurki(2,1,6,3,-2,lm,lmmax)
         elseif(s(10).and.isym.eq.6) then
            write(6,*) ' pointgroup is 6 (neg. iatnr!!)'
            write(6,*) ' axes should be: 6 || z'
               ixy=0
            call dirdeh(zvec,s,'  ',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=4
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(1,0,6,0,1,lm,lmmax)
            call kurki(1,0,6,0,-2,lm,lmmax)
         elseif((s(4).or.s(5).or.s(6).or.s(7)).and.isym.eq.6) then
            write(6,*) ' pointgroup is 32 (neg. iatnr!!)'
            write(6,*) ' axes should be: 3 || z, 2 || y'
            ixy=2
            call dirdeh(zvec,s,'  ',label)
 815        call dirdeh(yvec,s,'2y',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            if(label.ne.0) goto 815
            isplit=4
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(2,0,3,0,1,lm,lmmax)
            call kurki(2,1,3,0,-1,lm,lmmax)
         elseif((s(18).or.s(19).or.s(15).or.s(16)).and.isym.eq.6) then
            write(6,*) ' pointgroup is 3m (neg. iatnr!!)'
            write(6,*) ' axes should be: 3 || z, m n y'
            ixy=2
            call dirdeh(zvec,s,'  ',label)
 816        call dirdeh(yvec,s,'my',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            if(label.ne.0) goto 816
            isplit=4
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(1,0,3,0,1,lm,lmmax)
         elseif(isym.eq.3) then
            write(6,*) ' pointgroup is 3 (neg. iatnr!!)'
            write(6,*) ' axes should be: 3 || z'
               ixy=0
            call dirdeh(zvec,s,'  ',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=4
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(1,0,3,0,1,lm,lmmax)
            call kurki(1,0,3,0,-2,lm,lmmax)
         elseif(isym.eq.8) then
            write(6,*) ' pointgroup is mmm (neg. iatnr!!)'
            write(6,*) ' axes should be: m n z, m n y, m n x'
            ixy=2
            call dirdeh(zvec,s,'mz',label)
 817        call dirdeh(yvec,s,'my',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            if(label.ne.0) goto 817
            isplit=8
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(2,0,2,0,1,lm,lmmax)
         elseif(s(2).and.isym.eq.4) then
            write(6,*) ' pointgroup is 2/m (neg. iatnr!!)'
            write(6,*) ' axes should be: 2 || z, m n z'
            ixy=0
            call dirdeh(zvec,s,'2z',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=8
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(2,0,2,0,1,lm,lmmax)
            call kurki(2,0,2,0,-2,lm,lmmax)
         elseif((s(14).or.s(15).or.s(16).or.s(17).or.s(18).or.s(19) &
                      .or.s(20)).and.isym.eq.4) then
            write(6,*) ' pointgroup is mm2 (neg. iatnr!!)'
            write(6,*) ' axes should be: 2 || z, m n y'
            ixy=2
            call dirdeh(zvec,s,'2z',label)
 818        call dirdeh(yvec,s,'my',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            if(label.ne.0) goto 818
            isplit=8
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(1,0,2,0,1,lm,lmmax)
         elseif(isym.eq.4) then
            write(6,*) ' pointgroup is 222 (neg. iatnr!!)'
            write(6,*) ' axes should be: 2 || z, 2 || y'
            ixy=2
            call dirdeh(zvec,s,'2z',label)
 819        call dirdeh(yvec,s,'2y',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            if(label.ne.0) goto 819
            isplit=8
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(2,0,2,0,1,lm,lmmax)
            call kurki(2,1,2,0,-1,lm,lmmax)
         elseif(s(2).and.isym.eq.2) then
            write(6,*) ' pointgroup is -1 (neg. iatnr!!)'
            write(6,*) ' axes should be: any'
            ixy=0
            call dirdeh(zvec,s,'  ',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=8
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(2,0,1,0,1,lm,lmmax)
            call kurki(2,0,1,0,-2,lm,lmmax)
         elseif((s(14).or.s(15).or.s(16).or.s(17).or.s(18).or.s(19) &
                      .or.s(20)).and.isym.eq.2) then
            write(6,*) ' pointgroup is m (neg. iatnr!!)'
            write(6,*) ' axes should be: m n z'
            ixy=0
            call dirdeh(zvec,s,'mz',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=8
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            lm(1,2)=1
            lm(2,2)=1
            lm(1,3)=-1
            lm(2,3)=1
            lm(1,4)=2
            lm(1,5)=2 
            lm(2,5)=2
            lm(1,6)=-2
            lm(2,6)=2
            lm(1,7)=3 
            lm(2,7)=1
            lm(1,8)=-3
            lm(2,8)=1
            lm(1,9)=3 
            lm(2,9)=3
            lm(1,10)=-3
            lm(2,10)=3
            lm(1,11)=4
            lm(1,12)=4
            lm(2,12)=2
            lm(1,13)=-4
            lm(2,13)=2
            lm(1,14)=4 
            lm(2,14)=4
            lm(1,15)=-4
            lm(2,15)=4
            lm(1,16)=5
            lm(2,16)=1
            lm(1,17)=-5
            lm(2,17)=1
            lm(1,18)=5
            lm(2,18)=3
            lm(1,19)=-5
            lm(2,19)=3
            lm(1,20)=5
            lm(2,20)=5
            lm(1,21)=-5
            lm(2,21)=5
            lm(1,22)=6
            lm(1,23)=6
            lm(2,23)=2
            lm(1,24)=-6
            lm(2,24)=2
            lm(1,25)=6
            lm(2,25)=4
            lm(1,26)=-6
            lm(2,26)=4
            lm(1,27)=6
            lm(2,27)=6
            lm(1,28)=-6
            lm(2,28)=6
            lmmax=28
         elseif(isym.eq.2) then
            write(6,*) ' pointgroup is 2 (neg. iatnr!!)'
            write(6,*) ' axes should be: 2 || z'
            ixy=0
            call dirdeh(zvec,s,'2z',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=8
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(1,0,2,0,1,lm,lmmax)
            call kurki(1,0,2,0,-2,lm,lmmax)
         elseif(isym.eq.1) then
            call kurki(1,0,1,0,1,lm,lmmax)
            call kurki(1,0,1,0,-2,lm,lmmax)
            write(6,*) ' pointgroup is 1'
            write(6,*) ' axes should be: any'
            ixy=0
            call dirdeh(zvec,s,'  ',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=8
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
         else
            write(6,*) ' no hex-pointgroup found' 
         endif
!.....end of hexagonal lattic
      else
         if(isym.eq.48) then
            write(6,*) ' pointgroup is m3m (pos. iatnr!!)'
            write(6,*) ' axes should be: any'
            ixy=0
            call dirdef(zvec,s,'  ',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=2
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'lt')
            lm(1,2)=4
            lm(1,3)=4
            lm(2,3)=4
            lm(1,4)=6
            lm(1,5)=6
            lm(2,5)=4
            lmmax=5
       elseif(s(2).and.isym.eq.24) then
            write(6,*) ' pointgroup is m3 (pos. iatnr!!)'
            write(6,*) ' axes should be: any'
            ixy=0
            call dirdef(zvec,s,'  ',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=2
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'lt')
            lm(1,2)=4
            lm(1,3)=4
            lm(2,3)=4
            lm(1,4)=6
            lm(1,5)=6
            lm(2,5)=4
            lm(1,6)=6
            lm(2,6)=2
            lm(1,7)=6
            lm(2,7)=6
            lmmax=7
!            write(6,*) 'warning: second L=6 harmonics not programmed'
         elseif(s(9).and.isym.eq.24) then
            write(6,*) ' pointgroup is 432 (pos. iatnr!!)'
            write(6,*) ' axes should be: any'
            ixy=0
            call dirdef(zvec,s,'  ',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=2
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'lt')
            lm(1,2)=4
            lm(1,3)=4
            lm(2,3)=4
            lm(1,4)=6
            lm(1,5)=6
            lm(2,5)=4
            lmmax=5
         elseif(isym.eq.24) then
            write(6,*) ' pointgroup is -43m (pos. iatnr!!)'
            write(6,*) ' axes should be: any'
            ixy=0
            call dirdef(zvec,s,'  ',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=2
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'lt')
            lm(1,2)=4
            lm(1,3)=4
            lm(2,3)=4
            lm(1,4)=6
            lm(1,5)=6
            lm(2,5)=4
            lm(1,6)=-3
            lm(2,6)=2
            lmmax=6
         elseif(s(2).and.isym.eq.12) then
            write(6,*) ' pointgroup is -3m (neg. iatnr!!)'
            write(6,*) ' axes should be: -3 || z, m n y'
            call kurki(2,0,3,0,1,lm,lmmax)
            ixy=2
            call dirdef(zvec,s,'-3',label)
 800        call dirdef(yvec,s,'my',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            if(label.ne.0) goto 800
            isplit=4
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
         elseif(isym.eq.12) then
            write(6,*) ' pointgroup is 23 (pos. iatnr!!)'
            write(6,*) ' axes should be: any'
            ixy=0
            call dirdef(zvec,s,'  ',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=2
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'lt')
            lm(1,2)=4
            lm(1,3)=4
            lm(2,3)=4
            lm(1,4)=6
            lm(1,5)=6
            lm(2,5)=4
            lm(1,6)=6
            lm(2,6)=2
            lm(1,7)=6
            lm(2,7)=6
            lm(1,8)=-3
            lm(2,8)=2
            lmmax=8
!            write(6,*) 'warning: second L=6 harmonics not programmed'
         elseif(s(2).and.isym.eq.6) then
            write(6,*) ' pointgroup is -3 (neg. iatnr!!)'
            write(6,*) ' axes should be: -3 || z'
            ixy=0
            call dirdef(zvec,s,'-3',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=4
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(2,0,3,0,1,lm,lmmax)
            call kurki(2,0,3,0,-2,lm,lmmax)
         elseif((s(6).or.s(12).or.s(13)).and.isym.eq.6) then
            write(6,*) ' pointgroup is 3m (neg. iatnr!!)'
            write(6,*) ' axes should be: 3 || z, m n y'
            ixy=2
            call dirdef(zvec,s,'3z',label)
 801        call dirdef(yvec,s,'my',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            if(label.ne.0) goto 801
            isplit=4
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(1,0,3,0,1,lm,lmmax)
         elseif(isym.eq.6) then
            write(6,*) ' pointgroup is 32 (neg. iatnr!!)'
            write(6,*) ' axes should be: 3 || z, 2 || y'
            ixy=2
            call dirdef(zvec,s,'3z',label)
 802        call dirdef(yvec,s,'2y',label)
            if(label.ne.0) goto 802
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=4
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(2,0,3,0,1,lm,lmmax)
            call kurki(2,1,3,0,-1,lm,lmmax)
         elseif(isym.eq.3) then
            write(6,*) ' pointgroup is 3 (neg. iatnr!!)'
            write(6,*) ' axes should be: 3 || z'
            ixy=0
            call dirdef(zvec,s,'3z',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=4
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(1,0,3,0,1,lm,lmmax)
            call kurki(1,0,3,0,-2,lm,lmmax)
         elseif(isym.eq.16) then
            write(6,*) ' pointgroup is 4/mmm (neg. iatnr!!)'
            write(6,*) ' axes should be: 4 || z, m n z, m n x'
            ixy=2
            call dirdef(zvec,s,'4z',label)
 803        call dirdef(yvec,s,'my',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            if(label.ne.0) goto 803
            isplit=-2
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(2,0,4,0,1,lm,lmmax)
         elseif(s(2).and.(s(9).or.s(10).or.s(11)).and.isym.eq.8) then
            write(6,*) ' pointgroup is 4/m (neg. iatnr!!)'
            write(6,*) ' axes should be: 4 || z, m n z'
            ixy=0
            call dirdef(zvec,s,'4z',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=-2
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(2,0,4,0,1,lm,lmmax)
            call kurki(2,0,4,0,-2,lm,lmmax)
         elseif(s(2).and.isym.eq.8) then
            write(6,*) ' pointgroup is mmm (neg. iatnr!!)'
            write(6,*) ' axes should be: m n z, m n y, m n x'
            ixy=2
            call dirdef(zvec,s,'mz',label)
 804        call dirdef(yvec,s,'my',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            if(label.ne.0) goto 804
            isplit=8
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(2,0,2,0,1,lm,lmmax)
         elseif((s(32).or.s(33).or.s(34)).and.isym.eq.8) then
            write(6,*) ' pointgroup is -42m (neg. iatnr!!)'
            write(6,*) ' axes should be: -4 || z, 2 || x'
            ixy=1
            call dirdef(zvec,s,'-4',label)
 805        call dirdef(yvec,s,'2x',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            if(label.ne.0) goto 805
            isplit=-2
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(2,0,4,0,1,lm,lmmax)
            call kurki(2,1,4,2,-1,lm,lmmax)
         elseif((s(6).or.s(7).or.s(8).or.s(12).or.s(13)) &
                    .and.isym.eq.8) then
            write(6,*) ' pointgroup is 4mm (neg. iatnr!!)'
            write(6,*) ' axes should be: 4 || z, m n y'
            ixy=2
            call dirdef(zvec,s,'4z',label)
 806        call dirdef(yvec,s,'my',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            if(label.ne.0) goto 806
            isplit=-2
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(1,0,4,0,1,lm,lmmax)
         elseif(isym.eq.8) then
            write(6,*) ' pointgroup is 422 (neg. iatnr!!)'
            write(6,*) ' axes should be: 4 || z, 2 || y'
            ixy=2
            call dirdef(zvec,s,'4z',label)
 807        call dirdef(yvec,s,'2y',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            if(label.ne.0) goto 807
            isplit=-2
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(2,0,4,0,1,lm,lmmax)
            call kurki(2,1,4,0,-1,lm,lmmax)
         elseif((s(9).or.s(10).or.s(11)).and.isym.eq.4) then
            write(6,*) ' pointgroup is 4 (neg. iatnr!!)'
            write(6,*) ' axes should be: 4 || z'
            ixy=0
            call dirdef(zvec,s,'4z',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=-2
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(1,0,4,0,1,lm,lmmax)
            call kurki(1,0,4,0,-2,lm,lmmax)
         elseif((s(32).or.s(33).or.s(34)).and.isym.eq.4) then
            write(6,*) ' pointgroup is -4 (neg. iatnr!!)'
            write(6,*) ' axes should be: -4 || z'
            ixy=0
            call dirdef(zvec,s,'-4',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=-2
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(2,0,4,0,1,lm,lmmax)
            call kurki(2,0,4,0,-2,lm,lmmax)
            call kurki(2,1,4,2,1,lm,lmmax)
            call kurki(2,1,4,2,-2,lm,lmmax)
         elseif(s(2).and.isym.eq.4) then
            write(6,*) ' pointgroup is 2/m (neg. iatnr!!)'
            write(6,*) ' axes should be: 2 || z, m n z'
            ixy=0
            call dirdef(zvec,s,'2z',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=8
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(2,0,2,0,1,lm,lmmax)
            call kurki(2,0,2,0,-2,lm,lmmax)
         elseif((s(6).or.s(7).or.s(8).or.s(12).or.s(13).or.s(14) &
           .or.s(15).or.s(16).or.s(17)).and.isym.eq.4) then
            write(6,*) ' pointgroup is mm2 (neg. iatnr!!)'
            write(6,*) ' axes should be: 2 || z, m n y'
            ixy=2
            call dirdef(zvec,s,'2z',label)
 808        call dirdef(yvec,s,'my',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            if(label.ne.0) goto 808
            isplit=8
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(1,0,2,0,1,lm,lmmax)
         elseif(isym.eq.4) then
            write(6,*) ' pointgroup is 222 (neg. iatnr!!)'
            write(6,*) ' axes should be: 2 || z, 2 || y'
            ixy=2
            call dirdef(zvec,s,'2z',label)
 809        call dirdef(yvec,s,'2y',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            if(label.ne.0) goto 809
            isplit=8
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(2,0,2,0,1,lm,lmmax)
            call kurki(2,1,2,0,-1,lm,lmmax)
         elseif(s(2).and.isym.eq.2) then
            write(6,*) ' pointgroup is -1 (neg. iatnr!!)'
            write(6,*) ' axes should be: any'
            ixy=0
            call dirdef(zvec,s,'  ',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=8
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(2,0,1,0,1,lm,lmmax)
            call kurki(2,0,1,0,-2,lm,lmmax)
         elseif((s(6).or.s(7).or.s(8).or.s(12).or.s(13).or.s(14) &
            .or.s(15).or.s(16).or.s(17)).and.isym.eq.2) then
            write(6,*) ' pointgroup is m (neg. iatnr!!)'
            write(6,*) ' axes should be: m n z'
            ixy=0
            call dirdef(zvec,s,'mz',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=8
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            lm(1,2)=1
            lm(2,2)=1
            lm(1,3)=-1
            lm(2,3)=1
            lm(1,4)=2
            lm(1,5)=2 
            lm(2,5)=2
            lm(1,6)=-2
            lm(2,6)=2
            lm(1,7)=3 
            lm(2,7)=1
            lm(1,8)=-3
            lm(2,8)=1
            lm(1,9)=3 
            lm(2,9)=3
            lm(1,10)=-3
            lm(2,10)=3
            lm(1,11)=4
            lm(1,12)=4
            lm(2,12)=2
            lm(1,13)=-4
            lm(2,13)=2
            lm(1,14)=4 
            lm(2,14)=4
            lm(1,15)=-4
            lm(2,15)=4
            lm(1,16)=5
            lm(2,16)=1
            lm(1,17)=-5
            lm(2,17)=1
            lm(1,18)=5
            lm(2,18)=3
            lm(1,19)=-5
            lm(2,19)=3
            lm(1,20)=5
            lm(2,20)=5
            lm(1,21)=-5
            lm(2,21)=5
            lm(1,22)=6
            lm(1,23)=6
            lm(2,23)=2
            lm(1,24)=-6
            lm(2,24)=2
            lm(1,25)=6 
            lm(2,25)=4
            lm(1,26)=-6
            lm(2,26)=4
            lm(1,27)=6 
            lm(2,27)=6
            lm(1,28)=-6
            lm(2,28)=6
            lmmax=28
         elseif(isym.eq.2) then
            write(6,*) ' pointgroup is 2 (neg. iatnr!!)'
            write(6,*) ' axes should be: 2 || z'
            ixy=0
            call dirdef(zvec,s,'2z',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=8
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
            call kurki(1,0,2,0,1,lm,lmmax)
            call kurki(1,0,2,0,-2,lm,lmmax)
         elseif(isym.eq.1) then
            call kurki(1,0,1,0,1,lm,lmmax)
            call kurki(1,0,1,0,-2,lm,lmmax)
            write(6,*) ' pointgroup is 1 (neg. iatnr!!)'
            write(6,*) ' axes should be: any'
            ixy=0
            call dirdef(zvec,s,'  ',label)
            call matrot(zvec,yvec,ixy,rotloc,rotold,label)
            isplit=8
      call strwri(spos,mult,isplit,name,rotloc,iatnr(jatom),'gt')
         else
            write(6,*) ' no pointgroup found, isym =',isym 
         endif
      endif
      if (iatnr(jatom).lt.0) call lmsort(lm,lmmax)
      write(6,900) (lm(1,i),lm(2,i),i=1,lmmax)
      write(54,903) (lm(1,i),lm(2,i),i=1,lmmax)
 903  format(121(i3,i2))
      write(6,*) '=============================================='
      return
      end
