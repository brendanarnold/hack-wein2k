      SUBROUTINE write1(norb,nqn,nk,den,zel,valel,mult,lcore,ecore)
      IMPLICIT REAL*8(A-H,O-Z)
      dimension nqn(30),nk(30),den(30),zel(30,2)
      character*28 file(8)
      logical l0,l1,l2,l3,lcore(30)
!
!      ecore=-6.d0
      esc1=-3.d0
      esc2=-1.d0
      offset=0.3
!
      l0=.true.
      l1=.true.
      l2=.true.
      l3=.true.
      icor=0
      ival=0
!
      do 10 i=1,norb
      etest=den(i)*2.d0
!.....core-state below -6.ry
      lcore(i)=.false.
      if(etest.lt.ecore) then
         icor=icor+1
         lcore(i)=.true.
!.....semi-core state below -3.Ry, low local orbital 
      else if(etest.lt.esc1) then
         valel1=valel
         valel=valel+(zel(i,1)+zel(i,2))*mult
         if(nk(i).eq.(-1)) then
            ival=ival+2
            l0=.false.
            write(file(ival-1),100) 0,    offset,0.000,'CONT'             
            write(file(ival),100) 0,etest+offset,0.005,'STOP'             
         else if(nk(i).eq.(-2)) then
            if((nk(i-1).eq.1.and.den(i-1)*2.d0.lt.ecore).or. &
               (nk(i+1).eq.1.and.den(i+1)*2.d0.lt.ecore)) then
               icor=icor+1
               lcore(i)=.true.
               valel=valel1
               goto 10
            endif
            ival=ival+2
            l1=.false.
            write(file(ival-1),100) 1,    offset,0.000,'CONT'           
            write(file(ival),100) 1,etest+offset,0.005,'STOP'           
         else if(nk(i).eq.(-3)) then
            if((nk(i-1).eq.2.and.den(i-1)*2.d0.lt.ecore).or. &
               (nk(i+1).eq.2.and.den(i+1)*2.d0.lt.ecore)) then
               icor=icor+1
               lcore(i)=.true.
               valel=valel1
               goto 10
            endif
            ival=ival+2
            l2=.false.
            write(file(ival-1),100) 2,    offset,0.000,'CONT'              
            write(file(ival),100) 2,etest+offset,0.005,'STOP'            
         else if(nk(i).eq.(-4)) then
            if((nk(i-1).eq.3.and.den(i-1)*2.d0.lt.ecore).or. &
               (nk(i+1).eq.3.and.den(i+1)*2.d0.lt.ecore)) then
               icor=icor+1
               lcore(i)=.true.
               valel=valel1
               goto 10
            endif
            ival=ival+2
            l3=.false.
!            write(file(ival-1),100) 3,   offset,0.000,'CONT'            
            write(file(ival-1),100) 3,   offset,0.000,'CONT'            
            write(file(ival),100) 3,etest+offset,0.005,'STOP'
         end if            
!.....semi-core state below -1.Ry, high local orbital
      else if(etest.lt.esc2) then
         valel=valel+(zel(i,1)+zel(i,2))*mult
         if(nk(i).eq.(-1).and.l0) then
            ival=ival+2
            l0=.false.
            write(file(ival-1),100) 0,etest+offset,0.010,'CONT'           
            write(file(ival),100) 0,        offset,0.000,'CONT'          
         else if(nk(i).eq.(-2).and.l1) then
            ival=ival+2
            l1=.false.
            write(file(ival-1),100) 1,etest+offset,0.010,'CONT'         
            write(file(ival),100) 1,        offset,0.000,'CONT'         
         else if(nk(i).eq.(-3).and.l2) then
            ival=ival+2
            l2=.false.
            write(file(ival-1),100) 2,etest+offset,0.010,'CONT'         
            write(file(ival),100) 2,        offset,0.000,'CONT'         
         else if(nk(i).eq.(-4).and.l3) then
            ival=ival+2
            l3=.false.
            write(file(ival),100) 3,etest+offset,0.010,'CONT'         
            write(file(ival-1),100) 3,        offset,0.000,'CONT'         
!            write(file(ival-1),100) 3,        offset,0.000,'CONT'         
         end if            
!.....valence state
      else 
         valel=valel+(zel(i,1)+zel(i,2))*mult
         if(nk(i).eq.(-1).and.l0) then
            ival=ival+1
            write(file(ival),100) 0,     offset,0.000,'CONT'        
         else if((nk(i).eq.(-2).or.nk(i).eq.1).and.l1) then
            ival=ival+1
            l1=.false.
            write(file(ival),100) 1,     offset,0.000,'CONT'       
         else if((nk(i).eq.2.or.nk(i).eq.(-3)).and.l2) then
            l2=.false.
            ival=ival+1
            write(file(ival),100) 2,     offset,0.010,'CONT'        
         else if((nk(i).eq.3.or.nk(i).eq.(-4)).and.l3) then
            l3=.false.
            ival=ival+1
            write(file(ival),100) 3,     offset,0.010,'CONT'        
         end if            
      end if
 10   continue
!
!.....write in1 file
      write(15,102) offset,ival,0
      do 15 i=1,ival
      write(15,101) file(i)
 15   continue
!
!.....write inc file
      if(icor.eq.0) then
          write(16,103) 1,0.d0
          write(16,104) nqn(1),nk(1),0
      else 
      write(16,103) icor,0.d0
      do 20 i=1,norb
      if(lcore(i)) write(16,104) nqn(i),nk(i),abs(nk(i))*2
 20   continue
      end if
!
 100  format(i2,f8.2,f11.3,1x,a4,' 1')
 101  format(a28)
 102  format(f6.2,i5,i3,6x,'(GLOBAL E-PARAMETER WITH n OTHER CHOICES, ' &
           ,'global APW/LAPW)')
 103  format(i2,f5.2,5x,'NUMBER OF ORBITALS (EXCLUDING SPIN), SHIFT')
 104  format(i1,',',i2,',',i1,15x,'( N,KAPPA,OCCUP)')
      RETURN
      END

