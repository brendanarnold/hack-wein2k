      subroutine fermi5(etemp,e,efmod)
!---------------------------------------------------------------------
!
! Authors: Paola Alippi
!   and    Guido Vielsack
!          FHI(MPG)  Berlin
! Modified by L. D. Marks to include Mermin functional, and more accurate
! iteration
! 
!---------------------------------------------------------------------
!     input via vector-file
!     output: weigh(k,n), and ef
!
!             delef -> in common-block /mkef/ to control right
!                      output of fermi energy in lapw2.
!                      in lapw2 the source code has to look like follows:
!                        CALCULATE FERMI-ENERGY AND WEIGHT OF K-POINTS
!                        delef = 0.d0
!                        CALL FERMI
!                        WRITE(6,1060)  ELECN,EF + delef
!                        WRITE(21,1060) ELECN,EF + delef
!---------------------------------------------------------------------
        USE param; USE parallel
        USE bandm
        USE com,only : weigh,ef,elecn,xwt,jspin=>nspin,nat,nband,rel,nk,nb,minwav,maxwav
        USE kpp1
        USE xa2, only: weight,ne
      implicit real*8 (a-h, o-z)
!
      parameter ( eps = 1.d-10 )
      CHARACTER *10    KNAME,ACONT
      CHARACTER *67    ERRMSG                                            
      CHARACTER *80    FNAME,VECFN,VECFND,enefn,enefnd
      LOGICAL          EFERM
      REAL*8           e(2*nkpt,nume)
      Character*5 efmod
      common /celect/betha,k
      common /mkef/delef,ts2
!para begin
      COMMON /PROC/    VECFN,VECFND,enefn,enefnd
      COMMON /IPROC/   IPROC
!para end		
!
      external electr
!---------------------------------------------------------------------
      betha = 1.0d0 / etemp
!
      acont = 'CONT      ' 
      enorm = 2.0d0
      itap = 30
      jspin = 1
!
      minwav = 18000
      maxwav = 0        
!
      elecn = elecn - 1.0d-10
      index = 0
      sumw = 0.0d0
!  
      do 10 i=1,nume
         do 15 j=1,nkpt
            WEIGH(j,i) = 0.d0
15       continue
10    continue
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!.....read vector file:
!
!para begin
! ensure to mimick the proper vector file
! if running on multiple processors
      iloop=0
	iloop0=0
      k=0	
 777  continue
      if(IPROC.GT.0) then
         iloop=iloop+1
         iloop0=iloop0+1

         close(29)
         close(30)
         call mknam(FNAME,enefn,ILOOP)
         open(30,FILE=FNAME,STATUS='old',ERR=950)
         call mknam(FNAME,enefnd,ILOOP)
         open(29,FILE=FNAME,STATUS='unknown',ERR=950)
      endif
!para end

 1003 continue
      DO I=1,NAT 
         READ(30,'(f9.5)') EMIST
         READ(30,'(f9.5)') EMIST
      ENDDO
      DO I=1,NAT                                                  
         READ(29,'(f9.5)',END=1002) EMIST
         READ(29,'(f9.5)',END=1002) EMIST
      ENDDO
      JSPIN=2                                                          
 1002 CONTINUE
!
!para begin
!         k = 0
!para end
 4       k=k+1 
!para begin
         if(iloop0.ne.0) KPP(ILOOP0)=K
! testing
!      write(*,*)'reading k=',K
!para end
      if(k.gt.2*nkpt) goto 900
      READ(ITAP,'(3e19.12,a10,2i6,f5.1)',END=999) SS,T,ZZ,KNAME,N,NE(K),WEIGHT(K)               
      nmat=MAX(n,nmat)
            if(ne(k).gt.nume) GOTO 920
            if (kname.eq.acont) goto 1
            index=index+ne(k)
            sumw=sumw+weight(k)       
            if(n.gt.maxwav) maxwav=n   
            if(n.lt.minwav) minwav=n 
!            read (itap) (kx(i), ky(i), kz(i), i=1,n)
 14         read (itap,*) num, e1
!
               e(k,num) = e1    
               if(itap.eq.30.and.e1.gt.ebmax(num)) ebmax(num)=e1
               if(itap.eq.30.and.e1.lt.ebmin(num)) ebmin(num)=e1
!               read (itap) (a(i), i=1,n)
            if (num.lt.ne(k)) goto 14       
         goto 4
!
 1       k=k-1
         enorm=enorm+2.0d0
      goto 1003
!
!para begin
! we still have vector files to read

!  999 K=K-1
 999  CONTINUE
      K=K-1
      nk=k
!      write(6,*)'skipping last k point'
      IF(ILOOP.LT.IPROC) goto 777
	
   
!para end
!
      if(itap.eq.30.and.jspin.eq.2) then
         itap=29
         sumup=sumw
!para begin
         iloop=0
!cccc         K=0
         IF(IPROC.EQ.0) GOTO 4
         GOTO 777
!         goto 4
!para end

      else if(itap.eq.29.and.jspin.eq.2) then
         nk=k/2
         dif=2*sumup-sumw
         IF(ABS(DIF).GT.eps) GOTO 930
      end if
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  renormalization of weight(k,n) by SUM_k.n [weight(k,n)]  
!
      fact = enorm /sumw
!
      sumw = 0.d0
      do 30 kk=1,k
         weight(kk) = weight(kk) * fact
         sumw = sumw + weight(kk)
         nmax = ne(kk)
30    continue
!
      write(6,120)sumw
120   format('    sum over weight(k) (after  renom.) = ',f10.5)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!.....find emax and emin
!
      emax = e(1,1)
      emin = e(1,1)
      do 35 ik=1,k
         do 40 ie=1,ne(ik)
            emax = max(emax,e(ik,ie)) 
            emin = min(emin,e(ik,ie)) 
40       continue
35    continue
!
      write(6,'(9h  emin = ,f10.5,9h  emax = ,f10.5)')emin,emax
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!.....call minimizing function  and  thus calcutating chemical potential  
!     emu
!
      emu = rtbis(electr,emin,emax,eps,e)
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!.....ef is the upper energy value with weight > 0
!
      ef = emu + 30.d0*etemp
      if( ef .gt. emax ) goto 941
      delef = emu -ef
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!.....renomalize charge
!
      sum = 0.d0
      do 45 ik=1,k
         do 50 ie=1,ne(ik)
!
            de = e(ik,ie) - emu
            fint = dexp( de * betha )
            focc = 1.0d0 / (fint + 1.0d0)
!
            weigh(ik,ie) = weight(ik) * focc
            sum = sum + weigh(ik,ie)
!
50       continue
 45   continue
!
      fact = ELECN/sum
!
      do 55 ik=1,k
         do 60 ie=1,ne(ik)
            weigh(ik,ie) = weigh(ik,ie) * fact
 60      continue
55    continue
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!.....entropy calculation
!     ( for spinpol.systems divide by 2 
!
      entr = 0.0d0
!
      do 65 ik = 1, k
         do 70 ie = 1, ne(ik)
!
            de = e(ik,ie) - emu
            fint = dexp( de * betha )
            focc = 1.0d0 / (fint + 1.0d0)
!
            f1 = focc
            f2 = 1.0d0 - f1
            if (f1.gt.0.0d0.and.f1.lt.1.0d0) then
               eint = f1*dlog(f1) + f2*dlog(f2)
               entr = entr + weight(ik) * eint
            endif
70       continue
 65   continue
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!.....output:
!
      ts = etemp * entr / jspin
      if(efmod.eq.'TEMPS')then
!             Include full entropy term, Mermin Functional
              ts2= ts 
      else
!             Include only half as an approximate T=0 form
              ts2= ts / 2.d0
      endif
      if(myid.eq.0) then
      write (6,151) entr
      write (21,151) entr
      write(6,152) etemp
      write(6,153) ts
      if(efmod .eq. 'TEMPS')then
        write(21,157) ts
      else
        write(21,154) ts2
      endif
      write(6,156)emu
      write(21,156)emu
      endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 151   format (9x,' -S / Kb ',10X,'= ',F12.8)
 152   format (9X,' Kb * T  ',10X,'= ',F12.8)
 153   format (9X,' -T*Entr ',10X,'= ',F12.8)
 154   format (9X,' -(T*S)/2',10X,'= ',F12.8)
 156   format (9X,' Chem Pot',10X,'= ',F12.8)
 157   format (9X,' -(T*S)  ',10X,'= ',F12.8)

!
       return
!
!        Error messages
!
 900  WRITE (ERRMSG,9000) NKPT
      CALL OUTERR('FERMI5',ERRMSG)
      STOP  'FERMI5 - Error'
 920  WRITE (ERRMSG,9020) nume,ne(k) 
      CALL OUTERR('FERMI5',ERRMSG)
      STOP  'FERMI5 - Error'
 930  CALL OUTERR('FERMI5',' SUMUP AND SUMDN IN FERMI NOT EQUAL')
      CALL OUTERR('FERMI5',' CHECK UP AND DOWN INPUTS OF LAPW1 ')
      WRITE (ERRMSG,9030) SUMUP,SUMW,DIF             
      CALL OUTERR('FERMI5',ERRMSG)
      STOP  'FERMI5 - Error'
 941  CALL OUTERR('FERMI5',' emax too low in vector-file')
      STOP  'FERMI5 - Error'
!para begin
 950  WRITE (ERRMSG,9060)
      CALL OUTERR('LAPW2',ERRMSG)
      WRITE (ERRMSG,9070) FNAME
      CALL OUTERR('LAPW2',ERRMSG)
      STOP 'FERMI - Error'
!para end
 9000 FORMAT ('NUMBER OF K-POINTS .GT. NKPT =',i6)
 9020 FORMAT ('nume,ne(k)',2i6)
 9030 FORMAT ('SUMUP,SUMW,DIF', 3f10.7)
!para begin
 9060 FORMAT('can''t open unit: ',I2)
 9070 FORMAT('       filename: ',A50)
!para end
       end
!-----------------------------------------------------------------------
       double precision function electr(emu,e)
!
!......calculates the difference between the real number of electrons
!      and the actual (corresponding to emu)
!
        USE com,only : weigh,ef,elecn,xwt,jspin=>nspin,nat,nband,rel,nk,nb,minwav,maxwav
         USE xa2, only: weight,ne
         USE param
       implicit REAL*8 (a-h,o-z)
!
      REAL*8  e(2*nkpt,nume)
      common /celect/betha,k
!     Use higher precision to nail the fermi energy
!!!      real*16 de,buse,sum,focc,fint,dt,sum2,tmp

      buse=betha
      sum = 0.d0
      do 10 ik=1,k
         sum2=0.d0
         tmp=weight(ik)
         do 15 ie=1,ne(ik)
            de = e(ik,ie) - emu
            dt = de * buse
!.... fix for overflow
        if ( dt .lt. 200. ) then
            fint = exp( dt )
            focc = 1.0d0 / (fint + 1.0d0)
             sum2= sum2 + focc
        endif
15       continue
        sum = sum + sum2*tmp
 10   continue

!     Old version
!
!      sum = 0.d0
!      do 10 ik=1,k
!         do 15 ie=1,ne(ik)
!            de = e(ik,ie) - emu
!!.... fix for overflow
!            if( de * betha .gt. 80.d0 ) de = 80.d0 / betha
!!
!            fint = dexp( de * betha )
!            focc = 1.0d0 / (fint + 1.0d0)
!            sum = sum + focc*weight(ik)
!15       continue
! 10   continue
!!
      electr = sum - ELECN
!
      return
      end
!
!--------------------------------------------------------------------
      double precision FUNCTION RTBIS(FUNC,X1,X2,XACC,e)
      implicit double precision (a-h,o-z)
      PARAMETER (JMAX=50)
      FMID=FUNC(X2,e)
      F=FUNC(X1,e)
      IF(F*FMID.GE.0.D0) GOTO 900
      IF(F.LT.0.)THEN
        RTBIS=X1
        DX=X2-X1
      ELSE
        RTBIS=X2
        DX=X1-X2
      ENDIF
      DO 11 J=1,JMAX
        DX=DX*.5D0
        XMID=RTBIS+DX
        FMID=FUNC(XMID,e)
        IF(FMID.LE.0.D0)RTBIS=XMID
        IF(ABS(DX).LT.XACC .OR. FMID.EQ.0.D0) RETURN
 11   CONTINUE
      GOTO 910
 900  if(abs(FMID) .lt. 1D-9)then
        RTBIS=X2
        return
      else if(abs(F) .lt. 1D-9)then
        RTBIS=X1
        return
      endif
      CALL OUTERR('RTBIS','Root must be bracketed for bisection.')
      STOP  'FERMI5 - Error'
 910  CALL OUTERR('RTBIS','too many bisections in rtbis')
      STOP  'FERMI5 - Error'
      END 

