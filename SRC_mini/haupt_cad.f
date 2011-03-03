 
      subroutine haupt(errfn)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'param.inc'
!
!.....input-data (struct-file,tmpM-file) should be given with 
!     'high' numerical accuracy (at least 10 digits) to enshure
!     the internaly assumed accuracy of about 8 digits.
!
!
!.....tolx: tollerance of xyz-coordinates
      parameter (tolx = 0.5d-3 )
!

      dimension      br3(3,3),vecrot(3)
      dimension      weight(nato)
      dimension      f(3,nato),alat(3)

      dimension      xwert(mxpm,maxit), ewert(maxit),fwert(mxpm,maxit)
      dimension      q(nato,2,maxit)
      dimension      xnos(maxit),temp2(maxit)
      dimension      damp(mxpm)
      
      integer qparam
!      character*1     abo/"'"/
      CHARACTER*4     LATTIC,IREL,RELA,cform,minmod
      CHARACTER*10    ANAME
      CHARACTER*80    title1,errfn
      LOGICAL         ortho,ende,imsl,random
      COMPLEX*16     ROKNEW,ROKOLD,ROKOLD1,ROKOLD2,rokmix

      COMMON /ORTH/   ORTHO
      COMMON /STRUK/  POS(3,NDIF),AA,BB,CC,ALPHA(3),RMT(NATO),V(NATO), &
                      PIA(3),VOL,IATNR(NATO),MULT(NATO),ISPLIT(NATO)
      COMMON /CHAR/   LATTIC,ANAME(NATO)
      COMMON /GENER/  BR1(3,3),BR2(3,3)
      COMMON /ROTMAT/ ROTLOC(3,3,NATO),ROTIJ(3,3,NDIF),TAUIJ(3,NDIF)
      COMMON /SYMETR/ TAU(3,NSYM),KVCNEW(3,NWAV),IMAT(3,3,NSYM), &
                      INUM(NSYM),IORD
      common /fmin/   alat,xwert,ewert,fwert,nparam,iter,ifnc,igrdc    
      common /ctoler/ tolxx
      COMMON /XA/ CLMOLD(NRAD,NCOM,NATO,2),CLMOLD1(NRAD,NCOM,NATO,2), &
                CLMOLD2(NRAD,NCOM,NATO,2),CLMNEW(NRAD,NCOM,NATO,2), &
                ROKNEW(NWAV,2),ROKOLD(NWAV,2), &
                ROKOLD1(NWAV,2),ROKOLD2(NWAV,2),rokmix(NWAV,2), &
                R(NRAD),RNOT(NATO), &
                DX(NATO),Zatom(NATO),JRI(NATO), &
                LM(2,NCOM,NATO),LMMAX(NATO)

      save /xa/
      DATA             RELA/'RELA'/                                     
      DATA             RANDOM/.true./                                     

      tolxx = tolx
!------------------------------------------------------------------
      
      ifnc=0
      igrdc=0
      imsl=.false.
      dpi=4.d0 * datan(1.d0)
!.....DEFINE CONSTANTS                                                  
      DO 100 I=1,NATO                                                   
        MULT(I)=0
        V(I)=0.D0
  100 RMT(I)=0.D0
!       
!.....READ STRUCT
      READ(20,1000) TITLE1                                              
      READ(20,1010) LATTIC,NAT,cform,IREL
      if(nat.gt.nato) goto 911

!     READ IN LATTICE CONSTANTS                                         
      READ(20,1020) AA,BB,CC,ALPHA(1),ALPHA(2),ALPHA(3)
      alat(1)=aa
      alat(2)=bb
      alat(3)=cc
      IF(ABS(ALPHA(1)).LT.1.D-5) ALPHA(1)=90.0D0
      IF(ABS(ALPHA(2)).LT.1.D-5) ALPHA(2)=90.0D0
      IF(ABS(ALPHA(3)).LT.1.D-5) ALPHA(3)=90.0D0
!     INDeX COUNTS ALL ATOMS IN THE UNIT CELL, JATOM COUNTS ONLY THE    
!     NOTEQUIVALENT ATOMS                                               
      INDEX=0                                                           
      DO 50 JATOM = 1,NAT                                               
        INDEX=INDEX+1                                                  
        READ(20,1030) IATNR(JATOM),( POS(J,INDEX),J=1,3 ),MULT(JATOM),  &
          ISPLIT(jatom)
        IF (MULT(JATOM).EQ.0) THEN                                     
          WRITE(6,1040) JATOM,INDEX,MULT(JATOM)                       
          STOP ' MULT EQ 0'                                    
        ENDIF                                                          
        DO 55 M=1,MULT(JATOM)-1                                     
          INDEX=INDEX+1                                            
          READ(20,1031) IATNR(JATOM),( POS(J,INDEX),J=1,3)         
   55   CONTINUE                                                    
          READ(20,1050) ANAME(JATOM),JRI(JATOM),RNOT(JATOM),RMT(JATOM)      &
          ,zatom(jatom)     
         DX(JATOM)=LOG(RMT(jatom)/RNOT(JATOM)) / (JRI(JATOM)-1)
          READ(20,1051) ((ROTLOC(I,J,JATOM),I=1,3),J=1,3)                
 1051     FORMAT(20X,3F10.7)                                             
         V(JATOM)=4.D0*dPI/3.D0*RMT(JATOM)**3
   50 CONTINUE                                                          
      WRITE(6,800)                                                      
      WRITE(6,805)  TITLE1                                              
      WRITE(6,810)  LATTIC                                              
      WRITE(6,820)  AA,BB,CC                                            
      WRITE(6,840)  NAT                                                 
     
      
      call latgen(nat)
      
   43 FORMAT(3X,A77)                                                    
   44 FORMAT(I3,A77)                                                    
  700 FORMAT(I3,A77)                                                    
  800 FORMAT(////,30X,50(1H-),/,33X,'S T R U C T U R A L   ',            &
        'I N F O R M A T I O N',/,30X,50(1H-),//)                  
  805 FORMAT(3X,'SUBSTANCE',20X,'= ',A80,/)                             
  810 FORMAT(3X,'LATTICE',22X,'= ',A4)                                  
  820 FORMAT(3X,'LATTICE CONSTANTS ARE',8X,'= ',3F12.7)                 
  830 FORMAT(3X,'SYMMETRY ATTITUDE IS',9X,'= ',A4)                      
  840 FORMAT(3X,'NUMBER OF ATOMS IN UNITCELL  = ',I3)                   
 1000 FORMAT(A80)                                                       
 1010 FORMAT(A4,24X,I2,1X,A4,/,13X,A4,14X,A4)                                 
 1020 FORMAT(6F10.6,10X,F10.7) 
!kew 1030 FORMAT(5X,I3,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)                    
 1030 FORMAT('ATOM=',I3,': X=',F10.7,1X,'Y=',F10.7,1X,'Z=',F10.7,/, &
       10X,'MULT=',I2,10X,'ISPLIT=',I2)
!kew 1031 FORMAT(5X,I3,4X,F10.7,3X,F10.7,3X,F10.7)          
 1031 FORMAT(5X,I3,': X=',F10.7,1X,'Y=',F10.7,1X,'Z=',F10.7)
 1040 FORMAT(1X,'Q(U)  :',12(2X,F7.4) )
 1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)
 1055 FORMAT('LOCAL ROT MATRIX:',2X,3F10.7)
 2020 FORMAT(1X,I4,' < NMAT < ',I4,3X,'SPIN=',I1,3X,'NATO=',I2) 
      
!     test existence of *.tmpM file
!     if not -> initialize and start
      
      iter=0
      read(16,*,end=1200) nparam
      read(16,*,end=1200) iter
      if(iter.gt.maxit) goto 912
      do 200 i=1,iter+1
        read(16,*) ewert(i),j,temp2(i)
        read(16,*) (xwert(j,i),fwert(j,i),xnos(i), &
        j=1,nparam)
  200 CONTINUE
         rewind(16)
         read(16,*) j
         read(16,*) j
         do 205 i=1,iter
           read(16,*) ewert(i),j,temp2(i)
           read(16,*) (xwert(j,i),fwert(j,i),xnos(i), &
                       j=1,nparam)
  205    CONTINUE
 1200 iter=iter+1
      
!     read input file (modus, damp, speed-parameters)
!.....tolf: tollerance for forces
!      
      read(5,1220) minmod,tolf
 1220 format (a4,f5.2)
!      if (minmod.eq.'NOSE') then
!         read(16,*,end=1210) dummy,i,temp2(iter)
!         read(16,*,end=1210) dummy
!         read(16,*,end=1210) dummy,xnos(iter)
!         rewind(16)
!         read(16,*,end=1210) j
!         read(16,*,end=1210) j
!         do 205 i=1,iter-1
!           read(16,*) ewert(i),i,temp2(i)
!           read(16,*) (xwert(j,i),fwert(j,i),xnos(i),
!     &                 j=1,nparam)
!  205    CONTINUE
!      endif
 1210 write(6,*) 'Minimization Method: ', minmod
      do 210 ia=1,nat
        read(5,*,end=910) dummy1,dummy2,dummy3,dummy4
        if (minmod.eq.'NEWT') then
           fwert(3*(ia-1)+1,iter)=dummy1
           fwert(3*(ia-1)+2,iter)=dummy2
           fwert(3*(ia-1)+3,iter)=dummy3
           damp(3*(ia-1)+1)=dummy4
           damp(3*(ia-1)+2)=dummy4
           damp(3*(ia-1)+3)=dummy4
        endif
        if (minmod.eq.'BFGS') then
           if(dummy1.ne.0.d0) fwert(3*(ia-1)+1,iter)=dble(mult(ia)) &
           *dummy1
           if(dummy2.ne.0.d0) fwert(3*(ia-1)+2,iter)=dble(mult(ia)) &
           *dummy2
           if(dummy3.ne.0.d0) fwert(3*(ia-1)+3,iter)=dble(mult(ia)) &
           *dummy3
        endif
        if (minmod.eq.'MOLD') then
           fwert(3*(ia-1)+1,iter)=-1.0d0
           fwert(3*(ia-1)+2,iter)=-1.0d0
           fwert(3*(ia-1)+3,iter)=-1.0d0
           weight(ia)=dummy1*1836.15d0
           delta=dummy2
        endif
        if (minmod.eq.'NOSE') then
           fwert(3*(ia-1)+1,iter)=-1.0d0
           fwert(3*(ia-1)+2,iter)=-1.0d0
           fwert(3*(ia-1)+3,iter)=-1.0d0
           weight(ia)=dummy1*1836.15d0
           delta=dummy2
           temp=dummy3
           fnose=dummy4
        endif
  210 CONTINUE
!     read new results for energy and forces 
      ende=.true.
      qparam=0
 1400 format(13x,4d15.7)
      read(15,*) ewert(iter)
      fwert1=0.d0
      index=1
      do 230 ia=1,nat
         read(15,*) (f(ik,ia),ik=1,3)
!kew  17/11/00 start
!kew---> since unit 15 contains forces in units of the internal coordinate
!kew     system, the rotation is not necessary anymore!
!kew         write (6,*) 'force before rotate:',(f(ik,ia),ik=1,3)
!kew         CALL inv(ROTLOC(1,1,ia),br3,det)
!kew         CALL ROTATE (f(1,ia),br3,VECROT)
!ccc         CALL ROTATE (f(1,ia),rotloc(1,1,ia),VECROT)
!kew         f(1,ia)=vecrot(1)
!kew         f(2,ia)=vecrot(2)
!kew         f(3,ia)=vecrot(3)
!kew         write (6,*) 'force after rotate:',(f(ik,ia),ik=1,3)
!kew        if (.not.ortho) then
!            tpi = 8.d0 * datan(1.d0)
!            do 240 i1=1,3
!               do 250 i2 = 1,3
!                  br1(i1,i2) = br1(i1,i2) / tpi*alat(i1)
!  250          CONTINUE
!  240       CONTINUE
!            CALL inv(br1,br3,det)
!            CALL ROTATE (f(1,ia),br3,vecROT)
!kew         vecrot(1)=f(1,ia)*BR1(1,1)+f(2,ia)*BR1(2,1)+f(3,ia)*BR1(3,1)
!kew         vecrot(2)=f(1,ia)*BR1(1,2)+f(2,ia)*BR1(2,2)+f(3,ia)*BR1(3,2)
!kew         vecrot(3)=f(1,ia)*BR1(1,3)+f(2,ia)*BR1(2,3)+f(3,ia)*BR1(3,3)
!kew            f(1,ia)=vecrot(1)
!kew            f(2,ia)=vecrot(2)
!kew            f(3,ia)=vecrot(3)
!kew            write (6,*) 'force after ortho:',(f(ik,ia),ik=1,3)
!ad         IF(LATTIC(1:1).EQ.'R') then
!ad            alat(1)=DSQRT(1.0d0/3.0d0*aa*aa+1.0d0/9.0d0*cc*cc)
!ad            alat(2)=alat(1)
!ad            alat(3)=alat(1)
!ad         endif
!kew         endif
!kew  17/11/00 end
         do 260 ik=1,3
            qparam=qparam+1
         if(abs(f(ik,ia)).lt.1.e-6) f(ik,ia)=0.d0
            if(fwert(qparam,iter).gt.0.000001d0) then
              if (abs(f(ik,ia)).gt.tolf) ende=.false.
              write(6,*) 'force converg.:natom,i',ia,ik,f(ik,ia),tolf
            endif
            f(ik,ia)=f(ik,ia)/1000.0d0
            fff=-f(ik,ia)
            if (minmod.eq.'MOLD'.or.minmod.eq.'NOSE') ende=.false.
            if (iter.eq.1) then
               xwert(qparam,iter)=pos(ik,index) * alat(ik)
            endif
            fwert1=fwert1+fwert(qparam,iter)
            fwert(qparam,iter)=fff * fwert(qparam,iter)
  260    CONTINUE
         index=index+mult(ia)
  230 CONTINUE
      nparam=qparam
            tolf1=tolf/1000.d0*fwert1/qparam
      write(6,*) 'mean tolf value:',tolf1
!     re-write tape 16
      rewind 16 

      write(6,*) 'number of parameters', nparam
      ewert(iter)=dmod(ewert(iter),1.d2)
      write(6,*) 'energy',ewert(iter)
      write(17,'(2i3)') nparam
      write(17,'(i5)') iter
      do 270 i=1,iter
         write(17,'(f20.5,i10,33x,f12.7)') ewert(i),i,temp2(i)
         write(17,'(3d20.12)')  &
              (xwert(j,i),fwert(j,i), &
               xnos(i),j=1,nparam)
  270 CONTINUE

      if (ende) then
        endfile 21
!
!.......output of final result:
        call finish( &
        xwert(1,iter),fwert(1,iter),ewert(iter),mult,nato,alat,nparam, &
        ' S T O P  (FORCES SMALL)', &
        'best result in case.struct-file and case.finM-file', &
        'stop_forces_small')
!
      CALL ERRCLR(ERRFN)
      CALL ERRFLG(ERRFN,'STOP in MINI, FORCES small')
        stop '>>  (mini) forces = 0 -> exit'
      endif
!
!

      if(minmod.eq.'BFGS') then
        call dfpmin(xwert(1,1),nparam,iiter,fret,nat)
        goto 900
      else if (minmod.eq.'NEWT') then
        call nwtmin(xwert(1,1),nparam,iiter,fret,nat,damp(1),tolf1)
        goto 900
      else if (minmod.eq.'MOLD') then
        if (iter.eq.1) then
        call mold0(xwert(1,iter),xwert(1,iter+1), &
                  nparam,fwert(1,iter),delta,weight,temp)
        temp2(iter+1)=temp
        else
        call mold(xwert(1,iter-1),xwert(1,iter),xwert(1,iter+1), &
                  nparam,fwert(1,iter),delta,weight,temp)
        temp2(iter+1)=temp
        endif
        goto 1203
      else if (minmod.eq.'NOSE') then
        if (iter.eq.1) then
           if (random) then
           call nose(xwert(1,10),xwert(1,iter),xwert(1,iter+1), &
                     20000.0d0,0.0d0,xnos(iter+1), &
                  nparam,fwert(1,iter),fnose,temp,delta,weight)
           temp2(iter+1)=temp
           else 
           call nose0(xwert(1,iter),xwert(1,iter+1), &
                  xnos(iter),xnos(iter+1), &
                  nparam,fwert(1,iter),fnose,temp,delta,weight)
           temp2(iter+1)=temp
           endif
        else
        call nose(xwert(1,iter-1),xwert(1,iter),xwert(1,iter+1), &
                  xnos(iter-1),xnos(iter),xnos(iter+1), &
                  nparam,fwert(1,iter),fnose,temp,delta,weight)
        temp2(iter+1)=temp
        endif
        goto 1203
      else 
         stop
      endif

!      write(6,7893) iter,ifnc,igrdc
!      write(6,*) ' fret: ',fret
!
!.......output of final result:
!        call finish(
!     &  xwert(1,iter),fwert(1,iter),ewert(iter),mult,nato,alat,nparam,
!     &  ' S T O P  (minimization ok)',
!     &  'best result in case.struct-file and case.finM-file',
!     &  'stop_min_ok')
!
!      stop '>>  (mini) from minimization -> exit'
!      
!      entry outp
        
 900  do 280 j=1,nparam
        if (xwert(j,iter).gt.0d0) goto 1203
 280  CONTINUE
      endfile 21
      stop '>>> (mini) xwert = 0 -> exit'

 1203 continue
      
!     berechnen der neuen atompositionen
      
!     conversation position von iter+1 xwert in positionen
!     write new positions to tape16 too
      qparam=0
      index=0
      call latmix(nat)
      write(6,*) jri
      if (minmod.eq.'NOSE'.or.minmod.eq.'MOLD') then
         write(17,'(f20.7,i10,33x,f12.7)') 0.,iter+1,temp2(iter+1)
      else
         write(17,'(f20.7,i10,33x,f12.7)') 0.,iter+1,0.
      endif
      do 290 ia=1,nat
        index=index+1
        index1=index
        do 300 ir=1,3
          qparam=qparam+1
          pos(ir,index)= xwert(qparam,iter+1) / alat(ir) 
          if (minmod.eq.'NOSE') then
             write(17,'(3d20.12)')  &
             pos(ir,index)*alat(ir),0.,xnos(iter+1)
          else
             write(17,'(3d20.12)')  &
             pos(ir,index)*alat(ir),0.,0.
          endif
          pos(ir,index)=pos(ir,index)+1.d0
 305      IF (pos(ir,index).lt.0d0) THEN 
             pos(ir,index)=pos(ir,index)+1.d0
	     goto 305 
          endif
          pos(ir,index) =  &
          dmod(  pos(ir,index) + 1.d0 , 1.d0 )
 300    CONTINUE
        do 310 im=1,mult(ia)-1
          index=index+1
!     rotieren und translatieren (aber umgekehrt)
          do 320 ir =1,3
            xxx=0.0d0
            do 330 ir1=1,3
              xxx=xxx+rotij(ir1,ir,Index)*POS(ir1,INDEX1)
 330        CONTINUE
!cc old            xxx=xxx-tauij(ir,index)
            xxx=xxx+tauij(ir,index)
            if (xxx.lt.0d0) xxx=xxx+1.d0
            xxx=dmod( xxx + 1.d0, 1.d0 )
            pos(ir,index)=xxx
 320      CONTINUE
 310    CONTINUE
 290  CONTINUE
!    write a new scf-file
      jspin=1
      call wrtscf(nat,iter,jspin,minmod,temp2,q)

      call mixin(q,xwert,iter,nparam,nat,jspin)
      
!     new struct-file on tape 21

      WRITE(21,1000) TITLE1
      WRITE(21,1010) LATTIC,NAT,cform,IREL
      fac=180.d0/dpi
!     WRITE IN LATTICE CONSTANTS
      WRITE(21,1020) AA,BB,CC,ALPHA(1)*fac,ALPHA(2)*fac,ALPHA(3)*fac
!     INDEX COUNTS ALL ATOMS IN THE UNIT CELL, JATOM COUNTS ONLY THE    
!     NOTEQUIVALENT ATOMS                                               
      INDEX=0                                                           
      DO 530 JATOM = 1,NAT
        INDEX=INDEX+1                                                  
        WRITE(21,1030) IATNR(JATOM),(POS(J,INDEX),J=1,3 ),MULT(JATOM) &
          ,ISPLIT(JATOM)
        DO 535 M=1,MULT(JATOM)-1                                     
          INDEX=INDEX+1                                            
          WRITE(21,1031)IATNR(JATOM),(POS(J,INDEX),J=1,3)
  535   CONTINUE                                                    
        WRITE(21,1050) ANAME(JATOM),JRI(JATOM),RNOT(JATOM),RMT(JATOM), &
          zatom(jatom)     
        WRITE(21,1055) ((ROTLOC(I,J,JATOM),I=1,3),J=1,3)
  530 CONTINUE                                                          
      
      write(21,10) IORD                                                  
      DO 5 J=1,IORD                                                     
    5 write(21,11) ( (IMAT(J1,J2,J),J1=1,3),TAU(J2,J), J2=1,3),J          
!.....READ IN SYMMETRY OPERATIONS AND NONPRIMITIVE TRANSLATIONS    
   10 FORMAT(I4,6X,'NUMBER OF SYMMETRY OPERATIONS')                                            
   11 FORMAT(3(3I2,F10.8/),i4)                                              

!    a couple of informations on tape 6

      ende = .true.
      write(6,*)
      write(6,'(a9,a4,a4,2a10,2a10)')  &
          '(mini)','ia','ik','F','x(old)','x(new)','delta x'
      qparam=0
      do 350 ia=1,nat
!        write(6,7890) ia
         do 360 ir=1,3 
            qparam=qparam+1
!           write(6,7891) ir, f(ir,ia),fwert(qparam,iter)
            xdiff=xwert(qparam,iter+1)-xwert(qparam,iter)
            write(6,'(9x,i4,i4,2f10.5,2f10.5)')  &
                  ia,ir,f(ir,ia),xwert(qparam,iter), &
                  xwert(qparam,iter+1),xdiff
            if (dabs(xdiff).gt.tolx) ende=.false.
 360     CONTINUE
 350  CONTINUE
      write(6,*)
      write(6,7893) iter,ifnc,igrdc
 7890 format(1x,'Atom: ',i3)
 7891 format(1x,'Component: ',i2/1x, &
        'Force: ',f10.5/1x, &
        'Gradient: ',f10.5)
 7892 format(1x,'Move: ',f10.5,' --> ' &
        ,f10.5,'  Delta: ',f10.5)
 7893 format(1x,'iter, fcn-calls, grad-calls: ',3i5)

      if(ende) endfile 21
      return
!
!        error handling
!
  910 INFO = 1
!
!        error in case.inM
!
      CALL OUTERR('MINI','error in reading case.inM')
      CALL OUTERR('MINI','check if 4 parameters per atom given')
      GOTO 999
 912  CALL OUTERR('MINI','MAXIT too small')
      GOTO 999
!
 911  CALL OUTERR('MINI','NATO too small')
  999 STOP 'MINI - Error'
      end
