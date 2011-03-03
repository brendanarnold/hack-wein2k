!$hp9000_800 intrinsics on
      program spag
!
!.....Generates postscript file showing the band structure of the 
!     Kohn-Sham eigenvalues obtained from the Wien2k program.
!     Includes averaging over degenerate bands (SO) and multiple atoms
!
!     last updated:  Dec. 2001   Peter Blaha, Vienna Techn University
!                                f77->f90
!                    Mars 2002   Clas Persson, Uppsala University
!                                incorporated irreducible representations
!                                for drawing lines of the bands.
!                                Peter Blaha, testing
!
!                    May 2006:   Juergen Spitaler, University of Leoben 
!                                (juergen.spitaler@mu-leoben.at)
!                                incorporated the generation of the xmgrace
!                                file case.bands.agr, which can directly be
!                                opened with xmgrace.  
!                                In xmgrace, the size of all symbols can
!                                be changed in the menu
!                                Plot -> Graph appearance -> Special -> Z normalization

      use reallocate
      use irr_param
      IMPLICIT REAL*8 (A-H,O-Z)
      character  aline*80,fname*80
      CHARACTER*11      STATUS,FORM
      character  k_pat*7,ei_pat*28
      character  dummy*1,label1*12
      character  symbol*12,label*12
      character  title*80,lattice*4
      character  backslash*5
      logical    flg1
!
      integer,allocatable :: jatom_list(:)         ! NATO)
      real*8,pointer ::  eigen(:,:),charr(:,:)     ! NEVL,NKP)
      integer,pointer ::  n_ene(:),lines(:)        ! NKP)
      real*8,pointer ::  vk(:,:)                   ! 3,NKP)
      real*8,pointer ::  xval(:)                   ! NKP)
      logical,pointer ::    break(:)               ! nkp)
      character*12,pointer ::  k_name(:)           ! nkp
      dimension  qtl(13),icomma(15)
!...xmgrace: define variables
      CHARACTER*80   insp, case_name
      character*10 dat_dummy, string_dummy
      character  xmlabel1*50, xmtitle*12
      integer    ii1, jj1, index_shift
!...xmgrace end      
      data  k_pat     /'     K='/
      data  ei_pat /'EIGENVALUES BELOW THE ENERGY'/

!-----------------------------------------------------------------------
!
!.....INITIALIZE VARS
      n_kpt=0
      NEVL=1000
      NKP=1000
      allocate (  eigen(nevl,nkp),charr(nevl,nkp) )
      allocate (  n_ene(nkp),k_name(nkp),lines(nkp)    )
      allocate (  vk(3,NKP))
      allocate (  xval(NKP))
      allocate (  break(nkp))
!
      allocate (  eneik(NEVL,NKP) )
      allocate (  ngrp(4,NEVL,NKP), ngde(4,NEVL,NKP) ) !(4,NEVL,NKP) 
      allocate (  ndeg(NEVL,NKP) )
      allocate (  neik(nkp),ikg(NKP),ntab4(NKP) )
      allocate (  iscol(NEVL),iscol2(NEVL),iste(NEVL),ibnd(NEVL) )  !NEVL)
      allocate (  lkg(48,48,NKP) )
      allocate (  nlkg(48,NKP)   )
      allocate (  grpnam(NKP)    )
      allocate (  cnam(48,NKP)   )
!
      iarg=iargc()
      if(iarg.ne.1) goto 900
      call getarg(1,fname)
      OPEN(1,FILE=fname,STATUS='OLD',ERR=8000)
 8003 READ(1,*,END=8001) IUNIT,FNAME,STATUS,FORM,IRECL
      OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM, &
      ERR=8002)
!...xmgrace: get file name of the .insp file
      IF(IUNIT.EQ.5) INSP=FNAME
!...xmgrace end      
      GOTO 8003
 8000 WRITE(*,*) ' ERROR IN OPENING SPAG.DEF !!!!'
      STOP 'SPAG.DEF'
 8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
      WRITE(*,*) '       FILENAME: ',FNAME,'  STATUS: ',STATUS, &
      '  FORM:',FORM
      STOP 'OPEN FAILED'
 8001 CONTINUE
!.....write date 
      call wrtdate(6)
!
      label='            '
      label1='            '
      do i=1,12
      if(fname(i:i).ne.'.') then
         label(i:i)=fname(i:i)
      else
!         label(i:i)=','
         goto 100
      endif
!...xmgrace: define title of the plot
      xmtitle(1:12)=label(1:12)
!...xmgrace end            
      enddo
!
!.....READ K-VECTORS
!
 100  CONTINUE
      read(7,'(a80)',end=200) aline
      if  (aline(1:7).eq.k_pat)  then
!         THIS IS THE BEGINNING OF AN EIGENVALUE SECTION (output1)
          n_kpt=n_kpt + 1
          if(n_kpt.gt.nkp) then
          nkp=nkp*2
      call doreallocate(eigen,nevl,nkp)
      call doreallocate(charr,nevl,nkp)
      call doreallocate(n_ene,nkp)
      call doreallocate(k_name,nkp)
      call doreallocate(lines,nkp)
      call doreallocate(vk,3,nkp)
      call doreallocate(xval,nkp)
      call doreallocate(break,nkp)
      call doreallocate (  eneik,NEVL,NKP )
      call doreallocate (  ngrp,4,NEVL,NKP)
      call doreallocate (  ngde,4,NEVL,NKP)  
      call doreallocate (  ndeg,NEVL,NKP) 
      call doreallocate (  neik,nkp)
      call doreallocate (  ikg,NKP) 
      call doreallocate (  ntab4,NKP)
      call doreallocate (  lkg,48,48,NKP)
      call doreallocate (  nlkg,48,NKP)
      call doreallocate (  grpnam,NKP)
!
             write(6,*) ' parameter NKP increased:',nkp
          end if
          n_ene(n_kpt)=0
          call get_k (aline,vk(1,n_kpt),vk(2,n_kpt), &
                      vk(3,n_kpt),k_name(n_kpt))
          read(7,'(a1)') dummy
          read(7,'(a1)') dummy
!         READ EIGENVALUES OF K-VECTOR
 110      read(7,'(a80)') aline
          if (aline(15:42).eq.ei_pat)  then
            goto 100
          else
             if(n_ene(n_kpt).ge.nevl-4) then
      nevl=nevl*2
      call doreallocate(eigen,nevl,nkp)
      call doreallocate(charr,nevl,nkp)
      call doreallocate (  eneik,NEVL,NKP )
      call doreallocate (  ngrp,4,NEVL,NKP)
      call doreallocate (  ngde,4,NEVL,NKP)  
      call doreallocate (  ndeg,NEVL,NKP) 
      call doreallocate (  iscol,NEVL)
      call doreallocate (  iscol2,NEVL)
      call doreallocate (  iste,NEVL)
      call doreallocate (  ibnd,NEVL)
                   write(6,*) 'parameter nevl increased:',nevl
             endif
             call get_ei(aline,eigen(1,n_kpt), &
                             n_ene(n_kpt))
             goto 110
          endif
      else
          goto 100
      endif
!
!.....ALL K-VECTORS HAVE BEEN READ; SEARCH FOR K-POINT WITH SMALLEST
!     NUMBER OF EIGENVALUES
!
 200  continue
      write(*,*) 'number of k-points read in case.vector=',n_kpt
      nu_min=999
      do 205 j=1,n_kpt
         if (n_ene(j).lt.nu_min)  then
            nu_min=n_ene(j)
            k_min=j
         endif
 205  continue
      write(6,*) 'smallest number eigenvalues at k=',k_min,' (', &
       k_name(k_min),')'
      write(6,*) '         =',nu_min
!
!.....READ IRREDUCIBLE REPRESENTATIONS OF THE EIGENVECTOR
      nikir=0
      read(30,fmt=511,end=444,err=444) tjunk   
!
      read(30,511) tjunk
      write(6,*) 'using irred. representations'
 440  read(30,fmt=502,end=445) ikir, skir(1), skir(2), skir(3)
        nikir=nikir+1
        read(30,504) grpnam(nikir),ikg(nikir),ntab4(nikir) 
        if(grpnam(nikir).eq.'Nop') then
          write(*,*) 'k-point nr ',nikir,' not treated with irrep'
          write(6,*) 'k-point nr ',nikir,' cannot be characterized according to any point group.'
          write(6,*) 'here, we draw bands between closest energy values.'
          write(6,*) 'thus, we treat the k-point as point group C1.'
          grpnam(nikir)='C1 '
        endif 
        read(30,*)
        do 460 itab=1,ntab4(nikir)
          read(30,507) cnam(itab,nikir)
!
!.........find number of group elements
          flg1=.true.
          nlkg(itab,nikir)=1
          icnam=1
          do while(flg1) 
            select case ( cnam(itab,nikir)(icnam:icnam) )
               case ('E')
                 flg1=.false.
               case ('I')
                 flg1=.false.
               case ('C')
                 flg1=.false.
               case ('2')
                 flg1=.false.
                 nlkg(itab,nikir)=2
               case ('3')
                 flg1=.false.
                 nlkg(itab,nikir)=3
               case ('4')
                 flg1=.false.
                 nlkg(itab,nikir)=4
               case ('6')
                 flg1=.false.
                 nlkg(itab,nikir)=6
               case ('8')
                 flg1=.false.
                 nlkg(itab,nikir)=8
            end select
            icnam=icnam+1
            if(icnam.eq.6) then 
              
              write(*,*) 'ERROR spag: cannot find nr: ',cnam(itab,nikir)
              write(6,*) 'ERROR spag: cannot find nr: ',itab,nikir,icnam
              write(6,*) 'ERROR spag: cannot find nr: ',cnam(itab,nikir)
!              call flush(6)

              stop 'spag: cannot find nr of pntgrp '
            endif  
          enddo
          read(30,508)  (lkg(itab,inge,nikir), inge=1,nlkg(itab,nikir))
 460    continue
!
        read(30,503) neik(nikir)
        write(6,510) 'ik=',nikir,'; k=',skir(1), skir(2), skir(3), &
                     'point group= ',   grpnam(nikir), &
                     '; no of symm.ops=', ikg(nikir),    &
                     '; no of classes=',  ntab4(nikir),  &
                     'no of energies=', neik(nikir)
        do itab=1,ntab4(nikir)
          write(6,506) cnam(itab,nikir)
          write(6,508) (lkg(itab,inge,nikir),inge=1,nlkg(itab,nikir))
        enddo
        do 480 ie=1,neik(nikir)
          read(30,518) inum, ndeg(ie,nikir),eneik(ie,nikir), &
              (ngrp(ipgr,ie,nikir),ngde(ipgr,ie,nikir),ipgr=1,4)
 480    continue
!      call flush(6)
      goto 440
 444  write(6,*) 'do not use irred. representations'
 445  write(6,*) 'number of k-points read in case.irrep= ',nikir
!
 499  continue
 502  format(1i6,3f10.6)
 503  format(1i10)
 504  format(3x,a3,i3,i10)
 506  format(1x,a6,$)
 507  format(1x,a6)
 508  format(48i3)
 518  format(2i10,f10.6,4(2x,2i3))
 510  format(a3,i6,a4,3f10.6,/,a13,a3,a17,i5,a16,i5,/,a15,i5)
 511  format(a3)

!
!.....READ STRUCT FILE
!
      read(20,'(a80)') title
      read(20,'(a4,23x,i3)') lattice,natom
      write(6,*) 'lattice is=',lattice,natom,' atoms'
      allocate (jatom_list(natom+1))
      read(20,*)
      read(20,457) aa,bb,cc , alpha,beta,gamma 
        if(alpha.eq.0.d0) alpha=90.d0
        if(beta.eq.0.d0) beta=90.d0
        if(gamma.eq.0.d0) gamma=90.d0
 457  format(6f10.7)
      write(6,*) '         a=',aa
      write(6,*) '         b=',bb
      write(6,*) '         c=',cc
!
!.....TRANSFORM K-VECTORS TO CARTESIAN COORDINATES
!
      call cartco (vk,n_kpt,lattice,aa,bb,cc,alpha,beta,gamma)
      write(6,*) 'transforming k-points to cartesian coordinates...'
!
!.....SEARCH FOR THE ENDPOINTS OF THE STRAIGHT LINES THROUGH THE
!     BRILLOUIN-ZONE; CALCULATE THE CUMULATIVE DISTANCES BETWEEN
!     THE K-POINTS, WHICH ARE THE X-VALUES FOR THE BANDSTRUCTURE
!     PLOT
!
      write(6,*) 'checking for lines...'
      call bz_lin(vk,n_kpt,lines,n_lin,xval,nbreak,break)
      write(6,505) n_lin,nbreak
 505  format('numbers of BRILLOUIN-zone lines found:',i3,/, &
                  'numbers of breaks',i3)
!...xmgrace: write page properties to case.bands.agr
      write(40,*) '@ page size 595, 842'      
      write(40,*) '@ view 0.120000, 0.150000, 0.900000, 1.280000'
      write(40,*) '@ default linewidth 2.0'
      write(40,*) '@ xaxis  label char size 1.5'
      write(40,*) '@ xaxis  ticklabel char size 1.25'
      write(40,*) '@ yaxis  label char size 1.5'
      write(40,*) '@ yaxis  ticklabel char size 1.25'
      write(40,*) '@ xaxis  tick major grid on'
      write(40,*) '@ xaxis  tick spec type both'
      write(40,*) '@ xaxis  tick spec', n_lin+1
!...xmgrace end      
      do 210 j=1,n_lin+1
         if(j.le.n_lin) then
          write(6,*) 'line',j,' is: ',k_name(lines(j)),' to ', &
                     k_name(lines(j+1))
          write(6,*) 'xmax=',xval(lines(j))
         endif 
!...xmgrace: write xaxis tick labels to case.bands.agr;
!            translate Greek letters into xmgrace syntax
         backslash='" xG"'
         backslash(2:2)=achar(92)
         write(40,224) j-1, xval(lines(j))
224     format('@ xaxis  tick major ',i3,',',f8.5)
         if((k_name(lines(j))(1:5).eq.'GAMMA').or.(k_name(lines(j))(1:5).eq.'Gamma')) then
            backslash(4:4)='G'
            write(40,*) '@ xaxis  ticklabel ',j-1,',',backslash
!            write(40,*) '@ xaxis  ticklabel ',j-1,',','"\xG"'
            goto 210
         endif
         if((k_name(lines(j))(1:5).eq.'DELTA').or.(k_name(lines(j))(1:5).eq.'Delta')) then
            backslash(4:4)='D'
            write(40,*) '@ xaxis  ticklabel ',j-1,',',backslash
            goto 210
         endif
         if((k_name(lines(j))(1:6).eq.'LAMBDA').or.(k_name(lines(j))(1:6).eq.'Lambda')) then
            backslash(4:4)='L'
            write(40,*) '@ xaxis  ticklabel ',j-1,',',backslash
            goto 210
         endif
         if((k_name(lines(j))(1:5).eq.'SIGMA').or.(k_name(lines(j))(1:5).eq.'Sigma')) then
            backslash(4:4)='S'
            write(40,*) '@ xaxis  ticklabel ',j-1,',',backslash
            goto 210
         endif
         write(40,*) '@ xaxis  ticklabel ',j-1,',"',k_name(lines(j)),'"'
!...xmgrace end      
 210  continue
!
!.....READ VIEWING PARAMETERS
!
      call inview (ymin,ymax,xcm,ycm,xoffs,yoffs,height,efermi, &
           lyflag,nb_min,nb_max, &
           jatom,jatom_list,jtype,sizec,ticks,nticks) 
!...xmgrace: write xmin, ymin, xmax, ymax for the plot to case.bands.agr
      write(40,'(A)') '@ with g0'
      write(40,701) ymin,xval(lines(n_lin+1)),ymax
701     format('@ world 0, ',f12.5,',',f12.5,',',f12.5)
      write(40,*) "@ autoticks"
    !...xmgrace: label yaxis (written to case.bands.agr)
      if(lyflag.eq.1) then
         write(40,*) '@ yaxis  label "Energy (Ry)"'            
      else
         write(40,*) '@ yaxis  label "Energy(eV)"'
      endif
    ! ...xmgrace: insert line at Fermi energy (written to case.bands.agr)
      write(40,*) '@ with line'
      write(40,*) '@ line on'
      write(40,*) '@ line loctype world'
      write(40,222)  xval(lines(n_lin+1))
222     format('@ line 0, 0, ', f8.5,', 0')      
      write(40,*) '@ line linestyle 3'
      write(40,*) '@ line def'    
    ! ...xmgrace: label Fermi energy#
      write(40,*) '@ with string'
      write(40,*) '@ string on'
      write(40,*) '@ string loctype world'
      write(40,223) xval(lines(n_lin+1))+0.02
223     format('@ string ', f8.5,', -0.1')      
      write(40,*) '@ string char size 1.500000'
      write(40,*) '@ string def "E'//achar(92)//'sF"'
    !
!...xmgrace end
      xmin=0.d0
      xmax=xval(n_kpt)  
      if(lyflag.eq.2) call con_ev(n_ene,n_kpt, &
        efermi,eigen,nevl)
!
! ....READ QTL FILE
!
      char0=0.01d0
      do 208 i=1,n_kpt
      do 208 ii=1,n_ene(i)
 208  charr(ii,i)=char0
      if(jatom.eq.0) goto 206
 207  read(9,'(a80)',end=206) aline
      if(aline(2:5).eq.'JATO') then
         read(aline(8:9),'(i2)') iatom
         if(iatom.eq.jatom_list(1)) then
             i0=2
             icomma(1)=31
             do i=32,80
             if(aline(i:i).eq.',') then
               icomma(i0)=i
               i0=i0+1
             endif
             enddo
             icomma(i0)=icomma(i0-1)+2
             label1=aline(icomma(jtype)+1:icomma(jtype+1)-1)
             if(label1(1:1).eq.'0') label1(1:1)='s,'
             if(label1(1:1).eq.'1') label1(1:1)='p,'
             if(label1(1:1).eq.'2') label1(1:1)='d,'
             if(label1(1:1).eq.'3') label1(1:1)='f,'
!...xmgrace: dummy variable used 
             xmlabel1(:)=label1(:)
             index_shift=0
             do jj1=1,3
               do ii1=1,39
                 if(xmlabel1(ii1:ii1+1).eq.'PX') then
                    index_shift = index_shift+4
                    xmlabel1(ii1+6:12+index_shift) = xmlabel1(ii1+2:12+index_shift-4) 
                    xmlabel1(ii1:ii1+5) ='p'//achar(92)//'sx'//achar(92)//'N'
                 endif                    
                 if(xmlabel1(ii1:ii1+1).eq.'PY') then
                    index_shift = index_shift+4
                    xmlabel1(ii1+6:12+index_shift) = xmlabel1(ii1+2:12+index_shift-4) 
                    xmlabel1(ii1:ii1+5) ='p'//achar(92)//'sy'//achar(92)//'N'
                 endif                    
                 if(xmlabel1(ii1:ii1+1).eq.'PZ') then
                    index_shift = index_shift+4
                    xmlabel1(ii1+6:12+index_shift) = xmlabel1(ii1+2:12+index_shift-4) 
                    xmlabel1(ii1:ii1+5) ='p'//achar(92)//'sz'//achar(92)//'N'
                 endif                    
               enddo !ii1
               do ii1=1,39
                 if(xmlabel1(ii1:ii1+2).eq.'DZ2') then
                    index_shift = index_shift+6
                    xmlabel1(ii1+9:12+index_shift) = xmlabel1(ii1+3:12+index_shift-6) 
                    xmlabel1(ii1:ii1+8) ='d'//achar(92)//'sz'//achar(92)//'S2'//achar(92)//'N'
                 endif                    
                 if(xmlabel1(ii1:ii1+2).eq.'DXY') then
                    index_shift = index_shift+4
                    xmlabel1(ii1+7:12+index_shift) = xmlabel1(ii1+3:12+index_shift-4) 
                    xmlabel1(ii1:ii1+6) ='d'//achar(92)//'sxy'//achar(92)//'N'//achar(92)//'N'
                 endif                    
                 if(xmlabel1(ii1:ii1+2).eq.'DXZ') then
                    index_shift = index_shift+4
                    xmlabel1(ii1+7:12+index_shift) = xmlabel1(ii1+3:12+index_shift-4) 
                    xmlabel1(ii1:ii1+6) ='d'//achar(92)//'sxz'//achar(92)//'N'//achar(92)//'N'
                 endif                    
                 if(xmlabel1(ii1:ii1+2).eq.'DYZ') then
                    index_shift = index_shift+4
                    xmlabel1(ii1+7:12+index_shift) = xmlabel1(ii1+3:12+index_shift-4) 
                    xmlabel1(ii1:ii1+6) ='d'//achar(92)//'syz'//achar(92)//'N'//achar(92)//'N'
                 endif                    
               enddo !ii1
               do ii1=1,39
if(xmlabel1(ii1:ii1+4).eq.'DX2Y2') then
                    index_shift = index_shift+13
                    xmlabel1(ii1+18:12+index_shift) = xmlabel1(ii1+5:12+index_shift-13) 
      xmlabel1(ii1:ii1+17)='d'//achar(92)//'sx'//achar(92)//'S2'//achar(92)//'N'//achar(92)//'s-y'//achar(92)//'S2'//achar(92)//'N'
                 endif                    
               enddo !ii1
             enddo !jj1
         endif
!...xmgrace end      
      endif
      if  (aline(2:5).ne.'BAND')  goto 207  
      jband=0
 211  jband=jband+1     
      if(jband.ge.nb_min.and.jband.le.nb_max) then
         chmin=1.0d0
         chmax=0.0d0
         do 209 i=1,n_kpt
         do 209 ii=1,natom+1
           do jatom1=1,jatom
             if(ii.ne.jatom_list(jatom1)) goto 5210
             read(9,'(13X,F8.5,3X,12F8.5)',err=212,end=212) &
                 (qtl(j),j=1,13)
             charr(jband,i)=charr(jband,i)+qtl(jtype)*sizec
!            write(*,*) jband,i,charr(jband,i)
             if(chmin.gt.charr(jband,i)) chmin=charr(jband,i)
             if(chmax.lt.charr(jband,i)) chmax=charr(jband,i)
             goto 209
 5210      continue
           enddo
           read(9,*,err=212,end=212)            
 209     continue                            
         write(6,*) ' band',jband,'  char-min/max:',chmin,chmax   
      else
         do 213 i=1,n_kpt
         do 213 ii=1,natom+1
 213     read(9,*,err=212,end=212)
      end if
      if(jband.eq.nu_min) goto 206
      read(9,*,err=212,end=212)
      goto 211
 212  write(*,*) '        ERROR reading QTLs:'
      write(*,*) ' band:',jband,' k-point:',i
      write(*,*) ' execution continued ....'
 206  continue                         
!.....averaging character of degenerate bands (max 6 fold)
      do  jk=1,n_kpt
      je=1
 555  if(eigen(je,jk)-eigen(je+1,jk).gt.-5.d-6) then
        if(eigen(je,jk)-eigen(je+2,jk).gt.-5.d-6) then
          if(eigen(je,jk)-eigen(je+3,jk).gt.-5.d-6) then
            if(eigen(je,jk)-eigen(je+4,jk).gt.-5.d-6) then
              if(eigen(je,jk)-eigen(je+5,jk).gt.-5.d-6) then
              charr(je,jk)=(charr(je,jk)+charr(je+1,jk)+charr(je+2,jk)+ &
                     charr(je+3,jk)+charr(je+4,jk)+charr(je+5,jk))/6.d0
              charr(je+1,jk)=charr(je,jk)
              charr(je+2,jk)=charr(je,jk)
              charr(je+3,jk)=charr(je,jk)
              charr(je+4,jk)=charr(je,jk)
              charr(je+5,jk)=charr(je,jk)
              je=je+6
              goto 556
            endif
            charr(je,jk)=(charr(je,jk)+charr(je+1,jk)+charr(je+2,jk)+ &
                       charr(je+3,jk)+charr(je+4,jk))/5.d0
            charr(je+1,jk)=charr(je,jk)
            charr(je+2,jk)=charr(je,jk)
            charr(je+3,jk)=charr(je,jk)
            charr(je+4,jk)=charr(je,jk)
            je=je+5
            goto 556
          endif
          charr(je,jk)=(charr(je,jk)+charr(je+1,jk)+charr(je+2,jk)+ &
                       charr(je+3,jk))/4.d0
          charr(je+1,jk)=charr(je,jk)
          charr(je+2,jk)=charr(je,jk)
          charr(je+3,jk)=charr(je,jk)
          je=je+4
          goto 556
          endif
        charr(je,jk)=(charr(je,jk)+charr(je+1,jk)+charr(je+2,jk) &
                       )/3.d0
        charr(je+1,jk)=charr(je,jk)
        charr(je+2,jk)=charr(je,jk)
        je=je+3
        goto 556
        endif
      charr(je,jk)=(charr(je,jk)+charr(je+1,jk) &
                       )/2.d0
      charr(je+1,jk)=charr(je,jk)
      je=je+2
      goto 556
      endif
      je=je+1
 556      if(je.lt.n_ene(jk)) goto 555
      write(6,557) jk,(je,eigen(je,jk),charr(je,jk),je=1,n_ene(jk))
      enddo
 557  format(i4,(i5,2f10.6))
!
!.....initialize ps 
!
      pi=acos(-1.d0)
      write(6,*)xoffs,yoffs
      call psinit (xoffs,yoffs,xcm,ycm,dlwid)
      write(6,*) ' init done'
!     labelling:
         xval1=0.0d0
         call movet (xval1,ycm+1.d0)
         call writs (label,1.d0*height,'label ')
         if(jatom.gt.0) then
         write(label,878) (jatom_list(j),j=1,min(3,jatom))
 878     format(' atom ',3i2)
         call writs (label,1.d0*height,'label ')
         call writs (label1,1.d0*height,'label ')
!...xmgrace: write title to case.bands.agr
         if(jatom_list(1).gt.0) then
            write(40,123) xmtitle, jatom_list(1), xmlabel1
 123     format(' @ title "',a, ' atom ',i2,' ',a,'"')
         else
            write(40,124) xmtitle
 124        format(' @ title "',a,'"')
         endif
!...xmgrace end

         write(label,879) sizec
 879     format('  size',f5.2)
         call writs (label,1.d0*height,'label ')
         endif
!
!.....draw vertical lines
!
      write(6,*) n_lin
      write(11,*) '%vertical lines'
      do 220 j=1,n_lin+1
         xval1=xval(lines(j))*xcm/xmax
         write(11,140) dawid
         call movet (xval1,0.d0)
         call drawt (xval1,ycm)
 220  continue
      write(6,*) ' vertical lines done'
!
!.....draw horizontal lines and FERMI-energy
!
      write(11,*) '%horizontal lines and FERMI-energy'
      do 230 j=1,n_lin
      if (.not.break(lines(j+1))) then
	 xval1=xval(lines(j))*xcm/xmax
         xval2=xval(lines(j+1))*xcm/xmax
         call movet (xval1,0.d0)
         call drawt (xval2,0.d0)
         call movet (xval1,ycm)
         call drawt (xval2,ycm)
         write(11,*) 'stroke'
         write(11,140) dawid
!
         if (iprtf.ne.0)  then
!!         if (efermi.ne.999.d0)  then
            yval1=ycm/(ymax-ymin)*(efermi-ymin)
            write(11,'(A41)') 
            write(11,*) 'stroke'            
            write(11,140) dfwid
! 
            if(iprtf.eq.1) then
            call movet (xval1,yval1)
            call drawt (xval2,yval1)
            elseif(iprtf.eq.2) then
              xvaldel = (xval2-xval1)/18.5
              xval2   = xval1+xvaldel/2.d0
              do j1=1,19 
                call movet (xval1,yval1)
                call drawt (xval2,yval1)
                xval1 = xval1+xvaldel
                xval2 = xval1+xvaldel/2.d0
              enddo
            elseif(iprtf.eq.3) then
              xvaldel = (xval2-xval1)/41
              xval2   = xvaldel
              col='setgray'
              call movet (xval1,yval1)
              do j1=1,41
!c                call movet (xval1,yval1)
	        call transt(xval1,yval1)     
                call pointi (0,0.01d0,iprto) 
                call transt(-xval1,-yval1)
                xval1 = xval1+xvaldel                   
              enddo
            else
               stop 'error in input, incorrect line for Fermi'
            endif
!
            write(11,*) 'stroke'
            write(11,140) dawid
            write(6,*) ' EF-line done'
         endif
      endif
 230  continue
      write(6,*) ' horizontal lines and FERMI-energy done'
      write(11,*) 'stroke'
      write(11,140) dlwid
!
!.....label x-axis with names of k-points
!     write the names of the HIGH-SYMMETRY k-points
!
      do 240 j=1,n_kpt
      symbol=k_name(j) 
      if (symbol(1:1).eq.' ')  goto 240
      xval1=xval(j)*xcm/xmax
      call movet(xval1-0.2d0,-0.5d0)
      call writs (symbol,0.8d0*height,'labels')
 240  continue
      write(6,*) ' x-axis labeling done'
!
!.....label FERMI-energy
!
!      if (efermi.ne.999.d0)  then
      if (iprtf.ne.0)  then
         yval1=ycm/(ymax-ymin)*(efermi-ymin)
         call movet (xcm+0.1d0,yval1-0.2d0)
         call writs ('E           ',1.d0*height,'fermi ')
         call movet (xcm+0.5d0,yval1-0.4d0)
         call writs ('F           ',0.8d0*height,'fermi ')
         write(6,*) ' EF-labeling done'
      endif
!
!.....label y-axis with 'Energy [Ry] '
!
      yval1=ycm/2.d0-2.d0
      call movet (-2.7d0,yval1)
      write(11,*) '90. rotate'
      if (lyflag.eq.1)  then
         call writs ('Energy (Ry) ',1.d0*height,'label ')
      else if (lyflag.eq.2)  then
         call writs ('Energy (eV) ',1.d0*height,'label ')
      endif
      write(11,*) '-90. rotate'
      write(11,*) 'stroke'
      write(6,*) ' y-axis labeling done'
!
!.....label y-axis with ticks
!
      if (lyflag.eq.1)  then
         ymin1=ymin
 246     yval1=ycm/(ymax-ymin)*(ymin1-ymin)
         call movet (0.0d0,yval1)
         call drawt (-0.3d0,yval1)
         call writz (ymin1,0.7*height,yval1,lyflag)
         do 245 i=1,nticks
            if ((ymin1+ticks*i/(nticks+1)).ge.ymax) goto 247
            yval1=ycm/(ymax-ymin)*(ymin1+ticks*i/(nticks+1)-ymin)
            call movet (0.0d0,yval1)
            call drawt (-0.2d0,yval1)
 245     continue
         ymin1=ymin1+ticks
         if (ymin1.le.ymax) goto 246
 247     continue
      else if (lyflag.eq.2) then
         ymin1=efermi
 256     yval1=ycm/(ymax-ymin)*(ymin1-ymin)
         call movet (0.0d0,yval1)
         call drawt (-0.3d0,yval1)
         call movet (0.0d0,yval1)
         call writz (ymin1,0.7d0*height,yval1,lyflag)
         do 255 i=1,nticks
            if ((ymin1+ticks*i/(nticks+1)).ge.ymax) goto 257
            yval1=ycm/(ymax-ymin)*(ymin1+ticks*i/(nticks+1)-ymin)
            call movet (0.0d0,yval1)
            call drawt (-0.2d0,yval1)
 255     continue
         ymin1=ymin1+ticks
         if (ymin1.le.ymax) goto 256
 257     continue
         ymin1=efermi
 259     ymin1=ymin1-ticks
         if (ymin1.lt.ymin) goto 258
         yval1=ycm/(ymax-ymin)*(ymin1-ymin)
         call movet (0.0d0,yval1)
         call drawt (-0.3d0,yval1)
         call movet (0.0d0,yval1)
         call writz (ymin1,0.7d0*height,yval1,lyflag)
         do 265 i=1,nticks
            if ((ymin1+ticks*i/(nticks+1)).le.ymin) goto 258
            yval1=ycm/(ymax-ymin)*(ymin1+ticks*i/(nticks+1)-ymin)
            call movet (0.0d0,yval1)
            call drawt (-0.2d0,yval1)
 265     continue
         goto 259
 258     continue
      endif
!
!.....set clipboard
!
      write(11,*) 'stroke'
      if (nbreak.eq.0) then
         call clipin (0.d0,0.d0,xcm,ycm)
      else
      do 500 j=1,n_lin
         if(break(lines(j+1))) then
	    xval1=xval(lines(j))*xcm/xmax
            call clipin (0.d0,0.d0,xval1,ycm)
	    goto 501
         endif
 500  continue
 501  continue
      endif
!
!.....draw energies
!
      write(6,*) 'drawing energy-bands...'
      do 607 je2=1,NEVL
        iscol(je2) =je2
        iscol2(je2)=je2
 607  continue
!
!.....for each k-point
      do 250 jk=1,n_kpt
      if (break(jk)) then
         do 600 jjj=jk+1,n_kpt
	 if (break(jjj)) then
	    xval2=xval(jjj-1)*xcm/xmax
	    xval1=xval(jk)*xcm/xmax
            write(11,*) 'stroke'
            call clipin (xval1,0.d0,xval2,ycm)
	    goto 601
         endif
 600     continue
         xval1=xval(jk)*xcm/xmax
         write(11,*) 'stroke'
         call clipin (xval1,0.d0,xcm,ycm)
 601     continue
      endif
!........order of the lines
       do 608 je2=1,n_ene(jk)
 608      iste(je2)=je2
!
       if(nikir.ne.0.and.jk.lt.n_kpt) then
         if(n_ene(jk  ).ne.neik(jk  )) stop 'no of egval 1'
         if(n_ene(jk+1).ne.neik(jk+1)) stop 'no of egval 2'
         if(k_name(jk  )(1:1).eq.' '.or. &
            k_name(jk+1)(1:1).eq.' ')  call bndind(jk,nevl)
       endif
!
!.....for each energy state
      do 260 je=1,n_ene(jk)
!
!........define colors
         col='0.0 setgray'
         if(iprto(2).ne.0) then         
           if(iprto(2).eq.1) then
              nprtc = 1
              write(col,6001) nprtc
           elseif(iprto(2).eq.2) then
              nprtc = 2 + mod(iscol2(je),3)
              write(col,6001) nprtc
           elseif(iprto(2).eq.3) then
              nprtc = 10 + mod(iscol2(je),12)
              write(col,6010) nprtc
           elseif(iprto(2).eq.4) then
              nprtc = 10 + mod(iscol2(je),12)
              write(col,6010) nprtc
           else
              stop 'color switch > 4 not implemented'
           endif
           if(nprtc.gt.22) stop 'spag: need more color definitions'  
!
         endif  
 !        write(11,*) '0.0 ',col  
!
 6001    format('c0',I1,8X)
 6010    format('c', I2,8X)

         xval1=xval(jk)*xcm/xmax
         yval1=ycm/(ymax-ymin)*(eigen(je,jk)-ymin)
         if (yval1.gt.ycm) goto 260
         if (yval1.lt.0.d0) goto 260
         write(11,*) col  
         if(iprto(1).ne.1) then
	 call transt(xval1,yval1)
           if (charr(je,jk).le.char0) then
             call pointi (0,0.01d0,iprto)
           else
	     call pointi (1,charr(je,jk),iprto)
           endif
	 call transt(-xval1,-yval1)
         endif
!........draw lines of the energy bands
         if(k_name(jk  )(1:1).eq.' '.or. &
            k_name(jk+1)(1:1).eq.' ') then
         if((iprto(1).ne.0).and.(jk.ne.n_kpt).and. &
            (je.le.n_ene(jk+1)).and.iste(je).ne.0 ) then

           xval11=xval(jk  )*xcm/xmax
           yval11=ycm/(ymax-ymin)*(eigen(je,jk)-ymin) 
           xval22=xval(jk+1)*xcm/xmax
           yval22=ycm/(ymax-ymin)*(eigen(iste(je),jk+1)-ymin) 
	   call writln(xval11,yval11,xval22,yval22)    
         endif    
         endif
 260     continue
!
         if(iprto(2).ne.4) then
         do 609 je2=1,n_ene(jk)
 609       iscol2(iste(je2))=iscol(je2)
         do 610 je2=1,n_ene(jk)
 610       iscol(je2)=iscol2(je2)
         endif
!
 250  continue
!
!..... write energies bandwise on file
!
      write(6,*) 'writing energy-bands on file'

!...xmgrace: write header for data to case.bands.agr
      write(40,*) "#    k       ene     character"
!...xmgrace end   
      do 270 je=nb_min,min(nu_min,nb_max)
      write(10,*) ' bandindex:',je  
!...xmgrace: write command for autoscaling to case.bands.agr
        write(40,*) '@ autoscale onread none'
!...xmgrace end   
        write(dat_dummy,'(i4)') je
        read(dat_dummy,'(A)') string_dummy
        write(40,*) '@ target g0.s'//adjustl(string_dummy)
        write(40,*) '@ type xysize'
      if(jatom_list(1).gt.0) then
        write(40,*) '@ s'//adjustl(string_dummy), 'symbol 1'
      endif  
      write(40,225) je
!...xmgrace end   
      do 280 jk=1,n_kpt
!...xmgrace: write bands to case.bands.agr
      write(40,226) xval(jk),eigen(je,jk), charr(je,jk)*7
 225  format(/,' # bandindex:',1i3)
 226  format(3f10.5)   
 227  format(5f10.5)
!...xmgrace end   

 280  write(10,227)vk(1,jk),vk(2,jk),vk(3,jk),xval(jk),eigen(je,jk)
! 280  write(10,227) xval(jk),eigen(je,jk)
      write(40,'(A)') "&"
 270  continue 

  
  
!
      call psend
      write(6,*) ' plotcl done'
      stop 'SPAGH END'
!
!.....ERROR HANDLING FOR FILES
!
! 7    write(6,*) 'cannot open LAPW1 output file'
!      write(6,*) 'filename=',ofname
!      stop 'end NOT ok'
!
! 20   write(6,*) 'cannot open STRUCT file'
!      write(6,*) 'filename=',sfname
!      stop 'end NOT ok'
!
! 8    write(6,*) 'cannot open VIEWING-parameter file'
!      write(6,*) 'filename=',ifname
!      stop 'end NOT ok'
 140  format(f5.2,' setlinewidth')
 900 STOP 'GTFNAM - Exactly one commandline argument has to be given.'
      end

