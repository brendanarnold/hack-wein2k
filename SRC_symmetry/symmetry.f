       program symmetry
!
! ----------------------------------------------------------------------
!       gives the symmetry operations of all atoms of the unit cell      
!                                  
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      CHARACTER*4       LATTIC                                          
      CHARACTER*10      KNAME                                     
      CHARACTER*11      STATUS,FORM                                     
      CHARACTER*79      TITLE                                            
      CHARACTER*80      DEFFN, ERRFN ,fname                                   
      character*20      oname
!                                                                       
!-----------------------------------------------------------------------
       real*8,allocatable ::  pos(:,:),spos(:,:)
       integer,allocatable :: index1(:),index2(:),mult(:),mult0(:),iatom(:)
       dimension a(3),atrans(3,NOPER)
       CHARACTER*79,allocatable :: name(:)
       dimension alpha(3),rotold(3,3),amat(3,3,NOPER)
       logical numsym(48)
      COMMON /GENER  /  BR2(3,3)                                     
      integer,allocatable :: iatnr(:)
      logical tric,trici
      common /tricl/ tric,trici
!-----------------------------------------------------------------------
!                                                                       
!      call getarg(2,fname)
!      if(fname.eq.'      ') call getarg(1,fname)
!       fname='symmetry.def'
      CALL GTFNAM(DEFFN,ERRFN)
      OPEN(1,FILE=DEFFN,STATUS='OLD',ERR=8000)
 8003 READ(1,*,END=8001) IUNIT,FNAME,STATUS,FORM,IRECL
      OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM, &
      ERR=8002)
      GOTO 8003
 8000 WRITE(*,*) ' ERROR IN OPENING symmetry.def !!!!'
      STOP 'NN.DEF'
 8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
      WRITE(*,*) '       FILENAME: ',FNAME,'  STATUS: ',STATUS, &
      '  FORM:',FORM
      STOP 'OPEN FAILED'                                
 8001 CONTINUE
!.....check lattice and generate symmetry operations
      call latsym(nsym,amat,atrans,nat)
      allocate (iatnr(NAT),spos(3,nat*48*16))
      allocate ( pos(3,nat*48*16),index1(nat),index2(nat))
      allocate ( name(nat),mult(nat),mult0(nat),iatom(nat*48*16))
!.....read struct file
      rewind 20
      call rstruc(title,lattic,nat,a,alpha,index2,iatom, &
      index1,pos,mult,mult0,name,iatnr)
!
      gmax=12.d0
      rewind 20
!.....START READING/writing  FILE STRUCT (first 4 lines)                     
      do 5 i=1,4
      READ(20,1510) TITLE 
 5    write(21,1510) TITLE 
 1510 FORMAT(A79)                                                       
!
!      read(17,*,end=778)
!      read(17,*,end=778)
!      read(17,*,end=778)
 778  continue
      noper1=noper
      if(lattic(1:1).eq.'H') noper1=24
      if(lattic(1:1).eq.'R') noper1=24
      nbas=index2(nat)
!      nbas=index
!           write(*,*) alpha
      write(6,199) nbas
!          do 777 ib=1,nbas
!777       write(6,*) ib,iatom(ib),(pos(k,ib),k=1,3)
      tol=1.e-3
      pi=4.d0*atan(1.d0)
!.....loop over atoms
      do 21 jatom=1,nat
         write(6,*) 
         write(6,*) 'ATOM:',iatnr(jatom)
         do 25 isym=1,48
 25      numsym(isym)=.false.
         isym=0
         i1=index1(jatom)
!.....loop over sym.operations
         do 22 io=1,noper1 
!..... dirty fix for triclinic lattic
         if(tric.and.(io.eq.1)) then
           oname='1'
           goto 822
         endif
         if(tric.and.(io.eq.2)) then
           if(mult0(jatom).eq.1.and.trici) then
             oname='-1'
             goto 822
           else
             goto 22
           endif
         endif
         if(tric.and.(io>2)) then   ! GM
            goto 22                 ! GM
         endif                      ! GM
!....
            do 23 i2=1,nbas 
               ja=iatom(i2)
               x = pos(1,i2)-pos(1,i1)
               y = pos(2,i2)-pos(2,i1)
               z = pos(3,i2)-pos(3,i1)
!     if(i1.eq.1.and.io.eq.1) write(6,*) i1,io,i2,x/a(1),Y/a(2),z/a(3)
               iz=0
               nz=1
 124           call symop(io,x,y,z,x1,y1,z1,oname,nz,lattic)
               iz=iz+1
               x1 = x1+pos(1,i1)
               y1 = y1+pos(2,i1)
               z1 = z1+pos(3,i1)
!     if(jatom.eq.1.and.io.eq.1) write(6,*) i1,io,oname,i2,x1,y1,z1,x1/a(1),Y1/a(2),z1/a(3)
 1001          continue
               if(x1.lt.-0.00001d0) then
!     if(jatom.eq.1.and.io.eq.1) write(*,*) 'alarm0',x1,y1,z1
                     x1=x1+1.d0*a(1)*sin(alpha(3)/180.d0*pi)
                     y1=y1+1.d0*a(1)*cos(alpha(3)/180.d0*pi)
                     goto 1001
               endif
 1002          continue
               if(y1.lt.-0.00001d0) then
!     if(jatom.eq.1.and.io.eq.1) write(*,*) 'alarm1',x1,y1,z1
                     y1=y1+1.d0*a(2)*sin(alpha(1)/180.d0*pi)
                     z1=z1+1.d0*a(2)*cos(alpha(1)/180.d0*pi)
! cad                     y1=y1+1.*a(2)
                     goto 1002
               endif
 1003          continue
               if(z1.lt.-0.00001d0) then
!     if(jatom.eq.1.and.io.eq.1) write(*,*) 'alarm2',x1,y1,z1
                     z1=z1+1.d0*a(3)*sin(alpha(2)/180.d0*pi)
                     x1=x1+1.d0*a(3)*cos(alpha(2)/180.d0*pi)
! cad                     z1=z1+1.*a(3)
                     goto 1001
               endif
 1004          continue
               if(x1.gt. a(1)*sin(alpha(3)/180.d0*pi)-0.00001d0) then
!               if(x1.gt. a(1)-0.0001) then
!     if(jatom.eq.1.and.io.eq.1) write(*,*) 'alarm3',x1,y1,z1
                     x1=x1-1.d0*a(1)*sin(alpha(3)/180.d0*pi)
                     y1=y1-1.d0*a(1)*cos(alpha(3)/180.d0*pi)
                     goto 1001
               endif
 1005          continue
               if(y1.gt. a(2)-0.00001d0) then
!     if(jatom.eq.1.and.io.eq.1) write(*,*) 'alarm4',x1,y1,z1
                     y1=y1-1.d0*a(2)*sin(alpha(1)/180.d0*pi)
                     z1=z1-1.d0*a(2)*cos(alpha(1)/180.d0*pi)
! cad                     y1=y1-1.*a(2)
                     goto 1002
               endif
 1006          continue
               if(z1.gt. a(3)-0.00001d0) then
!     if(jatom.eq.1.and.io.eq.1) write(*,*) 'alarm5',x1,y1,z1
                     z1=z1-1.d0*a(3)*sin(alpha(2)/180.d0*pi)
                     x1=x1-1.d0*a(3)*cos(alpha(2)/180.d0*pi)
! cad                     z1=z1-1.*a(3)
                     goto 1001
               endif
               do 24 itest=index1(ja),index2(ja)
!      if(jatom.eq.1.and.io.eq.1) write(6,*)x1,y1,z1,itest,pos(1,itest),pos(2,itest),pos(3,itest)
                  x2=x1-pos(1,itest)
!          write(6,*) io,jatom,i2,itest,x2  
                  if(abs(x2).gt.tol) goto 24
                  y2=y1-pos(2,itest)
!          write(6,*) io,jatom,i2,itest,y2  
                  if(abs(y2).gt.tol) goto 24
                  z2=z1-pos(3,itest)
!          write(6,*) io,jatom,i2,itest,z2  
                  if(abs(z2).gt.tol) goto 24
                  goto 123
  24           continue
               goto 22  
 123           continue
!               write(6,*) 'iz-check: iz.lt.nz',iz,nz
               if(iz.lt.nz) goto 124
  23        continue 
 822        write(6,299) name(jatom)(1:10),io,oname
            isym=isym+1
            numsym(io)=.true.
  22     continue
!
!
!   read/write struct file
         READ(20,1011) IATNR1,( SPOS(J,1),J=1,3 )
           if(spos(1,1).ge.1.d0) spos(1,1)=spos(1,1)-1.d0
           if(spos(2,1).ge.1.d0) spos(2,1)=spos(2,1)-1.d0
           if(spos(3,1).ge.0.999999999d0) spos(3,1)=spos(3,1)-1.d0
           if(spos(1,1).lt.0.d0) spos(1,1)=spos(1,1)+1.d0
           if(spos(2,1).lt.0.d0) spos(2,1)=spos(2,1)+1.d0
           if(spos(3,1).lt.0.d0) spos(3,1)=spos(3,1)+1.d0
         READ(20,1013) mult1,isplit
 1011    FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
         DO 55 M=2,MULT1                                     
         READ(20,1011) IATNR1,( SPOS(J,M),J=1,3)
           if(spos(1,m).ge.1.d0) spos(1,m)=spos(1,m)-1.d0
           if(spos(2,m).ge.1.d0) spos(2,m)=spos(2,m)-1.d0
           if(spos(3,m).ge.0.999999999d0) spos(3,m)=spos(3,m)-1.d0
           if(spos(1,m).lt.0.d0) spos(1,m)=spos(1,m)+1.d0
           if(spos(2,m).lt.0.d0) spos(2,m)=spos(2,m)+1.d0
           if(spos(3,m).lt.0.d0) spos(3,m)=spos(3,m)+1.d0
  55     continue 
         READ(20,1510) title 
         if(title(1:2).eq.'H ') gmax=20.
         read(20,1012) ((rotold(i,j),j=1,3),i=1,3)          
 1012    format(20x,3f10.8)
 1013    format(15x,i2,17x,i2)
 1014    format(i4,6x,'NUMBER OF SYMMETRY OPERATIONS')
 1015    FORMAT(3I2,F11.8)                                             
 1016    format(i8)
!
!
      call class(jatom,isym,numsym,lattic,spos,mult1,isplit, &
      title,rotold,iatnr)
         if(mult1*isym.ne.nsym) then
         write(6,*) '---------- ERROR ------------------'
         write(6,*) ' The (multiplicity of this atom)*(number of pointgroup-operations) is NOT'
         write(6,*) ' = (number of spacegroup-operations)'
         write(6,*) 'MULT:',mult1,' ISYM:',isym,' NSYM',nsym
         write(6,*) 'Check your struct file with      sgroup -wien case.struct'
         write(6,*) '---------- ERROR ------------------'
         endif
  21  continue   
!
!....write struct_st 
      READ(20,*) IORD                                                
      if(iord.eq.0) then
        WRITE(21,1014) NSYM                                                
        DO  I=1,NSYM   
             write(21,1015) (nint(amat(1,j,i)),j=1,3),atrans(1,i)
             write(21,1015) (nint(amat(2,j,i)),j=1,3),atrans(2,i)
             write(21,1015) (nint(amat(3,j,i)),j=1,3),atrans(3,i)
             write(21,'(i8)') i
        ENDDO
      else
        WRITE(21,1014) IORD                                                
        DO  I=1,IORD   
          do j=1,3                                                 
          READ(20,1015) I1,I2,I3,T1               
          write(21,1015) I1,I2,I3,T1               
          enddo
        read(20,*)
        write(21,1016) i
        enddo
      endif
!
      write(17,'(f6.2,10x,"GMAX")') gmax
      write(17,'(a)') 'NOFILE        FILE/NOFILE  write recprlist'
!
      stop
 199  format('number of atoms: ',i3)                           
 299  format(a10,1x,'operation #',i3,5x,a20)
      end 
