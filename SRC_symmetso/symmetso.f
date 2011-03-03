       program symmetso
!  Modified program SYMMETRY of WIEN97 - analysis of symmetry in the
!  presence of spin-orbit coupling and magnetization. Program creates
!  modified structure file case.struct_so.
!  Analysis is made in two steps:
!  1/ File case.struct is read, check is made which operations are
!     retained with M in (tm1,tm2,tm3) direction. Distribution of atoms
!     into groups of magnetically equivalent atoms is performed and
!     intermediate case.struct file created (SUBROUTINE SYMSO).
!  2/ For each nonequivalent group of atoms local rotation matrix
!     local symmetry and nonzero (l,m) components are found - output
!     is then analogous as when using SYMMETRY package.
!
! ----------------------------------------------------------------------
!       gives the symmetry operations of all atoms of the unit cell      
!                                  
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      CHARACTER*4       LATTIC,type 
      CHARACTER*11      STATUS,FORM                                     
      CHARACTER*79      TITLE
      character*69      help2 
      CHARACTER*80      DEFFN, ERRFN ,fname                                   
      character*20      oname,onami,onamj,onaml
!                                                                       
!-----------------------------------------------------------------------
      CHARACTER*79, allocatable ::  name(:)
      real*8, allocatable :: pos(:,:),spos(:,:)
      integer,allocatable :: mult(:),iatom(:),index1(:),index2(:)
!       dimension pos(3,NDIF),index1(nato),index2(nato),spos(3,NDIF)
!       dimension name(nato),mult(nato),iatom(NDIF)
       dimension tm(3),iop(noper),a(3),atrans(3,NOPER)
       dimension alpha(3),rotold(3,3),amat(3,3,NOPER)
       dimension i0(noper),ij(noper,noper),indso(noper)
       logical numsym(noper)
      COMMON /GENER  /  BR2(3,3)                                     
      integer,allocatable :: iatnr(:)
!      common /iat/ iatnr(NATO)
!-----------------------------------------------------------------------
!                                                                       
      CALL GTFNAM(DEFFN,ERRFN)
      OPEN(1,FILE=DEFFN,STATUS='OLD',ERR=8000)
 8003 READ(1,*,END=8001) IUNIT,FNAME,STATUS,FORM,IRECL
      OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM, &
      ERR=8002)
      GOTO 8003
 8000 WRITE(*,*) ' ERROR IN OPENING symmetso.def !!!!'
      STOP 'NN.DEF'
 8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
      WRITE(*,*) '       FILENAME: ',FNAME,'  STATUS: ',STATUS, &
      '  FORM:',FORM
      STOP 'OPEN FAILED'                                
 8001 CONTINUE
      do 39 i=1,noper
39    i0(i)=i
!  read direction of magnetization in units of unit cell parameters
!  if ipot.lt.0 only struct_so and struct_interm are written
!  if ipot.gt.0 VNS_so and VSP_so files are also written
!cccc....use inso file
      read(5,*)
      read(5,*)
      read(5,*)
      read(5,*)(tm(i),i=1,3)
       ipot=2
!ccc....
      write(6,100)(tm(i),i=1,3)
100   format(' Magnetization is along',3f8.4,' direction')
      ifile=22
      call latsym(nsym,amat,atrans,ifile,nat)
      write(6,*)' LATSYM done'
      rewind 22
      call symso(tm,ipot,nA,nB,invers,indso)
      rewind 24
      write(6,*)' SYMSO done'
!
!  determine local symmetry of atom types in presence of s-o and M
      ifile=24
!.....check lattice and generate symmetry operations
      call latsym(nsym,amat,atrans,ifile,nat)
      write(6,*)' LATSYM done'
      allocate (iatnr(NAT),spos(3,nat*48*16))
      allocate ( pos(3,nat*48*16),index1(nat),index2(nat))
      allocate ( name(nat),mult(nat),iatom(nat*48*16))

!.....read new struct file
      rewind 24
      call rstruc(title,lattic,nat,a,alpha,index2,iatom, &
      index1,pos,mult,name,iatnr,tm)
!     tm for R lattice has now been changed from R to hex coordinates
      write(6,*)' RSTRUCT done'
      rewind 24
!.....START READING/writing  FILE STRUCT (first 4 lines)                     
      do 5 i=1,4
      READ(24,1510) TITLE 
 5    write(21,1510) TITLE 
 1510 FORMAT(A79)                                                       
      noper1=noper
      if(lattic(1:1).eq.'H') noper1=24
      if(lattic(1:1).eq.'R') noper1=24
      nbas=index2(nat)
      write(6,199) nbas
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
            do 23 i2=1,nbas 
               ja=iatom(i2)
               x = pos(1,i2)-pos(1,i1)
               y = pos(2,i2)-pos(2,i1)
               z = pos(3,i2)-pos(3,i1)
!      if(i1.eq.1.and.io.eq.1) write(6,'(i3,5x,3f10.6,5x,3f10.6,5x,3f10.6)') i2,(pos(j,i1),j=1,3),(pos(j,i2),j=1,3),x/a(1), &
!          Y/a(2),z/a(3)
               iz=0
               nz=1
 124           call symop(io,x,y,z,x1,y1,z1,oname,nz,lattic)
               iz=iz+1
               x1 = x1+pos(1,i1)
               y1 = y1+pos(2,i1)
               z1 = z1+pos(3,i1)
!      if(jatom.eq.1.and.io.eq.2) write(6,*) i1,io,oname,i2,x1,y1,z1
!     *          ,x1/a(1),Y1/a(2),z1/a(3)
 1001          continue
               if(x1.lt.-0.0001) then
                     x1=x1+1.*a(1)*sin(alpha(3)/180.d0*pi)
                     y1=y1+1.*a(1)*cos(alpha(3)/180.d0*pi)
                     goto 1001
               endif
 1002          continue
               if(y1.lt.-0.0001) then
                     y1=y1+1.*a(2)*sin(alpha(1)/180.d0*pi)
                     z1=z1+1.*a(2)*cos(alpha(1)/180.d0*pi)
! cad                     y1=y1+1.*a(2)
                     goto 1002
               endif
 1003          continue
               if(z1.lt.-0.0001) then
                     z1=z1+1.*a(3)*sin(alpha(2)/180.d0*pi)
                     x1=x1+1.*a(3)*cos(alpha(2)/180.d0*pi)
! cad                     z1=z1+1.*a(3)
                     goto 1003
               endif
 1004          continue
               if(x1.gt. a(1)-0.0001) then
                     x1=x1-1.*a(1)*sin(alpha(3)/180.d0*pi)
                     y1=y1-1.*a(1)*cos(alpha(3)/180.d0*pi)
                     goto 1004
               endif
 1005          continue
               if(y1.gt. a(2)-0.0001) then
                     y1=y1-1.*a(2)*sin(alpha(1)/180.d0*pi)
                     z1=z1-1.*a(2)*cos(alpha(1)/180.d0*pi)
! cad                     y1=y1-1.*a(2)
                     goto 1005
               endif
 1006          continue
               if(z1.gt. a(3)-0.0001) then
                     z1=z1-1.*a(3)*sin(alpha(2)/180.d0*pi)
                     x1=x1-1.*a(3)*cos(alpha(2)/180.d0*pi)
! cad                     z1=z1-1.*a(3)
                     goto 1006
               endif
               do 24 itest=index1(ja),index2(ja)
                  x2=x1-pos(1,itest)
                  if(abs(x2).gt.tol) goto 24
                  y2=y1-pos(2,itest)
                  if(abs(y2).gt.tol) goto 24
                  z2=z1-pos(3,itest)
                  if(abs(z2).gt.tol) goto 24
                  goto 123
  24           continue
               goto 22  
 123           if(iz.lt.nz) goto 124
  23        continue 
! find whether the operation is compatible with M direction
         call symop(io,tm(1),tm(2),tm(3),tx,ty,tz,oname,nz,lattic)
         dp=abs(tm(1)-tx)+abs(tm(2)-ty)+abs(tm(3)-tz)
         dm=abs(tm(1)+tx)+abs(tm(2)+ty)+abs(tm(3)+tz)
         if((dp.lt.tol).or.(dm.lt.tol))then
            isym=isym+1
            numsym(io)=.true.
            iop(isym)=io
            write(6,299) name(jatom)(1:10),isym,io,oname
         endif
  22     continue
       write(6,*)' check whether the operations form a group'
! general position point (xg,yg,zg)
         xg=0.151
         yg=0.337
         zg=0.478
         ngr=1
         do 30 i=1,isym
         do 30 j=1,isym
         call symop(iop(i),xg,yg,zg,xi,yi,zi,onami,ni,lattic)
         call symop(iop(j),xi,yi,zi,xj,yj,zj,onamj,nj,lattic)
         ng=0
         do 31 l=1,isym
         call symop(iop(l),xg,yg,zg,xl,yl,zl,onaml,nl,lattic)
         diff=(xl-xj)**2+(yl-yj)**2+(zl-zj)**2
         if(diff.lt.0.00001)then
         ij(i,j)=l
!        write(6,201)i,j,l
201      format(' G',i2,' x  G',i2,' =  G',i2)
         ng=1
         endif
31       continue
         if(ng.eq.0)then
         ngr=0
         endif
30       continue
         if(ngr.eq.0)then
         write(6,*)'!! SYMM. OP. DO NOT FORM A GROUP !!'
         else
         write(6,*)' SYMMETRY OPERATIONS FORM A GROUP' &
                  ,' ,GROUP MULTIPLICATION TABLE:'
         endif
         write(6,203)(i0(i),i=1,isym)
         write(6,*)'j'
         do 32 i=1,isym
         write(6,204)i,(ij(i,j),j=1,isym)
32       continue
203      format('   i',24i3)
204      format(i2,2x,24i3)
!
!
!   read/write struct file
         READ(24,1011) IATNR1,( SPOS(J,1),J=1,3 )
         READ(24,1013) mult1,isplit
 1011    FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
         DO 55 M=2,MULT1                                     
 55      READ(24,1011) IATNR1,( SPOS(J,M),J=1,3)
         READ(24,1510) title 
         read(24,1012) ((rotold(i,j),j=1,3),i=1,3)          
 1012    format(20x,3f10.8)
 1013    format(15x,i2,17x,i2)
 1014    format(i4,6x,'NUMBER OF SYMMETRY OPERATIONS')
 1015    FORMAT(3I2,F10.7)                                             
 1016    format(i8)
 1017    format(i8,a4,i4)
 1018    format(i8,a4,i4,' so. oper.  type  orig. index')
!
!
      call class(jatom,isym,numsym,lattic,spos,mult1,isplit, &
      title,rotold,iatnr)
  21  continue   
!
!....write struct_st 
      READ(24,*) IORD                                                
        WRITE(21,1014) IORD                                                
        DO  I=1,IORD   
          if(i.le.nA)then
           type='   A'
          else
           type='   B'
          endif
          do j=1,3                                                 
          READ(24,1015) I1,I2,I3,T1               
          write(21,1015) I1,I2,I3,T1               
          enddo
        read(24,*)
        if(i.eq.1)then
        write(21,1018) i,type,indso(i)
        else
        write(21,1017) i,type,indso(i)
        endif
        enddo
!....if system has no inversion center, write struct.kgen
      if(invers.eq.0)then
      rewind (21)
      call rewr(21,23)
      READ(21,*) IORD                                                
        WRITE(23,1014) IORD                                                
        DO  I=1,IORD   
          do j=1,3                                                 
          READ(21,1015) I1,I2,I3,T1               
          if(i.le.nA)then
          write(23,1015) I1,I2,I3,T1               
          else
          write(23,1015) -I1,-I2,-I3,T1               
          endif
          enddo
        read(21,*)
        write(23,1016) i
        enddo
      endif
!
      if(ipot.ge.0)then
      read(34,778,err=779)help2
      write(54,778)help2
      read(34,778,end=779)help2
      write(54,778)help2
779   continue
      endif
!
      stop
 778  format(a69)
 199  format('number of atoms: ',i3)                           
 299  format(a10,' G',i2,' oper. #',i3,5x,a20,' GM=',i2,'M')
      end 
