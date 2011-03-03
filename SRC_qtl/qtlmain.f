      PROGRAM QTLMAIN                                                      
! Program QTL calculates atomic-site projected densities of states. It is based
! on modified LAPW2U program for calculating density matrices. The program in fact
! calculates energy spectrum of density matrix for given atom and orbital L 
! (including terms non-diagonal in spin). Knowledge of such a density matrix allows
! for calculation of the densties of states of any symmetry within the given orbital
! including relativistic ones (for example eg, t2g symmetry or px,py,pz symmetry or
! (j,jz) symmetry). Those densities of states are obtained as diagonal elements 
! of the density matrix in appropriate basis. The density matrix is originally 
! calculated in Y(l,m)f(s) basis ( in XSPLT ) and then transformed to the appropriate
! basis by unitary tranformation ( in QTL ). The matrix of the unitary transformation
! is supplied from outside ( read in READC ). Program in its current version is
! designed for calculations with s-o coupling ( the prime aim of the program is to 
! calculate the realtivistic densities of states ).  
! Variable MODUS has values FULL or SUMA. If FULL projected densities of all symmetries
! are calculated. If SUMA only sums over given densities are on the output
! (i.e. total j=1/2,j=3/2 densities of states). The states over which to average
! are marked line * in the *.cf file. One can speed up the calculation considerably
! when calculating (j,jz) densities. In this case symmetrization is not necessary
! (due to transformation properties of spherical harmonics). Attention without 
! symmetrization only the diagonal elements of the density matrix (which are the
! densities of states) are correct. The symmetrization can be switched of by 
! setting NOSYM in input file.
      USE param
      USE struct
      USE case
      USE sym2
      USE com
      USE abc
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*4      MODUS,adum
      CHARACTER*5      MODSYM                                       
      CHARACTER*10    KNAME
      CHARACTER*11     STATUS,FORM                                      
      CHARACTER*67       ERRMSG
      CHARACTER*80       DEFFN, ERRFN
      CHARACTER*80     FNAME,FNAME1,enefn,enefnd  
      CHARACTER*80     VECFN
      INTEGER          IPROC

      LOGICAL         LOG_SYM,SO,cmplx                               
      COMMON /CHAR/   MODUS
      COMMON /GENER/  BR1(3,3),BR2(3,3)                                 
      COMMON /RADFU/  RRAD1(NRAD,0:LMAX2),RADE1(NRAD,0:LMAX2), &
                      RRAD2(NRAD,0:LMAX2),RADE2(NRAD,0:LMAX2)                  
      COMMON /PROC/    VECFN(2)
      COMMON /IPROC/   IPROC 
                                                                       
      DIMENSION  S(3)
      
      DATA LOG_SYM /.true./, SO /.false./,cmplx /.false./
!-----------------------------------------------------------------------  
      CALL GTFNAM(DEFFN,ERRFN,IPROC)
      CALL ERRFLG(ERRFN,'Error in LDAU')
      OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)
               iso=1
   10 CONTINUE
         READ (1,*,END=20,ERR=960) IUNIT,FNAME,STATUS,FORM,IRECL
         OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
         if(iunit.eq.9) VECFN(1)=FNAME
         if(iunit.eq.10)  VECFN(2)=FNAME
         if(iunit==59) enefnd=fname
         if(iunit==60) enefn=fname
         if(iunit.eq.7) then
         read(iunit,123,end=12)adum
 123     format(A5)
         cmplx=.true.
 12      continue
         end if
	 if(iunit.eq.9) then
           do i=80,1,-1
             if(fname(i:i).ne.' ') then
               if(fname(i-3:i).eq.'soup') then
               SO=.true.
               iso=2
               cmplx=.true.
	       endif
             goto 10
             endif
           enddo
         endif
      GOTO 10
   20 CONTINUE
      CLOSE (1)
      
      if (cmplx) write(6,*)'CALCUATION WITH COMPLEX VECTORS'
      if (so) write(6,*)'CALCULATION WITH SOC-UP AND DN VECTORS REQUIRED'

	write(6,*)'QTLMAIN***********************************'
      IF(IPROC.GT.0)THEN
         write(6,*)'Running QTL on ',IPROC,' processors'
         write(6,*)' '
	 do is=1,2
	 call mknam(FNAME,VECFN(is),1)
	 iunit=8+is
	 status='UNKNOWN'
	 form='UNFORMATTED'
	 OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
         end do
      ELSE
         write(6,*)'Running QTL in single processor mode'
         write(6,*)' '
      END IF

!.....READ STRUCT 
      CALL init_struct                                                      
                                                                       
!.....READ INPUT AND POTE  
      READ(5,1003)  MODUS 
      write(6,*)  ' Modus: ', MODUS
!     if (modus.eq.'SPIN'.AND.not.so) STOP 'this modus allowed with 
!    &    SOC only'
      READ(5,1004) MODSYM
      if (MODSYM.EQ.'NOSYM') then
      LOG_SYM=.false.
      iord=1
      write(6,*)'Symmetrization over eq. k-points is not performed'
      write(6,*)'allowed for rotationally invariant densities of states'
      end if
      write(6,*)'Atoms for projected density of states calculation:'
      READ(5,*)  EMIN,EMAX
      READ(5,*)  EF
      CALL init_coef(ndim2)
      CALL init_sym2(iord)
 1234 FORMAT(//,1A)
      WRITE(6,800)                                                      
      WRITE(6,805)  TITLE                                               
      WRITE(6,810)  LATTIC                                              
      WRITE(6,820)  AA,BB,CC                                            
      WRITE(6,840)  NAT                                                 
      WRITE(6,850)  IREL                                                
      CALL LATGEN   
	write(6,*)'IORD=',IORD

! if log_sym.eq.(.false.) then  the knowledge of magnetization direction
! is not necessary, but IR part of BZ must be of course correctly chosen
!....reading *.inso
	if (so) then
      do i=1,3
      read(4,*)
      end do
      read(4,*)s(1),s(2),s(3)
      CALL ANGLE(S,THETA,FI)   
      write(6,*)'direction of magnetic moment theta, phi:', theta,fi 
        end if
      write(6,*)'LOG_SYM:',LOG_SYM
      IF (LOG_SYM) THEN 
      if (so) CALL SYM(THETA,FI)
      ELSE
	do i=1,3
	 do j=1,3
	  iz(i,j,1)=0
	  if (i.eq.j) iz(i,j,1)=1
	 end do
        end do
!	idet(1)=1
       END IF

      if (MODUS.EQ.'FULL'.OR.MODUS.EQ.'SUMA') CALL READC(MODUS,iso)
!find nkpt, nmat and nume in energy file
     k=0
      iloop=0
 778  continue
      if(IPROC.GT.0) then
         iloop=iloop+1
         close(59)
         close(60)
         call mknam(FNAME,enefn,ILOOP)
         open(60,FILE=FNAME,STATUS='old',ERR=951)
         call mknam(FNAME,enefnd,ILOOP)
         open(59,FILE=FNAME,STATUS='unknown',ERR=951)
      endif

      DO I=1,NAT
         READ(60,'(f9.5)') EMIST
         READ(60,'(f9.5)') EMIST
      ENDDO
      DO I=1,NAT
         READ(59,'(f9.5)',END=1005) EMIST
         READ(59,'(f9.5)',END=1005) EMIST
      ENDDO
      JSPIN=2
 1005 CONTINUE
       DO
         READ(60,'(3e19.12,a10,2i6)',IOSTAT=ios) SS,T,Z,KNAME,N,NEn
         IF (ios /= 0) EXIT
         k=k+1
         nmat=MAX(n,nmat)
         nume=MAX(nen,nume)
         DO ii=1,nen
            READ(60,*) NUM,E1
         ENDDO
         IF(jspin.EQ.2) THEN
            READ(59,'(3e19.12,a10,2i6)',IOSTAT=ios) SS,T,Z,KNAME,N,NEn
            nmat=MAX(n,nmat)
            nume=MAX(nen,nume)
            DO ii=1,nen
               READ(59,*) NUM,E1
            ENDDO
         ENDIF
      ENDDO
      IF(ILOOP.LT.IPROC) goto 778
      nkpt=k
      REWIND(59)
      REWIND(60)
      CALL init_com(nkpt)
      CALL init_abc(nume,nmat,lmax2,ndim2,nrf)
      itape=10
      jtape=18
      call l2main(cmplx)
      CALL CPUTIM(TTIME)
                                                                       
!.....CALCULATE CPUTIME REQUIRED                                        
      WRITE(6,2000)                                                     
      WRITE(6,2010) TTIME                                        
      CALL ERRCLR(ERRFN)
      STOP ' QTL END'                                                 
!
!        error handling
!
  910 INFO = 1
!
!        'lapw2.def' couldn't be opened
!
      WRITE (ERRMSG,9000) FNAME
      CALL OUTERR('lapw2U',ERRMSG)
      GOTO 999
  920 INFO = 2
!
!        file FNAME couldn't be opened
!
      WRITE (ERRMSG,9010) IUNIT
      CALL OUTERR('lapw2U',ERRMSG)
      WRITE (ERRMSG,9020) FNAME
      CALL OUTERR('lapw2U',ERRMSG)
      WRITE (ERRMSG,9030) STATUS, FORM
      CALL OUTERR('lapw2U',ERRMSG)
      GOTO 999
  930 INFO = 3
!
!        illegal number of equivalent atoms
!
      CALL OUTERR('lapw2U','MULT .EQ. 0')
      GOTO 999
  951 INFO = 5
      write(6,*) 'error:',FNAME
      CALL OUTERR('lapw2U','file open error')
      GOTO 999
  955 INFO = 6
      CALL OUTERR('lapw2U','Too many atoms (NDIF too small).')
      GOTO 999
  956 INFO = 56
      CALL OUTERR('lapw2U','LXDOS must be 3 for ISPLIT=999.')
      GOTO 999
  960 INFO = 7
!
!        Error reading file 'lapw2.def'
!
      WRITE (ERRMSG,9040) FNAME
      CALL OUTERR('lapw2U',ERRMSG)
      GOTO 999
  999 STOP 'lapw2U - Error'
!                                                                       
!                                                                       
  43  FORMAT(3X,A77)                                                    
  44  FORMAT(I3,A77)                                                    
 700  FORMAT(I3,A77)                                                    
 800  FORMAT(////,30X,50(1H-),/,33X,'S T R U C T U R A L   ',            &
             'I N F O R M A T I O N',/,30X,50(1H-),//)                  
 805  FORMAT(3X,'SUBSTANCE',20X,'= ',A80,/)                             
 810  FORMAT(3X,'LATTICE',22X,'= ',A4)                                  
 820  FORMAT(3X,'LATTICE CONSTANTS ARE',8X,'= ',3F12.7)                 
 830  FORMAT(3X,'SYMMETRY ATTITUDE IS',9X,'= ',A4)                      
 840  FORMAT(3X,'NUMBER OF ATOMS IN UNITCELL  = ',I3)                   
 850  FORMAT(3X,'MODE OF CALCULATION IS',7X,'= ',A4)                    
 860  FORMAT(3X,'SELFCONSISTENT CYCLE-NUMBER  = ',I3,/)                 
 870  FORMAT(3X,'TYPE OF COORDINATES IN DSPLIT= ',A5)                   
 1000 FORMAT(A80)                                                       
 1010 FORMAT(A4,24X,I2,1x,a4,/,13X,A4,18X,A4)                                 
 1020 FORMAT(6F10.7,10X,F10.7)                                          
 1030 FORMAT(5X,I3,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
 1031 FORMAT(5X,I3,4X,F10.7,3X,F10.7,3X,F10.7)                          
 1040 FORMAT(///,3X,'ERROR IN lapw2U : MULT(JATOM)=0 ...',                &
             /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)           
 1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
 1002 FORMAT(3F10.5,I5)                                                 
 1003 FORMAT(A4)                                                        
 1004 FORMAT(A5)                                                        
 1060 FORMAT(//,':NOE  :',1X,'NUMBER OF ELECTRONS',10X,'= ',F7.3,            &
             //   ,':FER  :',1X,'F E R M I - ENERGY',11X,'= ',F9.5)            
 1061 FORMAT(//,':NOE  :',1X,'NUMBER OF ELECTRONS',10X,'= ',F7.3,            &
       //   ,':FER  :',1X,'F E R M I - ENERGY(TETRAH.M.)','= ',F9.5)           
 1062 FORMAT(//,':NOE  :',1X,'NUMBER OF ELECTRONS',10X,'= ',F7.3,            &
       //   ,':FER  :',1X,'F E R M I - ENERGY(GAUSS-.M.)','= ',F9.5)           
 1063 FORMAT(//,':NOE  :',1X,'NUMBER OF ELECTRONS',10X,'= ',F7.3,            &
       //   ,':FER  :',1X,'F E R M I - ENERGY(FERMI-SM.)','= ',F9.5)           
 2000 FORMAT(//,3X,'=====>>> CPU TIME SUMMARY',/)                       
 2010 FORMAT(12X,'TOTAL       : ',F8.1)       
 6000 FORMAT(///,3X,'ERROR IN lapw2U : MULT(JATOM)=0 ...', &
             /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
 9000 FORMAT('can''t open definition file ',A40)
 9010 FORMAT('can''t open unit: ',I2)
 9020 FORMAT('       filename: ',A50)
 9030 FORMAT('         status: ',A,'  form: ',A)
 9040 FORMAT('Error reading file: ',A47)
      END                                                               
