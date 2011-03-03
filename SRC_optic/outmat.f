      SUBROUTINE OUTMAT(OUTME)
      use opme
      use bindex
      use xa3
!ad
!ad ___________________ CALCULATION OF MATRIX ELEMENTS _________________
!ad

      INCLUDE 'param.inc'
      IMPLICIT REAL*8 (A-H,O-Z) 
!ole ##### Begin #####
      LOGICAL     INFO_FLAG
      LOGICAL     REL,IM,LSO,SPIN,MME_FLAG
!ole ##### End #####
      CHARACTER*3  OUTME
      CHARACTER*10 KNAME
      COMPLEX*16  O(3),OX1,OX2,OY1,OY2,OZ2,OZ1     
!      
      COMMON /outop/ Ncol,icol(9)          
      COMMON /KPOI/ S,T,Z,NEMIN,NEMAX,KKZ,N,NNLO,KNAME
      COMMON /LEAD/ KFIRST,KLAST,KEND,KSTEP,KOUT,KSTOP
      COMMON /COM/  EMIN,EMAX,ELECN,EULIMIT,EDLIMIT, &
                    NK,IOUT,NSPIN,NAT,NBAND,ix,NB(NKPT),MINWAV,MAXWAV
!      COMMON /OPME/ OPMATX(NUMEO),OPMATY(NUMEO),OPMATZ(NUMEO)  
      COMMON /MIM / MIMA(2) 
!ole ##### Begin #####
      COMMON /CLOGIC/ LSO,SPIN,REL,MME_FLAG
!ole ##### End #####
!      COMMON /BINDEX/ N_(numeo),NN_(numeo),NIN(NUME,NUME)

!ad   write(6,*)'lso,spin:  ',LSO,SPIN 
      NEMIN=MIMA(1)
      NEMAX=MIMA(2)
!ole ##### Begin #####
      INFO_FLAG=.FALSE.
      IF (MME_FLAG) THEN
         WRITE(24,9010) NK,NEMIN,NEMAX,EMIN,EMAX,KNAME
      ELSE
         WRITE(3,9010) NK,NEMIN,NEMAX,EMIN,EMAX,KNAME
      if(OUTME.EQ.'ON ') then
        WRITE(4,9010) NK,NEMIN,NEMAX,EMIN,EMAX,KNAME
      endif
         WRITE(9,9010) NK,NEMIN,NEMAX,EMIN,EMAX,KNAME
      END IF
!ole ##### End #####

      NBINDEX=0
      IF (LSO.AND.(.NOT.SPIN)) THEN
!ole ##### Begin #####
         IF (MME_FLAG) THEN
            write(6,*) 'SPIN-orbit coupling for systems without inversion requires spinpolarized setup. Contact authors'
            stop 'SO requires spinpol.setup'
            INFO_FLAG=.TRUE.
            DO  NB1=NEMIN,NEMAX
               DO  NB2=NB1,NEMAX
                  N1=NIN(NB1,NB2)
                  IF (NB2.NE.NB1+1) THEN
                     WRITE(24,9030) NB1,NB2, & 
                          DBLE(OPMATX(N1)),AIMAG(OPMATX(N1)), &
                          DBLE(OPMATY(N1)),AIMAG(OPMATY(N1)), &
                          DBLE(OPMATZ(N1)),AIMAG(OPMATZ(N1))  
                  ELSE
!ole case of the Selection rule
                     WRITE(24,9030) NB1,NB2,0.,0.,0.,0.,0.,0.
                  END IF
               END DO
            END DO
         ELSE
!ole ##### End #####  
      DO  NB1=NEMIN,NEMAX,2
        DO  NB2=NB1,NEMAX,2
          if(nb2.gt.nb1) then
          N1=NIN(NB1,NB2)
          N2=NIN(NB1+1,NB2+1)
          N3=NIN(NB1,NB2+1)
          N4=NIN(NB1+1,NB2)
          OX1=(OPMATX(n1)+conjg(OPMATX(n2)))
          OY1=(OPMATY(n1)+conjg(OPMATY(n2)))
          OZ1=(OPMATZ(n1)+conjg(OPMATZ(n2)))
          OX2=(OPMATX(n3)-conjg(OPMATX(n4)))
          OY2=(OPMATY(n3)-conjg(OPMATY(n4)))
          OZ2=(OPMATZ(n3)-conjg(OPMATZ(n4)))
!ad
          OPMATX(n1)=OX1
          OPMATX(n2)=OX1
          OPMATX(n3)=OX2
          OPMATX(n4)=OX2
          OPMATY(n1)=OY1
          OPMATY(n2)=OY1
          OPMATY(n3)=OY2
          OPMATY(n4)=OY2
          OPMATZ(n1)=OZ1
          OPMATZ(n2)=OZ1
          OPMATZ(n3)=OZ2
          OPMATZ(n4)=OZ2
!ad
          else
!ad
          N1=NIN(NB1,NB2)
          N2=NIN(NB1+1,NB2+1)
          N3=NIN(NB1,NB2+1)
          OX1=OPMATX(n1)+OPMATX(n2)
          OY1=OPMATY(n1)+OPMATY(n2)
          OZ1=OPMATZ(n1)+OPMATZ(n2)
!ad
          OPMATX(n1)=OX1
          OPMATX(n2)=OX1
          OPMATX(n3)=0.0d0
          OPMATY(n1)=OY1
          OPMATY(n2)=OY1
          OPMATY(n3)=0.0d0
          OPMATZ(n1)=OZ1
          OPMATZ(n2)=OZ1
          OPMATZ(n3)=0.0d0

          endif
          END DO
          END DO
         END IF
      END IF

       DO 119 NB1=NEMIN,NEMAX
        DO 119 NB2=NB1,NEMAX 
          NBINDEX=NBINDEX+1  
          O(1)=(OPMATX(nbindex))
          O(2)=(OPMATY(nbindex))
          O(3)=(OPMATZ(nbindex))
!ad
!ad   write out complex matrix elelemts before symmetrization
!ad
      if(OUTME.EQ.'ON ') write(4,9040) NB1,NB2,O(1),O(2),O(3),E(NB2)-E(NB1)

!ole ##### Begin #####
!ole          write(4,9030) NB1,NB2,O(1),O(2),O(3)

            IF ((MME_FLAG).AND.(.NOT.INFO_FLAG)) THEN
               WRITE(24,9030) NB1,NB2, &
                    (DBLE(O(i)),AIMAG(O(i)),i=1,Ncol/2)
            END IF
 
            IF (.NOT.MME_FLAG) THEN
!ole #####  End  #####
!ad
          call outsym(o,nb1,nb2)
!ole ##### Begin #####
            END IF 
!ole #####  End  #####
  119 CONTINUE 
      RETURN
 9010 FORMAT(/,2X,' KP:',I6,' NEMIN NEMAX : ',2I5, &
            ' dE:',2f5.2,' K:',a10 /)
 9030 FORMAT(3X,2I4,6E13.6)
 9040 FORMAT(3X,2I4,6E13.6,F13.8)

      END

      SUBROUTINE OUTSYM(O,nb1,nb2)
!ad
!ad  ________________ SYMMETRIZATION OF MATRIX ELEMENTS ________________                  
!ad
!ad   May 1999: updated by Jan Kunes
!ad

      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      COMMON /outop/  Ncol,icol(9)
      COMMON /SYM2/   TAU(3,NSYM),IORD,IMAT(3,3,NSYM)
      COMMON /SYMo/   opimat(3,3,NSYM)
      COMMON /SYMd/   det(NSYM)
      COMPLEX*16  O(3),o2(6),o1(3)
!      COMPLEX*16  iimag
      DIMENSION outm(9)
!ad
!        iimag=(0.0,1.0)
        syfac=1.D0/REAL(iord)
        do ii=1,6
        o2(ii)=0.0
        end do
!        o3=((cdabs(O(1)+iimag*O(2)))**2- &
!            (cdabs(O(1)-iimag*O(2)))**2)/4.0
!ad
        do 772 l=1,iord
           do ii=1,3
           o1(ii)=0.0
           end do
         do   i=1,3
         do  ii=1,3
         o1(i)=o1(i)+opIMAT(ii,i,l)*O(ii)
         end do
         end do
           o2(1)=o2(1)+o1(1)*conjg(o1(1))
           o2(2)=o2(2)+o1(2)*conjg(o1(2))
           o2(3)=o2(3)+o1(3)*conjg(o1(3))
         if (det(l).LT.0) then
           o2(4)=o2(4)+conjg(o1(2))*o1(1)
           o2(5)=o2(5)+conjg(o1(3))*o1(1)
           o2(6)=o2(6)+conjg(o1(3))*o1(2)
                          else
           o2(4)=o2(4)+o1(2)*conjg(o1(1))
           o2(5)=o2(5)+o1(3)*conjg(o1(1))
           o2(6)=o2(6)+o1(3)*conjg(o1(2))
         end if
  772   continue

!____________________________________________________

	outm(1)=Dble(o2(1))*syfac
	outm(2)=Dble(o2(2))*syfac
	outm(3)=Dble(o2(3))*syfac        
        outm(4)=Dble(o2(4))*syfac
	outm(5)=Dble(o2(5))*syfac
	outm(6)=Dble(o2(6))*syfac
	outm(7)=aIMAG(o2(4))*syfac
	outm(8)=aIMAG(o2(5))*syfac
	outm(9)=aIMAG(o2(6))*syfac
!____________________________________________________

!cad
        if(nb1.eq.nb2) then
        write(9,9030)  NB1,NB2,(outm(icol(i)),i=1,Ncol)
        endif
!cad
        WRITE(3,9030)  NB1,NB2,(outm(icol(i)),i=1,Ncol)
!cad
 9030 FORMAT(3X,2I4,9E13.6)
 23   continue
      END
