! modified: IIM
!........1.........2.........3.........4.........5.........6.........7
!234567890123456789012345678901234567890123456789012345678901234567890
!.......................................................................
      SUBROUTINE ARBDOS (NEMIN,NEMAX,NEOCC)
!     **   ADOS BUILDS FUNCTION A(K) AND SUMS OVER THE BANDS b,b`     **
!     **        AND OVER ALL TETRAEDERS                               **
!     **        AND READS THE TETRAEDER-POINTS FROM FILE 14           **
!     **                                                              **
!     **                                                              **
!     **   SUBROUTINES USED:                                          **
!     **    NOCULC                                                    **
!     **                                                              **
!     **  INPUT AND OUTPUT:                                           **
!     **    NFU: NUMBER OF DIFFERENT A(K) PER K-POINT                 **
!     **    NNOC: NUMBER OF ALLREADY FILLED TETRAHDRA                **
!     **    NEMAX: NUMBER OF LAST BAND                                **
!     **    NEOCC: NUMBER OF LAST OCCUPIED BAND                       **
!     **    NEMIN: NUMBER OF FIRST BAND                               **
!     **    NKPTK: NUMBER OF K-POINTS                                 **
!     **    NTT:  NUMBER OF DIFFERENT TETRAHEDRA                      **
!     **                                                              **
!      INCLUDE 'param.inc'
      use felder
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NITT = 505)
!      REAL*4  OPMAT
!ad
      COMMON /EME/ EEF,EMIN, EMAX, EFACTR, ESTEP, NFIRST, NFU, NLAST
      COMMON /NCOUNT/ NNOC
!      COMMON /BST/ EBS(NKPT,NUME), FC(NKPT,NUME)
      COMMON /EMICRO/ D(4), F(MG0,4), V , D1(4)
      COMMON /SWITCH/ISWITCH
      COMMON /TETRC/ RK(3,NEKPT),NEQ(NEKPT)
!      COMMON /OPME/ OPMAT(NKPT,INUME,MG), &
!                    EMINo,EMAXo,OML,OM1,MIMA(NKPT,2),NK,KRA
      COMMON /OPME/ EMINo,EMAXo,OML,OM1,NK,KRA
!      COMMON /SN/ DENSTY(INUMEden,MET,MG)
      DIMENSION ITTFL(nitt),KPN(4)
!ad
      EF=EEF
!ad
      NNOC=0
      READ(14,1234) Ndim,NTT,V1,NWRIT,NREC
      NTTW=min(int(NTT/NWRIT)+1,nrec)
!ad
!ad
      if(iswitch.eq.0.or.iswitch.eq.1) goto 100
      if(iswitch.eq.2.or.iswitch.eq.3) goto 200
      if(iswitch.eq.4.or.iswitch.eq.5) goto 400
      if(iswitch.eq.6.or.iswitch.eq.7) goto 600
!ad
  200 continue
!ad
!ad _____________________________ case DOS _____________________________
!ad
      DO 30 N=1,NTTW
      IF (N.EQ.NTTW) THEN
        KMAX=NTT-(N-1)*NWRIT
      ELSE
        KMAX=NWRIT
      ENDIF
!cad
      if(5*kmax.gt.NITT) then
      write(*,*) 'increase parameter NITT in arbdos!'
      write(*,*) 'the value should be at least',5*kmax
      STOP 'NITT in arbdos too small'
      endif
!cad
      READ(14,1235) (ITTFL(K),K=1,KMAX*5)
        DO 30 K=1,KMAX
        V=DBLE(ITTFL((K-1)*5+1))*V1
        IB=0
        DO 30 II=NEMIN,NEMAX
        IB=IB+1    

        DO 20 KP=1,4
        KPP=ITTFL((K-1)*5+KP+1)  
!ad
!.......no calculations if one of the energies isn't defined..........
!ad
        IF ((EBS(KPP,II).EQ.0)) GO TO 30
        D(KP)=EBS(KPP,II)  
        F(1,KP)=1  
   20 CONTINUE               

!ad
!ad............................BAND ANALYSIS............................
!ad
      if (iswitch.eq.2) then
        ibb=ib
      else
        ibb=1
      endif
!ad
      CALL NOCULC (IBB)  
!ad
   30 CONTINUE    
     
      goto 900
!ad
!ad  ___________________________ case DOS end __________________________
!ad
   
  100 continue
     
!ad
!ad ____________________________ case JDOS _____________________________
!ad
      DO 31 N=1,NTTW

      IF (N.EQ.NTTW) THEN
        KMAX=NTT-(N-1)*NWRIT
      ELSE
        KMAX=NWRIT
      ENDIF

      READ(14,1235) (ITTFL(K),K=1,KMAX*5)

      DO 31 K=1,KMAX
       V=DBLE(ITTFL((K-1)*5+1))*V1
       IB=0
!fb       DO 31 II=NEMIN,NEMAX-1
       DO 31 II=NEMIN,NEOCC
        DO 31 JJ=II+1,NEMAX  
         IB=IB+1    
         DO  KP=1,4
          KPP=ITTFL((K-1)*5+KP+1) 
          KPN(KP)=KPP
!
!....... no calculations if one of the energies isn't defined ..........
!
          IF ((EBS(KPP,II).EQ.0).OR.(EBS(KPP,JJ).EQ.0)) GO TO 31
          D(KP)=EBS(KPP,JJ)
          D1(KP)=EBS(KPP,II)           
          F(1,KP) = 1
         END DO
!ad
!ad............................BAND ANALYSIS............................
!ad
      if (iswitch.eq.0) then
        ibb=ib
      else
        ibb=1
      endif
!ad
!fb  
      socc = 0.0d0
      do KP=1,4
        KPP=KPN(KP)
        socc = socc + (FC(KPP,II)*(1-FC(KPP,JJ)))
      enddo
      if (abs(socc).gt.1.0d-10) then
        if (abs(socc-4.0d0).lt.1.0d-10) then
            CALL NOCULC (IBB)
        else
            CALL OPT1 (IBB)
        endif
      endif
!fb
!ad
   31 CONTINUE

      goto 900

!ad
!ad  __________________________ case JDOS end __________________________
!ad
!ad
!ad  ______________________ case DIELECTRIC TENSOR _____________________
!ad
 
  400 continue

!..................... TET ..........................
!.....NTETR=number of actual tetraeder !?...........
      NTETR=0
      DO 32 N=1,NTTW
      IF (N.EQ.NTTW) THEN
        KMAX=NTT-(N-1)*NWRIT
      ELSE
        KMAX=NWRIT
      ENDIF
      READ(14,1235) (ITTFL(K),K=1,KMAX*5)
        DO 32 K=1,KMAX
!..................... TET ..........................
      NTETR=NTETR+1
      mot=ITTFL((K-1)*5+1)
        V=DBLE(mot)*V1
        DO KP=1,4
         KPPc=ITTFL((K-1)*5+KP+1) 
         KPN(KP)=KPPc
        END DO
!ad
        IB=0
!fb        DO 132 II=NEMIN,NEMAX-1
        DO 132 II=NEMIN,NEOCC
        DO 132 JJ=II+1,NEMAX  
        IB=IB+1    
        DO KP=1,4
         KPP=KPN(KP)
             eii=EBS(KPP,II)
             ejj=EBS(KPP,JJ)
!
!.......no calculations if one of the energies isn't defined..........
         IF ((eii.EQ.0).OR.(ejj.EQ.0))  GOTO 132 
!ad
!.......no calculation if one of the matrixelements isn't defined.....
         NKMI=MIMA(KPP,1) 
         NKMA=MIMA(KPP,2)  
         IF ((II.GT.NKMA).OR.(II.LT.NKMI).OR. &
             (JJ.GT.NKMA).OR.(JJ.LT.NKMI)) GOTO 132    
!
!......................................................................
         III=II-NKMI+1
         JJJ=JJ-NKMI+1
!ad
!ad.. linear index for matrix elements (including diagonal elements) ..
!ad
         INDEX=((2*(NKMA-NKMI+1)-III)*(III-1))/2+JJJ
!ad
!ad........... to exclude diagonal elements add line below ............
!ad                                                 -III
!ad
         D(KP)=EBS(KPP,JJ)
         D1(KP)=EBS(KPP,II) 
!ad
!ad ........ stretching of bands (also below fermieenergy) ............
!ad
!ad      if (FC(KPP,JJ).eq.0.0d0) then
!ad         eshift=1.0d0 
!ad         DEE=EF+(EBS(KPP,JJ)-EF)*eshift-EBS(KPP,II)
!ad      end if 
!ad.................... end of streching bands ........................
!ad
!ad............... matrix elements are already squared ................
!ad
         DO I=1,NFU
         F(I,KP)=OPMAT(KPP,INDEX,I)
!ad
!ad      f(i,kp)=f(i,kp)/(dee+1.d-8)
!ad
         END DO 
         END DO
!ad
!ad............................BAND ANALYSIS............................
!ad
      if (iswitch.eq.5) then
         iib=ib
      else
         iib=1
      end if
!fb  
      socc = 0.0d0
      do KP=1,4
        KPP=KPN(KP)
        socc = socc + (FC(KPP,II)*(1-FC(KPP,JJ)))
      enddo
      if (abs(socc).gt.1.0d-10) then
        if (abs(socc-4.0d0).lt.1.0d-10) then
            CALL NOCULC (IIB)
        else
            CALL OPT1 (IIB)
        endif
      endif
!fb
  132 CONTINUE 
   32 CONTINUE
      goto 900

!ad
!ad  ____________________ case DIELECTRIC TENSOR end ___________________
!ad
     
  600 continue
  
!ad
!ad  ____________________ case INTRABAND CONTRIBUTIONS _________________
!ad
!.....NTETR=number of actual tetraeder !?...........
      NTETR=0
      DO 33 N=1,NTTW
      IF (N.EQ.NTTW) THEN
        KMAX=NTT-(N-1)*NWRIT
      ELSE
        KMAX=NWRIT
      ENDIF
      READ(14,1235) (ITTFL(K),K=1,KMAX*5)
        DO 33 K=1,KMAX
!ad
      NTETR=NTETR+1
      mot=ITTFL((K-1)*5+1)
        V=DBLE(mot)*V1
        DO KP=1,4
         KPPc=ITTFL((K-1)*5+KP+1) 
         KPN(KP)=KPPc
!ad
        END DO
!ad
        IB=0
!ad
        DO 133 II=NEMIN,NEMAX
        IB=IB+1    
        DO KP=1,4
         KPP=KPN(KP)
         eii=EBS(KPP,II)
!ad
!.......no calculations if one of the energies isn't defined..........
         IF (eii.EQ.0)  GOTO 133 
!ad
!.......no calculation if one of the matrixelements isn't defined.....
         NKMI=MIMA(KPP,1) 
         NKMA=MIMA(KPP,2)  
         IF ((II.GT.NKMA).OR.(II.LT.NKMI)) GOTO 133    
!ad ...................................................................
!ad
         DEE=EBS(KPP,II) 
!ad
         DO I=1,NFU              
           F(I,KP)= OPMAT(KPP,II,I)     
         END DO 
         D(KP)=DEE
        END DO
!ad
!ad............................BAND ANALYSIS............................
!ad
      if (iswitch.eq.7) then
           iib=ib
      else
           iib=1
      endif
!ad
      CALL NOCULC (IIB)
  133 CONTINUE 
   33 CONTINUE
!ad
!ad  __________________ case INTRABAND CONTRIBUTIONS end _______________
!ad
  
  900 continue
!
      RETURN
 1234 format(2i10,e20.12,2i10)
 1235 format(6i10)    
 1771 format(i10,i3,3x,9e12.4)
      END
