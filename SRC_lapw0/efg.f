      SUBROUTINE EFG(V20,V21,V22,V21M,V22M,V20SRF,V21SRF,V22SRF, &
       V21MSR,V22MSR,VZZ,qmat,tens,eivec,winkel,eta,RNOT,jatom)
      IMPLICIT NONE
      INTEGER CNT (3),EIGEN3,BASIS,RES,LOW,HIGH,jatom,JT,jmax,jmin
      INTEGER  :: jmiddl,i,j
      REAL*8   :: QMAT(3,3),tens(3),WERTI(3),EIVEC(3,3),SKAL(3)
      REAL*8   :: V20,V22,V22M,V21,V21M,V20SRF,V22SRF,V22MSR, &
        V21SRF,V21MSR,RNOT,eta,mat(3,3)
      REAL*8   :: vxx,vyy,vzz
      REAL*8   :: PI, EFACT, FAC,anorm,winkel,absmin,absmax

      DATA BASIS/2/
!        ELECTRIC FIELDGRADIENT: TRACELESS TENSOR    
!                                                                       
!              VXX : = ABSMIN( TENS(1),TENS(2),TENS(3) )                
!              VZZ : = ABSMAX( TENS(1),TENS(2),TENS(3) )                
!                                                                       
!              EFG = VZZ                                                
!              ETA = (VXX-VYY)/VZZ                                      
!                                                                       
!        FAC *   0.5     BECAUSE OF CONVERSION RYDBERG TO HARTREE       
!        FAC * 324.13909 BECAUSE OF CONVERSION A.U.    TO ESU/CM**3     
!        FAC * .02997925 EFG IN  10**17 V / (m*m)
!  
            PI=4.D0*DATAN(1.D0)                                             
            EFACT=0.02997925d0        
            FAC=0.5d0/RNOT**2*SQRT(15.d0/4.d0/PI)*324.13909d0*efact            
            V20   =V20*FAC                                              
           V20SRF=V20SRF*FAC                                           
            V22   =V22*FAC                                              
            V22SRF=V22SRF*FAC                                           
            V22M   =V22M*FAC                                              
            V22MSR=V22MSR*FAC                                           
            V21   =V21*FAC                                              
            V21SRF=V21SRF*FAC                                           
            V21M   =V21M*FAC                                              
            V21MSR=V21MSR*FAC                                           
!                                                                       
            qmat(1,1)=V22-V20/SQRT(3.d0)                                   
            qmat(2,2)=-V22-V20/SQRT(3.d0)                                  
            qmat(3,3)=2.d0*V20/SQRT(3.d0) 
            qmat(1,2)=V22M 
            qmat(2,1)=V22M
            qmat(1,3)=V21
            qmat(3,1)=V21
            qmat(2,3)=V21M
            qmat(3,2)=V21M
!
            do i=1,3
               do j=1,3
                  mat(i,j)=qmat(i,j)
               enddo
            enddo
!
      RES =EIGEN3(BASIS,3,3,MAT,SKAL,EIVEC,tens,WERTI,CNT,LOW,HIGH)
!
      IF (RES .GT. 0) THEN
      RES = RES-400
!     AUSGABE EINER ABBRUCHMELDUNG GEMAESS RES (0/401/402/403)
      GOTO (900,40,910),RES
  40  write(*,*)  'ABBRUCH: DIE EFG-MATRIX IST DIE NULLMATRIX !'
      ELSE
!
!     AUSGABE DER LOESUNG
            anorm=sqrt(eivec(1,1)**2+eivec(1,2)**2+eivec(1,3)**2)
            winkel=acos(eivec(1,1)/anorm) * 180.d0 / acos(-1.d0)
!                   
                                                    
!...........DEFINE VXX,VYY,VZZ ACCORDING ABS( TENS)                     
            ABSMAX=-1.                                                  
            ABSMIN=1000000.                                             
            DO 62 JT=1,3                                                
               IF( ABS(TENS(JT)).GT.ABSMAX ) THEN                       
                   JMAX=JT                                              
                   ABSMAX=ABS(TENS(JT))                                 
               ENDIF                                                    
               IF( ABS(TENS(JT)).LT.ABSMIN ) THEN                       
                   JMIN=JT                                              
                   ABSMIN=ABS(TENS(JT))                                 
               ENDIF                                                    
 62         CONTINUE                                                    
            JMIDDL=6-JMAX-JMIN                                          
            VXX=TENS(JMIN)                                              
            VYY=TENS(JMIDDL)                                            
            VZZ=TENS(JMAX)                                              
!                                                                       
!...........CALC EFG AND ASYMMETRIC PARAMETER ETA                          
            IF(ABS(VZZ).GT.0.000001) ETA=(VXX-VYY)/VZZ 
!                                                                       
!...........WRITE OUT EFG RESULTS
!
            WRITE(6,1340)   jatom,VZZ                                         
            WRITE(6,1341)   V20,V20SRF                                  
            WRITE(6,1342)   V22,V22SRF                                  
            WRITE(6,1343)   V22M,V22MSR                                  
            WRITE(6,1344)   V21,V21SRF                                  
            WRITE(6,1345)   V21M,V21MSR                                 
            WRITE(6,1350)   qmat(1,1),qmat(1,2),qmat(1,3),TENS(1),0.,0.
            WRITE(6,1350)   qmat(2,1),qmat(2,2),qmat(2,3),0.,TENS(2),0.
            WRITE(6,1350)   qmat(3,1),qmat(3,2),qmat(3,3),0.,0.,TENS(3)
            WRITE(6,1355)   ((eivec(i,j),j=1,3),i=1,3),jatom,winkel
            WRITE(6,1360)   jatom,ETA                                         
      ENDIF                                                       
      return
!
!        Error messages
!
  900 CALL OUTERR('EFG','order of matrix lt. 1')
      STOP 'EFG - Error'
  910 CALL OUTERR('EFG','max order of QR-alghorithm exceeded')
      STOP 'EFG - Error'
!
!
!
 1340 FORMAT(':EFG',i3.3,':',23X,'EFG         =',F12.5, &
             '   *10**21  V / m**2')         
 1341 FORMAT(31X,'V20  TOT/SRF=',2F12.5)                                 
 1342 FORMAT(31X,'V22  TOT/SRF=',2F12.5)                               
 1343 FORMAT(31X,'V22M TOT/SRF=',2F12.5)                               
 1344 FORMAT(31X,'V21  TOT/SRF=',2F12.5)                               
 1345 FORMAT(31X,'V21M TOT/SRF=',2F12.5,/)                               
 1350 FORMAT(10X,3F11.5,4X,3F11.5)                                            
 1355 FORMAT(/,9X,'MAIN DIRECTIONS OF THE EFG ',3F8.4,/,36X,3F8.4,/ &
      ,36X,3F8.4,/, &
      ':ANG',i3.3,':',1X,'ANGLE WITH OLD X-AXIS =    ',f8.1)
 1360 FORMAT(/,':ETA',i3.3,':',24X,'ASYMM. ETA =',F12.5,//) 
      end
      INTEGER FUNCTION EIGEN3 &
         (BASIS,ND,N,MAT,SKAL,EIVEC,WERTR,WERTI,CNT,LOW,HIGH)
!
!*****************************************************************
!                                                                *
!     ZWECK DES PROGRAMMS:                                       *
!     ====================                                       *
!     DIESES FUNKTIONSUNTERPROGRAMM VOM TYP INTEGER BERECHNET    *
!     ALLE EIGENWERTE UND EIGENVEKTOREN EINER BELIEBIGEN         *
!     REELLEN MATRIX MAT.                                        *
!     DIE EIGENWERTE WERDEN IN DEN FELDERN WERTR(1:N)            *
!     (REALTEIL) UND WERTI(1:N) (IMAGINAERTEIL) UND DIE          *
!     EIGENVEKTOREN IM FELD EIVEC(1:N,1:N) ABGESPEICHERT.        *
!                                                                *
!     EINGABEPARAMETER:                                          *
!     =================                                          *
!     N:           ORDNUNG DER MATRIX MAT                        *
!     ND:          FUEHRENDE DIMENSION DER FELDER MAT UND        *
!                  EIVEC, WIE IM RUFENDEN PROGRAMM VEREINBART    *
!     MAT:         EIN (N,N)-FELD VOM TYP DOUBLE PRECISION,      *
!                  DAS DIE MATRIX ENTHAELT, DEREN EIGENWERTE     *
!                  UND EIGENVEKTOREN BERECHNET WERDEN SOLLEN     *
!     BASIS:       DIE BASIS DER GLEITKOMMADARSTELLUNG DES       *
!                  VERWENDETEN RECHNERS (MEISTENS 2 ODER 16)     *
!                                                                *
!     AUSGABEPARAMETER:                                          *
!     =================                                          *
!     MAT:         DER OBERE TEIL DIESES (N,N)-FELDES ENTHAELT   *
!                  DIE EIGENVEKTOREN DER QUASI-DREIECKSMATRIX,   *
!                  DIE VOM QR-VERFAHREN ERZEUGT WIRD.            *
!     WERTR,WERTI: ZWEI (N,1)-FELDER VOM TYP DOUBLE              *
!                  PRECISION, DIE REALTEIL UND IMAGINAER-        *
!                  TEIL DER EIGENWERTE AUFNEHMEN                 *
!     EIVEC:       EIN (N,N)-FELD VOM TYP DOUBLE PRECISION,      *
!                  DAS DIE NORMALISIERTEN EIGENVEKTOREN DER      *
!                  URSPRUENGLICHEN VOLLEN MATRIX AUFNIMMT.       *
!                  FALLS DER I-TE EIGENWERT REELL IST, DANN      *
!                  IST DIE I-TE SPALTE VON EIVEC DER DAZUGE-     *
!                  HOERIGE REELLE EIGENVEKTOR. FALLS DIE         *
!                  EIGENWERTE I UND I+1 EIN KOMPLEXES            *
!                  PAAR BILDEN, GEBEN I-TE UND (I+1)-TE          *
!                  SPALTE REAL- UND IMAGINAERTEIL DESJENIGEN     *
!                  EIGENVEKTORS AN, DER ZU DEM EIGENWERT MIT     *
!                  POSITIVEM IMAGINAERTEIL GEHOERT.              *
!     CNT:         EIN (N,1)-FELD VOM TYP INTEGER, DAS DIE       *
!                  ZAHL DER ITERATIONSSCHRITTE FUER JEDEN        *
!                  EIGENWERT AUFNIMMT. FALLS ZWEI EIGENWERTE     *
!                  ALS PAAR GLEICHZEITIG GEFUNDEN WERDEN, WIRD   *
!                  DIE ZAHL DER ITERATIONSSCHRITTE MIT EINEM     *
!                  POSITIVEN VORZEICHEN FUER DEN ERSTEN UND      *
!                  EINEM NEGATIVEN VORZEICHEN FUER DEN           *
!                  ZWEITEN EIGENWERT EINGETRAGEN.                *
!     SKAL:        EIN (N,1)-FELD VOM TYP DOUBLE PRECISION, DAS  *
!                  DIE INFORMATION UEBER DIE DURCHGEFUEHRTEN     *
!                  VERTAUSCHUNGEN UND DIE SKALIERUNGSFAKTOREN    *
!                  ENTHAELT.                                     *
!                                                                *
!     RUECKGABEWERTE:                                            *
!     ===============                                            *
!     0:      KEIN FEHLER                                        *
!     401:    DIE ORDNUNG N DER MATRIX MAT IST KLEINER ALS 1.    *
!     402:    MAT IST DIE NULLMATRIX.                            *
!     403:    DIE MAXIMALE SCHRITTZAHL VON FUER DAS              *
!             QR-VERFAHREN IST UEBER- SCHRITTEN, OHNE DASS       *
!             ALLE EIGENWERTE BERECHNET WERDEN KONNTEN.          *
!                                                                *
!     LOKALE GROESSEN:                                           *
!     ================                                           *
!     ONE,TWO,HALF: GLEITKOMMAKONSTANTEN                         *
!     EPS:          MASCHINENGENAUIGKEIT                         *
!     TEMP:         HILFSVARIABLE                                *
!                                                                *
!----------------------------------------------------------------*
!                                                                *
!  BENOETIGTE UNTERPROGRAMME: BALAN, ELMHES, ELMTRA, HQR2,       *
!                             BALBAK, NORMAL, SWAP, COMDIV,      *
!                             COMABS                             *
!                                                                *
!                                                                *
!  QUELLEN : 1. MARTIN, R. S. UND WILKINSON, J. H.,              *
!               SIEHE [MART70].                                  *
!            2. PARLETT, B. N. UND REINSCH, C., SIEHE [PARL69].  *
!            3. PETERS, G. UND WILKINSON, J. H., SIEHE [PETE70]. *
!                                                                *
!*****************************************************************
!                                                                *
!  AUTOR     : JUERGEN DIETEL                                    *
!  DATUM     : 10.04.1987                                        *
!  QUELLCODE : FORTRAN 77                                        *
!                                                                *
!*****************************************************************
!
      INTEGER BASIS,ND,N,CNT(N),LOW,HIGH
      DOUBLE PRECISION MAT(ND,N),SKAL(N),EIVEC(ND,N),WERTR(N), &
                       WERTI(N)
      DOUBLE PRECISION ONE,TWO,HALF
      PARAMETER (ONE = 1.0D0,TWO = 2.0D0,HALF = 0.5D0)
      INTEGER RES,BALAN,ELMHES,ELMTRA,HQR2,BALBAK,NORMAL
      DOUBLE PRECISION EPS,TEMP
!
!     BERECHNUNG DER MASCHINENGENAUIGKEIT EPS (D. H. DER KLEINSTEN
!     POSITIVEN MASCHINENZAHL, FUER DIE AUF DEM RECHNER GILT: 
!     1 + EPS > 1): 
!
      TEMP = TWO
      EPS = ONE
   10 IF (ONE .LT. TEMP) THEN
          EPS = EPS * HALF
          TEMP = ONE + EPS
          GOTO 10
          ENDIF
      EPS = TWO * EPS
      RES = BALAN(ND,N,MAT,SKAL,LOW,HIGH,BASIS)
         IF (RES .NE. 0) THEN
            EIGEN3 = RES + 100
            RETURN
         ENDIF
      RES = ELMHES(ND,N,LOW,HIGH,MAT,CNT)
         IF (RES .NE. 0) THEN
            EIGEN3= RES + 200
            RETURN
         ENDIF
      RES = ELMTRA(ND,N,LOW,HIGH,MAT,CNT,EIVEC)
         IF (RES .NE. 0) THEN
            EIGEN3= RES + 300
            RETURN
         ENDIF
      RES = HQR2(ND,N,LOW,HIGH,MAT,WERTR,WERTI,EIVEC,CNT,EPS)
         IF (RES .NE. 0) THEN
            EIGEN3= RES + 400
            RETURN
         ENDIF
      RES = BALBAK(ND,N,LOW,HIGH,SKAL,EIVEC)
         IF (RES .NE. 0) THEN
            EIGEN3= RES + 500
            RETURN
         ENDIF
      RES = NORMAL(ND,N,EIVEC,WERTI)
         IF (RES .NE. 0) THEN
            EIGEN3= RES + 600
            RETURN
         ENDIF
      EIGEN3= 0
      RETURN
      END
!
!
      INTEGER FUNCTION BALAN (ND,N,MAT,SKAL,LOW,HIGH,BASIS)
!
!*****************************************************************
!                                                                *
!     ZWECK DES PROGRAMMS:                                       *
!     ====================                                       *
!     DIE PROZEDUR BALAN BALANCIERT EINE GEGEBENE REELLE MATRIX  *
!     IN DER 1-NORM AUS.                                         *
!                                                                *
!     EINGABEPARAMETER:                                          *
!     =================                                          *
!     ND:       FUEHRENDE DIMENSION DER MATRIZEN, WIE SIE IM     *
!               HAUPTPROGRAMM VEREINBART WURDEN                  *
!     N:        DIE ORDNUNG DER GEGEBENEN MATRIX                 *
!     MAT:      EIN (1:N,1:N)-FELD, DAS DIE KOMPONENTEN DER      *
!               GEGEBENEN MATRIX ENTHAELT                        *
!     BASIS:    DIE BASIS DER GLEITKOMMADARSTELLUNG AUF          *
!               DER MASCHINE                                     *
!                                                                *
!     AUSGABEPARAMETER:                                          *
!     =================                                          *
!     MAT:      DIE AUSBALANCIERTE MATRIX                        *
!     LOW,HIGH: ZWEI INTEGERZAHLEN, FUER DIE MAT(I,J)            *
!               GLEICH NULL IST, FALLS GILT:                     *
!               1. I>J UND                                       *
!               2. J=1,...LOW-1 ODER I=HIGH+1,...N               *
!     SKAL:     EIN (1:N)-FELD, DAS DIE INFORMATIONEN UEBER      *
!               DIE DURCHGEFUEHRTEN VERTAUSCHUNGEN UND DIE       *
!               SKALIERUNGSFAKTOREN ENTHAELT.                    *
!                                                                *
!     RUECKGABEWERT:                                             *
!     ==============                                             *
!     0:             KEIN FEHLER                                 *
!                                                                *
!     LOKALE GROESSEN:                                           *
!     ================                                           *
!     ZERO,ONE,PT95: GLEITKOMMAKONSTANTEN                        *
!     I,J,K,L:       ZAEHLVARIABLEN                              *
!     B2:            QUADRAT VON BASIS                           *
!     R,C,F,G,S:     HILFSVARIABLEN ZUR BERECHNUNG VON           *
!                    ZEILENNORMEN, IHREN KEHRWERTEN UND          *
!                    AEHNLICHEM                                  *
!                                                                *
!----------------------------------------------------------------*
!                                                                *
!  BENOETIGTE UNTERPROGRAMME: SWAP                               *
!                                                                *
!                                                                *
!  QUELLEN : PARLETT, B. N. UND REINSCH, C., SIEHE [PARL69].     *
!                                                                *
!*****************************************************************
!                                                                *
!  AUTOR     : JUERGEN DIETEL                                    *
!  DATUM     : 10.04.1987                                        *
!  QUELLCODE : FORTRAN 77                                        *
!                                                                *
!*****************************************************************
!
      INTEGER ND,N,LOW,HIGH,BASIS
      DOUBLE PRECISION SKAL(N),MAT(ND,N)
      DOUBLE PRECISION ZERO,ONE,PT95
      PARAMETER (ZERO = 0.0D0,ONE = 1.0D0,PT95 = 0.95D0)
      INTEGER I,J,K,L,B2
      DOUBLE PRECISION R,C,F,G,S
!
!     DIE NORM VON MAT(1:N,1:N) REDUZIEREN DURCH EXAKTE AEHNLICH-
!     KEITSTRANSFORMATIONEN, DIE IN SKAL(1:N) ABGESPEICHERT WERDEN
!
      B2 = BASIS*BASIS
      L = 1
      K = N
!
!     NACH ZEILEN MIT EINEM ISOLIERTEN EIGENWERT SUCHEN UND SIE
!     NACH UNTEN SCHIEBEN
!
   10 DO 50 J=K,1,-1
         R = ZERO
         DO 20 I=1,K
   20       IF (I .NE. J) R = R+ABS(MAT(J,I))
         IF (R .EQ. ZERO) THEN
            SKAL(K) = J
            IF (J .NE. K) THEN
               DO 30 I=1,K
   30             CALL SWAP(MAT(I,J),MAT(I,K))
               DO 40 I=L,N
   40             CALL SWAP(MAT(J,I),MAT(K,I))
            ENDIF
            K = K-1
            GOTO 10
         ENDIF
   50    CONTINUE
!
!     NACH SPALTEN MIT EINEM ISOLIERTEN EIGENWERT SUCHEN UND SIE
!     NACH LINKS SCHIEBEN
!
   60 DO 100 J=L,K
         C = ZERO
         DO 70 I=L,K
   70       IF (I .NE. J) C = C+ABS(MAT(I,J))
         IF (C .EQ. ZERO) THEN
            SKAL(L) = J
            IF (J .NE. L) THEN
               DO 80 I=1,K
   80             CALL SWAP(MAT(I,J),MAT(I,L))
               DO 90 I=L,N
   90             CALL SWAP(MAT(J,I),MAT(L,I))
            ENDIF
            L = L+1
            GOTO 60
         ENDIF
  100    CONTINUE
!
!     NUN DIE TEILMATRIX IN DEN ZEILEN L BIS K AUSBALANCIEREN
!
      LOW = L
      HIGH = K
      DO 110 I=L,K
  110    SKAL(I) = ONE
  120 DO 180 I=L,K
         C = ZERO
         R = ZERO
         DO 130 J=L,K
            IF (J .NE. I) THEN
               C = C+ABS(MAT(J,I))
               R = R+ABS(MAT(I,J))
            ENDIF
  130       CONTINUE
         G = R/BASIS
         F = ONE
         S = C+R
  140    IF (C .LT. G) THEN
            F = F*BASIS
            C = C*B2
            GOTO 140
         ENDIF
         G = R*BASIS
  150    IF (C .GE. G) THEN
            F = F/BASIS
            C = C/B2
            GOTO 150
         ENDIF
         IF ((C+R)/F .LT. PT95*S) THEN
            G = ONE/F
            SKAL(I) = SKAL(I)*F
            DO 160 J=L,N
  160          MAT(I,J) = MAT(I,J)*G
            DO 170 J=1,K
  170          MAT(J,I) = MAT(J,I)*F
            GOTO 120
         ENDIF
  180    CONTINUE
      BALAN = 0
      END
!
!
      INTEGER FUNCTION BALBAK(ND,N,LOW,HIGH,SKAL,EIVEC)
!
!*****************************************************************
!                                                                *
!     ZWECK DES PROGRAMMS:                                       *
!     ====================                                       *
!     BALBAK FUEHRT EINE RUECKTRANSFORMATION ALLER RECHTSEIGEN-  *
!     VEKTOREN EINER AUSBALANCIERTEN MATRIX IN DIE EIGENVEKTOREN *
!     DER ORIGINALMATRIX DURCH, VON DER DIE BALANCIERTE MATRIX   *
!     DURCH AUFRUF DER PROZEDUR BALAN ABGELEITET WURDE.          *
!                                                                *
!     EINGABEPARAMETER:                                          *
!     =================                                          *
!     ND:       FUEHRENDE DIMENSION DER MATRIZEN, WIE SIE IM     *
!               HAUPTPROGRAMM VEREINBART WURDEN                  *
!     N:        DIE ORDNUNG DER EIGENVEKTOREN (ZAHL DER          *
!               KOMPONENTEN)                                     *
!     LOW,HIGH: ZWEI INTEGERZAHLEN, DIE VON DER PROZEDUR         *
!               BALAN STAMMEN                                    *
!     SKAL:     AUSGABEVEKTOR DER PROZEDUR BALAN                 *
!     EIVEC:    EIN (1:N,1:N)-FELD, VON DEM JEDE SPALTE EINEN    *
!               EIGENVEKTOR (ODER SEINEN REALTEIL ODER SEINEN    *
!               IMAGINAERTEIL) DER AUSBALANCIERTEN MATRIX        *
!               DARSTELLT                                        *
!                                                                *
!     AUSGABEPARAMETER:                                          *
!     =================                                          *
!     EIVEC:    DIE ENTSPRECHENDEN EIGENVEKTOREN (ODER REAL-     *
!               TEILE ODER IMAGINAERTEILE) DER                   *
!               URSPRUENGLICHEN MATRIX                           *
!                                                                *
!     RUECKGABEWERT:                                             *
!     ==============                                             *
!     0:        KEIN FEHLER                                      *
!                                                                *
!     LOKALE VARIABLEN:                                          *
!     =================                                          *
!     I,J,K:    HILFSVARIABLEN ZUR INDEXBILDUNG                  *
!     S:   :    SKALIERUNGSWERT                                  *
!                                                                *
!----------------------------------------------------------------*
!                                                                *
!  BENOETIGTE UNTERPROGRAMME: SWAP                               *
!                                                                *
!                                                                *
!  QUELLEN : PARLETT, B. N. UND REINSCH, C., SIEHE [PARL69].     *
!                                                                *
!*****************************************************************
!                                                                *
!  AUTOR     : JUERGEN DIETEL                                    *
!  DATUM     : 10.04.1987                                        *
!  QUELLCODE : FORTRAN 77                                        *
!                                                                *
!*****************************************************************
!
      INTEGER ND,N,LOW,HIGH
      DOUBLE PRECISION SKAL(N),EIVEC(ND,N)
      INTEGER I,J,K
      DOUBLE PRECISION S
!
      DO 20 I=LOW,HIGH
         S = SKAL(I)
!
!        LINKSEIGENVEKTOREN WERDEN ZURUECKTRANSFORMIERT, INDEM MAN
!        DIE VORIGE ANWEISUNG ERSETZT DURCH: 'S = 1.0D0/SKAL(I)'
!
         DO 10 J=1,N
   10       EIVEC(I,J) = EIVEC(I,J)*S
   20    CONTINUE
      DO 40 I=LOW-1,1,-1
         K=SKAL(I)
         DO 30 J=1,N
   30       CALL SWAP(EIVEC(I,J),EIVEC(K,J))
   40    CONTINUE
      DO 60 I=HIGH+1,N
         K=SKAL(I)
         DO 50 J=1,N
   50       CALL SWAP(EIVEC(I,J),EIVEC(K,J))
   60    CONTINUE
      BALBAK = 0
      END
!
!
      INTEGER FUNCTION ELMHES(ND,N,LOW,HIGH,MAT,PERM)
!
!*****************************************************************
!                                                                *
!     ZWECK DES PROGRAMMS:                                       *
!     ====================                                       *
!     GEGEBEN IST EINE UNSYMMETRISCHE MATRIX A(1:N,1:N). DANN    *
!     REDUZIERT DIESE PROZEDUR DIE TEILMATRIX DER ORDNUNG        *
!     HIGH-LOW+1, DIE BEIM ELEMENT A(LOW,LOW) BEGINNT UND BEIM   *
!     ELEMENT A(HIGH,HIGH) ENDET, AUF HESSENBERGFORM H DURCH     *
!     NICHTORTHOGONALE ELEMENTARTRANSFORMATIONEN. DIE TEILMATRIX *
!     WIRD MIT H UEBERSCHRIEBEN, WOBEI DIE EINZELHEITEN DER      *
!     TRANSFORMATIONEN IN DEM UEBRIGBLEIBENDEN DREIECK UNTERHALB *
!     VON H UND IN DEM FELD PERM ABGESPEICHERT WERDEN.           *
!                                                                *
!     EINGABEPARAMETER:                                          *
!     =================                                          *
!     ND:       FUEHRENDE DIMENSION DER MATRIZEN, WIE SIE IM     *
!               HAUPTPROGRAMM VEREINBART WURDEN                  *
!     N:        ORDNUNG DER VOLLEN MATRIX A                      *
!     LOW,HIGH: AUSGABEPARAMETER EINER PROZEDUR, DIE             *
!               A AUFBEREITET (SIEHE 'PARLETT, B. N., AND C.     *
!               REINSCH: BALANCING A MATRIX FOR CALCULATION      *
!               OF EIGENVALUES AND EIGENVECTORS. NUMERISCHE      *
!               MATHEMATIK 13 (1969), SEITEN 293 - 304,          *
!               PROZEDUR BALANCE'). FALLS A NICHT DERART         *
!               AUFBEREITET IST, SETZE LOW:=1, HIGH:=N.          *
!     MAT:      DIE (N,N)-MATRIX A, NORMALERWEISE IN             *
!               AUFBEREITETER FORM (SIEHE OBEN)                  *
!                                                                *
!     AUSGABEPARAMETER:                                          *
!     =================                                          *
!     MAT:      EIN (N,N)-FELD, DAS ZU EINEM TEIL AUS DER        *
!               ABGELEITETEN OBEREN HESSENBERGMATRIX BESTEHT;    *
!               DIE GROESSE N(I,R+1), DIE BEI DER REDUKTION      *
!               EINE ROLLE SPIELT, WIRD IM (I,R)-ELEMENT         *
!               ABGESPEICHERT.                                   *
!     PERM:     EIN INTEGERFELD, DAS DIE BEI DER REDUKTION       *
!               AUSGEFUEHRTEN ZEILEN- UND SPALTENVERTAU-         *
!               SCHUNGEN BESCHREIBT                              *
!                                                                *
!     RUECKGABEWERT:                                             *
!     ==============                                             *
!     0:        KEIN FEHLER                                      *
!                                                                *
!     LOKALE GROESSEN:                                           *
!     ================                                           *
!     ZERO,ONE: GLEITKOMMAKONSTANTEN                             *
!     I,J,M:    ZAEHLVARIABLEN                                   *
!     X,Y:      HILFSVARIABLEN ZUR AUFNAHME VON MATRIXELEMENTEN  *
!               UND ZWISCHENERGEBNISSEN                          *
!                                                                *
!----------------------------------------------------------------*
!                                                                *
!  BENOETIGTE UNTERPROGRAMME: SWAP                               *
!                                                                *
!                                                                *
!  QUELLEN : MARTIN, R. S. UND WILKINSON, J. H., SIEHE [MART70]. *
!                                                                *
!*****************************************************************
!                                                                *
!  AUTOR     : JUERGEN DIETEL                                    *
!  DATUM     : 10.04.1987                                        *
!  QUELLCODE : FORTRAN 77                                        *
!                                                                *
!*****************************************************************
!
      INTEGER N,LOW,HIGH,PERM(N)
      DOUBLE PRECISION MAT(ND,N)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO = 0.0D0,ONE = 1.0D0)
      INTEGER I,J,M
      DOUBLE PRECISION X,Y
!
      DO 70 M=LOW+1,HIGH-1
         I = M
         X = ZERO
         DO 10 J=M,HIGH
            IF (ABS(MAT(J,M-1)) .GT. ABS(X)) THEN
               X = MAT(J,M-1)
               I=J
            ENDIF
   10       CONTINUE
         PERM(M) = I
         IF (I .NE. M) THEN
!
!           ZEILEN UND SPALTEN VON MAT VERTAUSCHEN
!
            DO 20 J=M-1,N
   20          CALL SWAP (MAT(I,J),MAT(M,J))
            DO 30 J=1,HIGH
   30          CALL SWAP (MAT(J,I),MAT(J,M))
         ENDIF
         IF (X .NE. ZERO) THEN
            DO 60 I=M+1,HIGH
               Y = MAT(I,M-1)
               IF (Y .NE. ZERO) THEN
                  Y = Y/X
                  MAT(I,M-1) = Y
                  DO 40 J=M,N
   40                MAT(I,J) = MAT(I,J)-Y*MAT(M,J)
                  DO 50 J=1,HIGH
   50                MAT(J,M) = MAT(J,M)+Y*MAT(J,I)
               ENDIF
   60          CONTINUE
         ENDIF
   70    CONTINUE
      ELMHES = 0
      END
!
!
      INTEGER FUNCTION ELMTRA(ND,N,LOW,HIGH,MAT,PERM,H)
!
!*****************************************************************
!                                                                *
!     ZWECK DES PROGRAMMS:                                       *
!     ====================                                       *
!     DIEJENIGE MATRIX IN DEM FELD H(1:N,1:N) ABSPEICHERN, DIE   *
!     SICH AUS DEN INFORMATIONEN ERGIBT, DIE DIE PROZEDUR ELMHES *
!     HINTERLASSEN HAT IM UNTEREN DREIECK DER HESSENBERGMATRIX   *
!     H, UND ZWAR IM FELD MAT(1:N,1:N) UND IM INTEGERFELD        *
!     PERM(1:N)                                                  *
!                                                                *
!     EINGABEPARAMETER:                                          *
!     =================                                          *
!     ND:       FUEHRENDE DIMENSION DER MATRIZEN, WIE SIE IM     *
!               HAUPTPROGRAMM VEREINBART WURDEN                  *
!     N:        ORDNUNG DER HESSENBERGMATRIX H                   *
!     LOW,HIGH: INTEGERZAHLEN, DIE VON DER PROZEDUR BALAN        *
!               ERZEUGT WURDEN (FALLS SIE VERWANDT WURDE;        *
!               ANDERNFALLS SETZE LOW:=1, HIGH:=N.)              *
!     PERM:     EIN VON ELMHES ERZEUGTES (N,1)-INTEGERFELD       *
!     MAT:      EIN (N,N)-FELD, DAS VON ELMHES ERZEUGT WURDE     *
!               UND DIE HESSENBERGMATRIX H UND DIE               *
!               MULTIPLIKATOREN ENTHAELT, DIE BENUTZT WURDEN,    *
!               UM ES AUS DER ALLGEMEINEN MATRIX                 *
!               A ZU ERZEUGEN                                    *
!                                                                *
!     AUSGABEPARAMETER:                                          *
!     =================                                          *
!     H:        DASJENIGE (N,N)-FELD, DAS DIE AEHNLICHKEITS-     *
!               TRANSFORMATION VON A IN H DEFINIERT              *
!                                                                *
!     RUECKGABEWERT:                                             *
!     ==============                                             *
!     0:             KEIN FEHLER                                 *
!                                                                *
!     LOKALE GROESSEN:                                           *
!     ================                                           *
!     ZERO,ONE: GLEITKOMMAKONSTANTEN                             *
!     I,J,K:    INDEXVARIABLEN                                   *
!                                                                *
!----------------------------------------------------------------*
!                                                                *
!  BENOETIGTE UNTERPROGRAMME: KEINE                              *
!                                                                *
!                                                                *
!  QUELLEN : PETERS, G. UND WILKINSON, J. H., SIEHE [PETE70].    *
!                                                                *
!*****************************************************************
!                                                                *
!  AUTOR     : JUERGEN DIETEL                                    *
!  DATUM     : 10.04.1987                                        *
!  QUELLCODE : FORTRAN 77                                        *
!                                                                *
!*****************************************************************
!

      INTEGER ND,N,LOW,HIGH,PERM(N)
      DOUBLE PRECISION MAT(ND,N),H(ND,N)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO = 0.0D0,ONE = 1.0D0)
      INTEGER I,J,K
      DO 20 I=1,N
         DO 10 J=1,N
   10       H(I,J) = ZERO
         H(I,I) = ONE
   20    CONTINUE
!
      DO 50 I=HIGH-1,LOW+1,-1
         J=PERM(I)
         DO 30 K=I+1,HIGH
   30       H(K,I) = MAT(K,I-1)
         IF (I .NE. J) THEN
            DO 40 K=I,HIGH
               H(I,K)=H(J,K)
               H(J,K)=ZERO
   40          CONTINUE
            H(J,I) = ONE
         ENDIF
   50    CONTINUE
      ELMTRA = 0
      END
!
!
      INTEGER FUNCTION HQR2 (ND,N,LOW,HIGH,H,WERTR,WERTI, &
                             EIVEC,CNT,EPS)
!
!*****************************************************************
!                                                                *
!     ZWECK DES PROGRAMMS:                                       *
!     ====================                                       *
!     FINDET DIE EIGENWERTE UND EIGENVEKTOREN EINER REELLEN      *
!     MATRIX, DIE, AUF OBERE HESSENBERGFORM REDUZIERT, IM FELD   *
!     H(1:N,1:N) STEHT, WOBEI DAS PRODUKT DER BISHER DURCHGE-    *
!     FUEHRTEN TRANSFORMATIONEN IM FELD EIVEC(1:N,1:N) STEHT.    *
!     DIE REAL- UND DIE IMAGINAERTEILE DER EIGENWERTE WERDEN IN  *
!     DEN FELDERN WERTR(1:N),WERTI(1:N) UND DIE EIGENVEKTOREN IM *
!     FELD EIVEC(1:N,1:N) GEBILDET, WO NUR EIN KOMPLEXER VEKTOR, *
!     NAEMLICH DER ZU DER WURZEL MIT POSITIVEM IMAGINAERTEIL GE- *
!     HOERIGE, FUER JEDES KOMPLEXE PAAR VON EIGENWERTEN ERZEUGT  *
!     WIRD. LOW UND HIGH SIND ZWEI INTEGERZAHLEN, DIE BEIM AUS-  *
!     BALANCIEREN ENTSTEHEN, WO EIGENWERTE IN DEN POSITIONEN 1   *
!     BIS LOW-1 UND HIGH+1 BIS N ISOLIERT WERDEN.  FALLS KEINE   *
!     AUSBALANCIERUNG DURCHGEFUEHRT WURDE, SETZE LOW:=1,HIGH:=N. *
!     DAS UNTERPROGRAMM BRICHT MIT EINER FEHLERMELDUNG AB, FALLS *
!     IRGENDEIN EIGENWERT MEHR ALS MAXSTP ITERATIONSSCHRITTE     *
!     BENOETIGT.                                                 *
!                                                                *
!     EINGABEPARAMETER:                                          *
!     =================                                          *
!     N:        ORDNUNG DER HESSENBERGMATRIX H                   *
!     LOW,HIGH: VON BALAN ERZEUGTE INTEGERZAHLEN, FALLS          *
!               BALAN BENUTZT WURDE. ANSONSTEN SETZE             *
!               LOW:=1, HIGH:=N.                                 *
!     EPS:      DIE KLEINSTE ZAHL AUF DEM COMPUTER, FUER         *
!               DIE GILT: 1 + EPS > 1.                           *
!     H:        EIN (N,N)-FELD, DAS DIE MATRIX H IN IHREN        *
!               RELEVANTEN TEILEN ENTHAELT                       *
!     EIVEC:    EIN (N,N)-FELD, DAS DIE MATRIX ENTHAELT, DIE     *
!               DIE AEHNLICHKEITSTRANSFORMATION VON A IN H       *
!               DEFINIERT. (ES WIRD VON ELMTRA ERZEUGT.)         *
!               FALLS H DIE URSPRUENGLICHE MATRIX IST,           *
!               SETZE EIVEC := EINHEITSMATRIX.                   *
!                                                                *
!     AUSGABEPARAMETER:                                          *
!     =================                                          *
!     H:           DER OBERE TEIL DIESES (N,N)-FELDES ENT-       *
!                  HAELT DIE EIGENVEKTOREN DER QUASI-            *
!                  DREIECKSMATRIX, DIE VOM QR-VERFAHREN          *
!                  ERZEUGT WIRD.                                 *
!     WERTR,WERTI: ZWEI (N,1)-FELDER, DIE REAL- UND IMAGI-       *
!                  NAERTEIL DER EIGENWERTE AUFNEHMEN             *
!     CNT:         EIN (N,1)-INTEGERFELD, DAS DIE ZAHL DER       *
!                  ITERATIONSSCHRITTE FUER JEDEN EIGENWERT       *
!                  AUFNIMMT. FALLS ZWEI EIGENWERTE ALS PAAR      *
!                  GLEICHZEITIG GEFUNDEN WERDEN, DANN WIRD       *
!                  DIE ZAHL DER ITERATIONSSCHRITTE MIT EINEM     *
!                  POSITIVEN VORZEICHEN FUER DEN ERSTEN UND      *
!                  EINEM NEGATIVEN VORZEICHEN FUER DEN           *
!                  ZWEITEN EIGENWERT EINGETRAGEN.                *
!     EIVEC:       EIN (N,N)-FELD, DAS DIE NICHTNORMALISIER-     *
!                  TEN EIGENVEKTOREN DER URSPRUENGLICHEN         *
!                  VOLLEN MATRIX AUFNIMMT (FALLS H NICHT DIE     *
!                  AUSGANGSMATRIX WAR). FALLS DER I-TE EIGEN-    *
!                  WERT REELL IST, DANN IST DIE I-TE SPALTE      *
!                  VON EIVEC DER DAZUGEHOERIGE REELLE EIGEN-     *
!                  VEKTOR. FALLS DIE EIGENWERTE I UND I+1 EIN    *
!                  KOMPLEXES PAAR BILDEN, GEBEN I-TE UND         *
!                  (I+1)-TE SPALTE REAL- UND IMAGINAERTEIL       *
!                  DESJENIGEN EIGENVEKTORS AN, DER ZU DEM        *
!                  EIGENWERT MIT POSITIVEM IMAGINAERTEIL         *
!                  GEHOERT.                                      *
!                                                                *
!     RUECKGABEWERTE:                                            *
!     ===============                                            *
!     0:           KEIN FEHLER                                   *
!     1:           DIE PARAMETER N, LOW ODER HIGH HABEN          *
!                  UNERLAUBTE WERTE.                             *
!     2:           ALLE EIGENVEKTOREN SIND DER NULLVEKTOR.       *
!     3:           DIE MAXIMALE SCHRITTZAHL IST UEBERSCHRITTEN.  *
!                                                                *
!     LOKALE GROESSEN:                                           *
!     ================                                           *
!     ZERO,ONE,TWO,PT75,PT4375: WICHTIGE GLEITKOMMAKONSTANTEN    *
!     MAXSTP:                   KONSTANTE FUER DIE MAXIMALE      *
!                               SCHRITTZAHL                      *
!     I,J,K,L,M,N,NA,EN:        INDEXVARIABLEN                   *
!     ITER:                     SCHRITTZAEHLER                   *
!     P,Q,R,S,T,W,X,Y,Z,NORM,                                    *
!       RA,SA,VR,VI:            HILFSVARIABLEN FUER GLEITKOMMA-  *
!                               BERECHNUNGEN                     *
!                                                                *
!----------------------------------------------------------------*
!                                                                *
!  BENOETIGTE UNTERPROGRAMME: COMDIV                             *
!                                                                *
!                                                                *
!  QUELLEN : PETERS, G. UND WILKINSON, J. H., SIEHE [PETE70].    *
!                                                                *
!*****************************************************************
!                                                                *
!  AUTOR     : JUERGEN DIETEL                                    *
!  DATUM     : 10.04.1987                                        *
!  QUELLCODE : FORTRAN 77                                        *
!                                                                *
!*****************************************************************
!
      INTEGER ND,N,LOW,HIGH,CNT(N)
      DOUBLE PRECISION H(ND,N),EIVEC(ND,N),WERTR(N),WERTI(N),EPS
      DOUBLE PRECISION ZERO,ONE,TWO,PT75,PT4375
      INTEGER MAXSTP
      PARAMETER (ZERO = 0.0D0,ONE = 1.0D0,TWO = 2.0D0, &
                 PT75 = 0.75D0,PT4375 = 0.4375D0,MAXSTP = 100)
      INTEGER I,J,K,L,M,NA,EN,ITER
      DOUBLE PRECISION P,Q,R,S,T,W,X,Y,Z,NORM,RA,SA,VR,VI
!
!     FEHLER 1: DIE PARAMETER N, LOW ODER HIGH HABEN UNERLAUBTE
!     WERTE:
!
      IF (N .LT. 1 .OR. LOW .LT. 1 .OR. HIGH .GT. N) THEN
         HQR2 = 1
         RETURN
      ENDIF
!
!     VORBESETZUNG FUER DIE BEI DER AUSBALANCIERUNG GEFUNDENEN
!     ISOLIERTEN EIGENWERTE: 
!
      DO 10 I=1,N
         IF (I .LT. LOW .OR. I .GT. HIGH) THEN
            WERTR(I) = H(I,I)
            WERTI(I) = ZERO
            CNT(I) = 0
         ELSE
            CNT(I) = -1
         ENDIF
   10    CONTINUE
!
      EN = HIGH
      T = ZERO
   15    IF (EN .LT. LOW) GOTO 333
         ITER = 0
         NA = EN-1
!
!           NACH EINEM EINZELNEN KLEINEN SUBDIAGONALELEMENT
!           SUCHEN: 
!
   20       DO 30 L=EN,LOW+1,-1
   30          IF (ABS(H(L,L-1)) .LE. EPS* &
                   (ABS(H(L-1,L-1))+ABS(H(L,L)))) GOTO 40
   40       X = H(EN,EN)
            IF (L .EQ. EN) THEN
!
!              EINE WURZEL GEFUNDEN: 
!
               WERTR(EN) = X + T
               H(EN,EN) = WERTR(EN)
               WERTI(EN) = ZERO
               CNT(EN) = ITER
               EN = NA
               GOTO 15
            ENDIF
!
            Y = H(NA,NA)
            W = H(EN,NA)*H(NA,EN)
            IF (L .EQ. NA) THEN
!
!              ZWEI WURZELN GEFUNDEN: 
!
               P = (Y-X)/TWO
               Q = P*P+W
               Z = SQRT(ABS(Q))
               H(EN,EN) = X+T
               X = H(EN,EN)
               H(NA,NA) = Y+T
               CNT(EN) = -ITER
               CNT(NA) = ITER
               IF (Q .GE. ZERO) THEN
!
!                 EIN REELLES PAAR GEFUNDEN: 
!
                  IF (P .LT. ZERO) Z = -Z
                  Z = P+Z
                  WERTR(NA) = X+Z
                  WERTR(EN) = X-W/Z
                  WERTI(NA) = ZERO
                  WERTI(EN) = ZERO
                  X = H(EN,NA)
                  R = SQRT(X*X+Z*Z)
                  P = X/R
                  Q = Z/R
!
!                 ZEILENMODIFIKATION: 
!
                  DO 50 J=NA,N
                     Z = H(NA,J)
                     H(NA,J) = Q*Z+P*H(EN,J)
                     H(EN,J) = Q*H(EN,J)-P*Z
   50                CONTINUE
!
!                 SPALTENMODIFIKATION: 
!
                  DO 60 I=1,EN
                     Z = H(I,NA)
                     H(I,NA) = Q*Z+P*H(I,EN)
                     H(I,EN) = Q*H(I,EN)-P*Z
   60                CONTINUE
!
!                 AKKUMULATION: 
!
                  DO 70 I=LOW,HIGH
                     Z = EIVEC(I,NA)
                     EIVEC(I,NA) = Q*Z+P*EIVEC(I,EN)
                     EIVEC(I,EN) = Q*EIVEC(I,EN)-P*Z
   70                CONTINUE
               ELSE
!
!                 KOMPLEXES PAAR: 
!
                  WERTR(NA) = X+P
                  WERTR(EN) = WERTR(NA)
                  WERTI(NA) = Z
                  WERTI(EN) = -Z
               ENDIF
               EN = EN-2
               GOTO 15
            ENDIF
!
            IF (ITER .EQ. MAXSTP) THEN
!
!              FEHLER 3: MAXIMALE SCHRITTZAHL UEBERSCHRITTEN:
!
               CNT(EN) = MAXSTP+1
               HQR2 = 3
               RETURN
            ENDIF
            IF (MOD(ITER,10) .EQ. 0 .AND. ITER .NE. 0) THEN
!
!              EINEN UNGEWOEHNLICHEN SHIFT DURCHFUEHREN: 
!
               T = T+X
               DO 80 I=LOW,EN
   80             H(I,I) = H(I,I)-X
               S = ABS(H(EN,NA))+ABS(H(NA,EN-2))
               X = PT75*S
               Y = X
               W = -PT4375*S*S
            ENDIF
            ITER = ITER+1
!
!           NACH ZWEI AUFEINANDERFOLGENDEN KLEINEN
!           SUBDIAGONALELEMENTEN SUCHEN: 
!
            DO 90 M=EN-2,L,-1
               Z = H(M,M)
               R = X-Z
               S = Y-Z
               P = (R*S-W)/H(M+1,M)+H(M,M+1)
               Q = H(M+1,M+1)-Z-R-S
               R = H(M+2,M+1)
               S = ABS(P)+ABS(Q)+ABS(R)
               P = P/S
               Q = Q/S
               R = R/S
               IF (M .EQ. L) GOTO 100
               IF (ABS(H(M,M-1))*(ABS(Q)+ABS(R)) .LE. EPS*ABS(P)* &
                   (ABS(H(M-1,M-1))+ABS(Z)+ABS(H(M+1,M+1)))) &
                   GOTO 100
   90          CONTINUE
  100       DO 110 I=M+2,EN
  110          H(I,I-2) = ZERO
            DO 120 I=M+3,EN
  120          H(I,I-3) = ZERO
!
!           EIN DOPPELTER QR-SCHRITT, DER DIE ZEILEN L BIS EN UND
!           DIE SPALTEN M BIS EN DES GANZEN FELDES BETRIFFT: 
!
            DO 200 K=M,NA
               IF (K .NE. M) THEN
                  P = H(K,K-1)
                  Q = H(K+1,K-1)
                  IF (K .NE. NA) THEN
                     R = H(K+2,K-1)
                  ELSE
                     R = ZERO
                  ENDIF
                  X = ABS(P)+ABS(Q)+ABS(R)
                  IF (X .EQ. ZERO) GOTO 200
                  P = P/X
                  Q = Q/X
                  R = R/X
               ENDIF
               S = SQRT(P*P+Q*Q+R*R)
               IF (P .LT. ZERO) S = -S
               IF (K .NE. M) THEN
                  H(K,K-1) = -S*X
               ELSEIF (L .NE. M) THEN
                  H(K,K-1) = -H(K,K-1)
               ENDIF
               P = P+S
               X = P/S
               Y = Q/S
               Z = R/S
               Q = Q/P
               R = R/P
!
!              ZEILENMODIFIKATION: 
!
               DO 130 J=K,N
                  P = H(K,J)+Q*H(K+1,J)
                  IF (K .NE. NA) THEN
                     P = P+R*H(K+2,J)
                     H(K+2,J) = H(K+2,J)-P*Z
                  ENDIF
                  H(K+1,J) = H(K+1,J)-P*Y
                  H(K,J) = H(K,J)-P*X
  130             CONTINUE
               J = MIN(K+3,EN)
!
!              SPALTENMODIFIKATION: 
!
               DO 140 I=1,J
                  P = X*H(I,K)+Y*H(I,K+1)
                  IF (K .NE. NA) THEN
                     P = P+Z*H(I,K+2)
                     H(I,K+2) = H(I,K+2)-P*R
                  ENDIF
                  H(I,K+1) = H(I,K+1)-P*Q
                  H(I,K) = H(I,K)-P
  140             CONTINUE
!
!              TRANSFORMATIONEN AKKUMULIEREN: 
!
               DO 150 I=LOW,HIGH
                  P = X*EIVEC(I,K)+Y*EIVEC(I,K+1)
                  IF (K .NE. NA) THEN
                     P = P+Z*EIVEC(I,K+2)
                     EIVEC(I,K+2) = EIVEC(I,K+2)-P*R
                  ENDIF
                  EIVEC(I,K+1) = EIVEC(I,K+1)-P*Q
                  EIVEC(I,K) = EIVEC(I,K)-P
  150             CONTINUE
  200          CONTINUE
            GOTO 20
!
!
!     ALLE WURZELN GEFUNDEN, NUN WIRD RUECKTRANSFORMIERT: 
!
!     1-NORM VON H BESTIMMEN:
!
  333 NORM = ZERO
      K=1
      DO 201 I=1,N
         DO 101 J=K,N
  101       NORM = NORM+ABS(H(I,J))
  201    K = I
      IF (NORM .EQ. ZERO) THEN
!        FEHLER 2: 1-NORM VON H IST GLEICH 0: 
         HQR2 = 2
         RETURN
      ENDIF
!
!     RUECKTRANSFORMATION: 
!
      DO 207 EN=N,1,-1
         P = WERTR(EN)
         Q = WERTI(EN)
         NA = EN - 1
         IF (Q .EQ. ZERO) THEN
!
!           REELLER VEKTOR: 
!
            M = EN
            H(EN,EN) = ONE
            DO 63 I=NA,1,-1
               W = H(I,I)-P
               R = H(I,EN)
               DO 38 J=M,NA
   38             R = R+H(I,J)*H(J,EN)
               IF (WERTI(I) .LT. ZERO) THEN
                  Z = W
                  S = R
               ELSE
                  M = I
                  IF (WERTI(I) .EQ. ZERO) THEN
                     IF (W .NE. ZERO) THEN
                        H(I,EN) = -R/W
                     ELSE
                        H(I,EN) = -R/(EPS*NORM)
                     ENDIF
                  ELSE
!
!                    LOESE DAS GLEICHUNGSSYSTEM: 
!                    [ W   X ] [ H(I,EN)   ]   [ -R ]
!                    [       ] [           ] = [    ]
!                    [ Y   Z ] [ H(I+1,EN) ]   [ -S ]
!
                     X = H(I,I+1)
                     Y = H(I+1,I)
                     Q = (WERTR(I)-P)*(WERTR(I) - P)+ &
                         WERTI(I)*WERTI(I)
                     T = (X*S-Z*R)/Q
                     H(I,EN) = T
                     IF (ABS(X) .GT. ABS(Z)) THEN
                        H(I+1,EN) = (-R-W*T)/X
                     ELSE
                        H(I+1,EN) = (-S-Y*T)/Z
                     ENDIF
                  ENDIF
               ENDIF
   63          CONTINUE
         ELSEIF (Q .LT. ZERO) THEN
!
!           KOMPLEXER VEKTOR, DER ZU LAMBDA = P - I * Q GEHOERT: 
!
            M = NA
            IF (ABS(H(EN,NA)) .GT. ABS(H(NA,EN))) THEN
               H(NA,NA) = -(H(EN,EN)-P)/H(EN,NA)
               H(NA,EN) = -Q/H(EN,NA)
            ELSE
               CALL COMDIV(-H(NA,EN),ZERO,H(NA,NA)-P,Q, &
                         H(NA,NA),H(NA,EN))
            ENDIF
            H(EN,NA) = ONE
            H(EN,EN) = ZERO
            DO 190 I=NA-1,1,-1
               W = H(I,I)-P
               RA = H(I,EN)
               SA = ZERO
               DO 75 J=M,NA
                  RA = RA+H(I,J)*H(J,NA)
                  SA = SA+H(I,J)*H(J,EN)
   75             CONTINUE
               IF (WERTI(I) .LT. ZERO) THEN
                  Z = W
                  R = RA
                  S = SA
               ELSE
                  M = I
                  IF (WERTI(I) .EQ. ZERO) THEN
                      CALL COMDIV(-RA,-SA,W,Q,H(I,NA),H(I,EN))
                  ELSE
!
!                    LOESE DIE KOMPLEXEN GLEICHUNGEN:
!           [ W+Q*I   X   ] [H(I,NA)+H(I,EN)*I    ]   [-RA-SA*I]
!           [             ] [                     ] = [        ]
!           [   Y   Z+Q*I ] [H(I+1,NA)+H(I+1,EN)*I]   [-R-S*I  ]
!
                     X = H(I,I+1)
                     Y = H(I+1,I)
                     VR = (WERTR(I)-P)*(WERTR(I)-P)+ &
                          WERTI(I)*WERTI(I)-Q*Q
                     VI = TWO*Q*(WERTR(I)-P)
                     IF (VR .EQ. ZERO .AND. VI .EQ. ZERO) VR = &
                        EPS*NORM* &
                        (ABS(W)+ABS(Q)+ABS(X)+ABS(Y)+ABS(Z))
                     CALL COMDIV(X*R-Z*RA+Q*SA,X*S-Z*SA-Q*RA, &
                                 VR,VI,H(I,NA),H(I,EN))
                     IF (ABS(X) .GT. ABS(Z)+ABS(Q)) THEN
                        H(I+1,NA) = (-RA-W*H(I,NA)+Q*H(I,EN))/X
                        H(I+1,EN) = (-SA-W*H(I,EN)-Q*H(I,NA))/X
                     ELSE
                       CALL COMDIV(-R-Y*H(I,NA),-S-Y*H(I,EN),Z,Q, &
                                 H(I+1,NA),H(I+1,EN))
                     ENDIF
                  ENDIF
               ENDIF
  190          CONTINUE
         ENDIF
  207    CONTINUE
!
!     ZU DEN ISOLIERTEN WURZELN GEHOERIGE VEKTOREN: 
!
      DO 230 I=1,N
         IF (I .LT. LOW .OR. I .GT. HIGH) THEN
            DO 220 J=I+1,N
  220          EIVEC(I,J) = H(I,J)
         ENDIF
  230    CONTINUE
!
!     MIT DER TRANSFORMATIONSMATRIX MULTIPLIZIEREN, UM DIE
!     VEKTOREN DER URSPRUENGLICHEN VOLLEN MATRIX ZU ERHALTEN: 
!
      DO 300 J=N,LOW,-1
         IF (J .LE. HIGH) THEN
            M = J
         ELSE
            M = HIGH
         ENDIF
         L = J-1
         IF (WERTI(J) .LT. ZERO) THEN
            DO 330 I=LOW,HIGH
               Y = ZERO
               Z = ZERO
               DO 320 K=LOW,M
                  Y = Y+EIVEC(I,K)*H(K,L)
                  Z = Z+EIVEC(I,K)*H(K,J)
  320             CONTINUE
               EIVEC(I,L) = Y
               EIVEC(I,J) = Z
  330          CONTINUE
         ELSE
            IF (WERTI(J) .EQ. ZERO) THEN
               DO 350 I=LOW,HIGH
                  Z = ZERO
                  DO 340 K=LOW,M
  340                Z =Z+EIVEC(I,K)*H(K,J)
  350             EIVEC(I,J) = Z
            ENDIF
         ENDIF
  300    CONTINUE
!
!     RUECKGABE VON 0: KEIN FEHLER: 
!
      HQR2 = 0
      END
!
!
      SUBROUTINE COMDIV (AR,AI,BR,BI,RESR,RESI)
!
!*****************************************************************
!                                                                *
!     ZWECK DES PROGRAMMS:                                       *
!     ====================                                       *
!     KOMPLEXE DIVISION: RESR+I*RESI := (AR+I*AI)/(BR+I*BI).     *
!     (DIESE PROZEDUR SOLLTE NICHT MIT                           *
!      BR=BI=0 AUFGERUFEN WERDEN.)                               *
!                                                                *
!     EINGABEPARAMETER:                                          *
!     =================                                          *
!     AR,AI:     REAL- UND IMAGINAERTEIL DES DIVIDENDEN          *
!     BR,BI:     REAL- UND IMAGINAERTEIL DES DIVISORS            *
!                                                                *
!     AUSGABEPARAMETER:                                          *
!     =================                                          *
!     RESR,RESI: REAL- UND IMAGINAERTEIL DES QUOTIENTEN          *
!                                                                *
!     LOKALE GROESSEN:                                           *
!     ================                                           *
!     ZERO:              GLEITKOMMAKONSTANTE 0                   *
!     TEMP1,TEMP2,TEMP3: HILFSVARIABLEN ZUR SPEICHERUNG VON      *
!                        ZWISCHENERGEBNISSEN                     *
!                                                                *
!----------------------------------------------------------------*
!                                                                *
!  BENOETIGTE UNTERPROGRAMME: KEINE                              *
!                                                                *
!                                                                *
!  QUELLEN : MARTIN, R. S. UND WILKINSON, J. H., SIEHE [MART70]. *
!                                                                *
!*****************************************************************
!                                                                *
!  AUTOR     : JUERGEN DIETEL                                    *
!  DATUM     : 10.04.1987                                        *
!  QUELLCODE : FORTRAN 77                                        *
!                                                                *
!*****************************************************************
!
      DOUBLE PRECISION AR,AI,BR,BI,RESR,RESI
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO = 0.0D0)
      DOUBLE PRECISION TEMP1,TEMP2,TEMP3
!
      IF (BR .EQ. ZERO .AND. BI .EQ. ZERO) THEN
         RESR = ZERO
         RESI = ZERO
         RETURN
      ENDIF
      IF (ABS(BR) .GT. ABS(BI)) THEN
         TEMP1 = BI/BR
         TEMP2 = TEMP1*BI+BR
         TEMP3 = (AR+TEMP1*AI)/TEMP2
         RESI = (AI-TEMP1*AR)/TEMP2
         RESR = TEMP3
      ELSE
         TEMP1 = BR/BI
         TEMP2 = TEMP1*BR+BI
         TEMP3 = (TEMP1*AR+AI)/TEMP2
         RESI = (TEMP1*AI-AR)/TEMP2
         RESR = TEMP3
      ENDIF
      END
!
!
      DOUBLE PRECISION FUNCTION COMABS(AR,AI)
!
!*****************************************************************
!                                                                *
!     ZWECK DES PROGRAMMS:                                       *
!     ====================                                       *
!     BERECHNUNG DES BETRAGES DER KOMPLEXEN ZAHL AR+I*AI:        *
!     COMABS:=SQRT(AR*AR+AI*AI)                                  *
!                                                                *
!     EINGABEPARAMETER:                                          *
!     =================                                          *
!     AR,AI: REAL- UND IMAGINAERTEIL DER KOMPLEXEN ZAHL, DEREN   *
!            BETRAG ZU BERECHNEN IST                             *
!                                                                *
!     AUSGABEPARAMETER:                                          *
!     =================                                          *
!     KEINE                                                      *
!                                                                *
!     RUECKGABEWERT:                                             *
!     ==============                                             *
!     BETRAG DES KOMPLEXEN PARAMETERS                            *
!                                                                *
!     LOKALE GROESSEN:                                           *
!     ================                                           *
!     ZERO,ONE:    KONSTANTEN                                    *
!     TEMP1,TEMP2: HILFSVARIABLEN ZUR SPEICHERUNG VON            *
!                  ZWISCHENERGEBNISSEN                           *
!                                                                *
!----------------------------------------------------------------*
!                                                                *
!  BENOETIGTE UNTERPROGRAMME: SWAP                               *
!                                                                *
!                                                                *
!  QUELLEN : MARTIN, R. S. UND WILKINSON, J. H., SIEHE [MART70]. *
!                                                                *
!*****************************************************************
!                                                                *
!  AUTOR     : JUERGEN DIETEL                                    *
!  DATUM     : 10.04.1987                                        *
!  QUELLCODE : FORTRAN 77                                        *
!                                                                *
!*****************************************************************
!
      DOUBLE PRECISION AR,AI
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO = 0.0D0,ONE = 1.0D0)
      DOUBLE PRECISION TEMP1,TEMP2
      TEMP1 = ABS(AR)
      TEMP2 = ABS(AI)
      IF (AR .EQ. ZERO .OR. AI .EQ. ZERO) THEN
         COMABS = ZERO
         RETURN
      ENDIF
      IF (TEMP2 .GT. TEMP1) CALL SWAP (TEMP1,TEMP2)
      IF (TEMP2 .EQ. ZERO) THEN
         COMABS = TEMP1
      ELSE
         COMABS = TEMP1*SQRT(ONE+(TEMP2/TEMP1)**2)
      ENDIF
      END
!
!
      INTEGER FUNCTION NORMAL (ND,N,V,WI)
!
!*****************************************************************
!                                                                *
!     ZWECK DES PROGRAMMS:                                       *
!     ====================                                       *
!     NORMAL NORMALISIERT DIE EIGENVEKTOREN IN DER MAXIMUMNORM.  *
!                                                                *
!     EINGABEPARAMETER:                                          *
!     =================                                          *
!     ND:     FUEHRENDE DIMENSION DER MATRIX V, WIE SIE IM       *
!             HAUPTPROGRAMM VEREINBART WURDE                     *
!     N:      DIE ORDNUNG DER MATRIX V                           *
!     V:      EIN (N,N)-FELD VOM TYP DOUBLE PRECISION, DAS       *
!             SPALTENWEISE DIE EIGENVEKTOREN ENTHAELT            *
!             (SIEHE EIGEN)                                      *
!     WI:     EIN FELD MIT N KOMPONENTEN VOM TYP                 *
!             DOUBLE PRECISION, DESSEN KOMPONENTEN DIE           *
!             IMAGINAERTEILE DER EIGENWERTE SIND                 *
!                                                                *
!     AUSGABEPARAMETER:                                          *
!     =================                                          *
!     V:      MATRIX DER NORMALISIERTEN EIGENVEKTOREN            *
!                                                                *
!     LOKALE GROESSEN:                                           *
!     ================                                           *
!     ZERO,ONE: GLEITKOMMAKONSTANTEN 0 UND 1                     *
!     I,J:      INDEXVARIABLEN                                   *
!     MAXI:     HILFSVARIABLE ZUR BERECHNUNG DER REELLEN         *
!               VEKTORNORM                                       *
!     TR,TI:    HILFSVARIABLEN ZUR BERECHNUNG DER KOMPLEXEN      *
!               VEKTORNORM                                       *
!                                                                *
!----------------------------------------------------------------*
!                                                                *
!  BENOETIGTE UNTERPROGRAMME: COMABS, COMDIV                     *
!                                                                *
!*****************************************************************
!                                                                *
!  AUTOR     : JUERGEN DIETEL                                    *
!  DATUM     : 10.04.1987                                        *
!  QUELLCODE : FORTRAN 77                                        *
!                                                                *
!*****************************************************************
!
      INTEGER ND,N
      DOUBLE PRECISION V(ND,N),WI(N)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO = 0.0D0,ONE = 1.0D0)
      INTEGER I,J
      DOUBLE PRECISION MAXI,TR,TI,COMABS
!
      J = 1
   10 IF (J .GT. N) GOTO 80
         IF (WI(J) .EQ. ZERO) THEN
             MAXI = V(1,J)
             DO 15 I=2,N
   15           IF (ABS(V(I,J)) .GT. ABS(MAXI)) MAXI = V(I,J)
             IF (MAXI .NE. ZERO) THEN
                MAXI = ONE/MAXI
                DO 20 I=1,N
   20              V(I,J) = V(I,J)*MAXI
             ENDIF
             J = J+1
         ELSE
            TR = V(1,J)
            TI = V(1,J+1)
            DO 30 I=2,N
               IF (COMABS(V(I,J),V(I,J+1)) .GT. COMABS(TR,TI)) &
               THEN
                  TR = V(I,J)
                  TI = V(I,J+1)
               ENDIF
   30          CONTINUE
            IF (TR .NE. ZERO .OR. TI .NE. ZERO) THEN
               DO 40 I=1,N
   40             CALL COMDIV (V(I,J),V(I,J+1),TR,TI, &
                               V(I,J),V(I,J+1))
            ENDIF
            J = J+2
         ENDIF
         GOTO 10
   80 NORMAL = 0
      END
!
!
      SUBROUTINE SWAP(X,Y)
!
!*****************************************************************
!                                                                *
!     ZWECK DES PROGRAMMS:                                       *
!     ====================                                       *
!     SWAP VERTAUSCHT DIE WERTE DER BEIDEN DOUBLE PRECISION-     *
!     VARIABLEN X UND Y.                                         *
!                                                                *
!*****************************************************************
!
      DOUBLE PRECISION X,Y,TEMP
      TEMP = X
      X = Y
      Y = TEMP
      END
