      PROGRAM getepar                                                     
!
      DOUBLE PRECISION E
      INTEGER N,I
      implicit real*8 (a-h,o-z)                 
      CHARACTER*40 FNAME
      CHARACTER*20 FFORM
      call NAMEFI(FNAME,FFORM)
      call fiopen(fname,8)
      I=1
      N=0
      DO
         N=N+1
         READ(IO,123)B,A
         WRITE(*,*)A
         IF(I.EQ.1) then
            
            IF(A.GT.0) THEN
            I=I+1   
            WRITE(*,*)A
            ENDIF
         ENDIF

         IF(I.EQ.2) then
            IF(A.EQ.0) THEN
            I=I+1

            WRITE(*,*)A
            ENDIF
         ENDIF
         IF(I.EQ.3) then
            IF(A.GT.0) THEN
            WRITE(*,*)A
            STOP'END'
            ENDIF
         ENDIF
      ENDDO
 123  FORMAT(f10.5,f10.6)
 3000 END                                              
!



*********************************************************************
      SUBROUTINE NAMEFI(FNAME,FFORM)
!     
!     EINGABE DES FILENAMENS UND DESSEN FORMATS
!
      logical ex
      CHARACTER*40 FNAME,HELPN
      CHARACTER*20 FFORM,HELPF 
 101  WRITE(*,'(A)')' NEW FILENAME: '
      READ(7,'(A40)',ERR=100) HELPN
      IF(HELPN.NE.' ') THEN
        FNAME= HELPN
        inquire(file=fname,exist=ex)
        if(.not.ex) then
          write(*,'(a)')' file does not exist!!'
          goto 101
        end if
      end if

      RETURN
100   WRITE(*,'(A/A/A)') ' I THINK, YOU MADE A MISTAKE !', &
                         ' PLEASE, TRY AGAIN.',  &
                         ' AS USUAL ENTER RETURN TO CONTINUE.'
      READ(7,'(F2.0)') XX
      goto 101
      END
************************************************************************
      SUBROUTINE FIOPEN(FNAME,IO)
!
!     OPEN FILE
!
      CHARACTER*40 FNAME
      CLOSE(IO,ERR=9)

9     CONTINUE
      OPEN(IO,FILE=FNAME,STATUS='UNKNOWN',ERR=10)

      RETURN
10    WRITE(*,*) 'ERROR IN OPENING FILE. ENTER RETURN TO CONTINUE'
      READ (7,'(F2.0)') XX
      RETURN
      END
************************************************************************





