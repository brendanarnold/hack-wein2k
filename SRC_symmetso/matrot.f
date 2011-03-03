      subroutine matrot(znor,xpri,ifotro,cb,rotold,label)
      implicit real*8 (a-h,o-z)
!  This program determines the MATRIX that rotates the axis
!        x1=(1,0,0)  
!        y1=(0,1,0)     
!        z1=(0,0,1)     
!  into the ones x2,y2,z2 where z2 is read and x2 and y2 are 
!  determined in the following way, depending on the value of
!  the option parameter IFOTRO:
!
!      if IFOTRO=1: x2 is read, y2=z2^x2
!      if IFOTRO=2: y2 is read, x2=y2^z2
!      if IFOTRO=0: x2 is determined by the program to be ortogonal
!                   to z2, y2=z2^x2
!
!----------------------------------------------------------------
!  The Loc.Rot.Mat. is in the format required by *.struct file
! 
!  Revised 12/11/96. 
!
!
      DIMENSION ICI(5),CB(3,3),CB2(3,3),rotold(3,3) &
               ,ZNOR(3),XPRI(3),YPRI(3)
      logical equal

      ICI(1)=1
      ICI(2)=2
      ICI(3)=3
      ICI(4)=1
      ICI(5)=2         

! ##################################################   

!      open(unit=01,file='inp',status='old')
!      open(unit=02,file='outp',status='new')

!      READ(1,*)
      write(6,15)       (ZNOR(ix),ix=1,3)
      write(6,16)       (XPRI(ix),ix=1,3),IFotro
 15   format('  z-rotation vector:',3f8.4)
 16   format('  y-rotation vector:',3f8.4,i5)
!  Defino un nuevo sistema S2 de ejes: ZNOR, XPRI, YPRI

      ZZ=XXMOD(ZNOR)
      DO I=1,3
        ZNOR(I)=ZNOR(I)/ZZ
      END DO
!...................................IFotro=0                         
      IF(IFotro.EQ.0)THEN
         LCUEN=0
         DO I=1,3
           QQ=ABS(ZNOR(I))
           IF(QQ.GT.0.000001) THEN
              IZNONU=I
              LCUEN=LCUEN+1
           END IF
         END DO          

         I1=ICI( IZNONU )
         I2=ICI( IZNONU + 1 )
         I3=ICI( IZNONU + 2 )

! Si ZNOR tiene una sola componente no nula...
         IF(LCUEN.EQ.1) THEN
            XPRI(I1)=0.0
            XPRI(I2)=1.0
            XPRI(I3)=0.0
! Si hay mas de una comp. de ZNOR no nula, elijo XPRI de modo que
! sea ortogonal a ZNOR:
                        ELSE
            XPRI(I1)= 0.0
            XPRI(I2)= ZNOR(I3)
            XPRI(I3)=-ZNOR(I2)
         END IF

         XX=XXMOD(XPRI)
         DO I=1,3
           XPRI(I)=XPRI(I)/XX
         END DO
         
         DO I=1,3
           YPRI(I)=V(I,ZNOR,XPRI)
         END DO
                     ELSE
         PESC=0.
         DO I=1,3
            PESC=PESC+ZNOR(I)*XPRI(I)
         END DO
         IF(PESC.GT.0.000001)THEN
            IF(IFotro.EQ.1) WRITE(6,*)' x2 is not ortogonal to z2'
            IF(IFotro.EQ.2) WRITE(6,*)' y2 is not ortogonal to z2'
            write(6,*)' COULD NOT DETERMINE LOCAL ROTATION MATRIX!!!'
            label=label+1
            return
!ccc            STOP
         END IF
      END IF                 
!...................................IFotro=0                         

!...................................IFotro=1                         
      IF(IFotro.EQ.1)THEN
      XX=XXMOD(XPRI)
      DO I=1,3
        XPRI(I)=XPRI(I)/XX
      END DO
         
      DO I=1,3
        YPRI(I)=V(I,ZNOR,XPRI)
      END DO       
      END IF                      
!...................................IFotro=1                         


!...................................IFotro=2                         
      IF(IFotro.EQ.2)THEN
      do i=1,3
        YPRI(i)=XPRI(i)
      end do
      XX=XXMOD(YPRI)
      DO I=1,3
        YPRI(I)=YPRI(I)/XX
      END DO
         
      DO I=1,3
        XPRI(I)=V(I,YPRI,ZNOR)
      END DO
      END IF                             
!...................................IFotro=2                         

      DO I2=1,3
        CB(I2,1)=XPRI(I2)
        CB(I2,2)=YPRI(I2)
        CB(I2,3)=ZNOR(I2)      
      END DO

      CALL INVA(CB,CB2,IERR)

      IF(IERR.EQ.1)THEN
         WRITE(6,*)'WARNING !!! DET(CB)=0'
            write(6,*)' COULD NOT DETERMINE LOCAL ROTATION MATRIX!!!'
            label=label+1
            return
!ccc         STOP
      END IF
      equal=.true.
      do i=1,3
      do j=1,3
      if(abs(cb(i,j)-rotold(i,j)).gt.1.d-4) equal=.false.
      enddo
      enddo
      if(.not.equal) write(6,*)' WARNING: LOCAL ROTATION MATRIX CHANGED'      
      write(6,12)
      write(6,10) (CB(1,j),j=1,3),(rotold(1,j),j=1,3)
      do i=2,3
      write(6,11) (CB(i,j),j=1,3),(rotold(i,j),j=1,3)
      end do

   10 FORMAT(10x,3F10.7,5x,3f10.7)                        
   11 FORMAT(10X,3F10.7,5x,3f10.7)                             
 12   FORMAT('LOCAL ROT MATRIX:     ','  NEW',30x,'  OLD')
!      close (unit=01)
!      close (unit=02)
!      STOP
      END

!$23
      FUNCTION XXMOD(VEC)
      implicit real*8 (a-h,o-z)
      DIMENSION VEC(3)
      
      XXMOD=SQRT( VEC(1)**2 + VEC(2)**2 + VEC(3)**2)

      RETURN
      END   


!$24
      SUBROUTINE INVA(A,AINV,IERR)
      implicit real*8 (a-h,o-z)
! Este programa calcula la inversa de una matriz de 3*3 
                                                        
      DIMENSION A(3,3),AEXT(5,5),ADJA(3,3),AINV(3,3),ICI(5)

      IERR=0

      ICI(1)=1
      ICI(2)=2
      ICI(3)=3
      ICI(4)=1
      ICI(5)=2
! Extiendo A auna matriz de 5x5 
      DO I=1,5
      DO J=1,5
        AEXT(I,J)=A(ICI(I),ICI(J))
      END DO
      END DO
! Calculo la Adjunta de A
      DO I=1,3 
      DO J=1,3
      ADJA(I,J)=AEXT(I+1,J+1)*AEXT(I+2,J+2)-AEXT(I+1,J+2)*AEXT(I+2,J+1)
      END DO
      END DO
! Calculo el determinante de A
      DETA=0.0
      DO I=1,3
        DETA=DETA + A(1,I)*ADJA(1,I)
      END DO

      QQ=ABS(DETA)
      IF(QQ.LT.0.000001)THEN
        IERR=1
        RETURN
      END IF
! Calculo la inversa de A :
      DO I=1,3
      DO J=1,3
        AINV(I,J)=ADJA(J,I)/DETA
      END DO
      END DO
 
      RETURN
      END        
      




      FUNCTION V(I,A,B)
      implicit real*8 (a-h,o-z)
!  V(i,A,B)=(AxB)coord i
      DIMENSION A(3),B(3)  
      IF(I.EQ.1) V = A(2)*B(3) - A(3)*B(2)
      IF(I.EQ.2) V = A(3)*B(1) - A(1)*B(3)
      IF(I.EQ.3) V = A(1)*B(2) - A(2)*B(1)
      RETURN
      END




