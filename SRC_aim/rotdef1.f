      SUBROUTINE ROTDEF(NAT,MULT,IOP,POS,LATTIC)                               
!     
!     FINDS MATRICES, WICH TRANSFORM "index-ATOM" ONTO "main-ATOM"
!     WARNING:   opposite to usual rotdef !!!!      
      IMPLICIT REAL*8 (A-H,O-Z)
	character*4 lattic
      real*8 invoiz
      INCLUDE 'param.inc'
      COMMON /SYM2/ tau(3,NSYM),otau(3,NSYM),oiz(3,3,NSYM), &
         invoiz(3,3,NSYM),iz(3,3,NSYM),iord
      DIMENSION POS(3,*),IOP(*),MULT(NAT)                         
!      DATA TOLER/ 1.D-4/,ONE/1.0D0/                                       
      parameter (one=1.D0, toler=1.D-4)
      toler2=toler/2.d0
      INDEX=0                                                           
      NCOUNT=0                                                          

      DO 20 JATOM=1,NAT                                                 
         INDEX1=INDEX+1                                                 
         DO 30 M=1,MULT(JATOM)                                          
            INDEX=INDEX+1                                               
            X=0.D0
            Y=0.D0
            Z=0.D0
            DO 25 I=1,IORD
               X=0.D0                                                        
               DO 26 J=1,3
 26               X=X+IZ(J,1,I)*POS(J,INDEX)
                  X=X+TAU(1,I) + 5.D0
                  X= MOD (X+toler2,ONE)-toler2
                  Y=0.D0                                                        
                  DO 27 J=1,3
 27                  Y=Y+IZ(J,2,I)*POS(J,INDEX)
                     Y=Y+TAU(2,I) + 5.D0
                     Y= MOD (Y+toler2,ONE)-toler2
                     Z=0.D0          
                     DO 28 J=1,3
 28                     Z=Z+IZ(J,3,I)*POS(J,INDEX)
                        Z=Z+TAU(3,I) + 5.D0
                        Z= MOD (Z+toler2,ONE)-toler2
                        X1=ABS(X-POS(1,INDEX1))
                        Y1=ABS(Y-POS(2,INDEX1))
                        Z1=ABS(Z-POS(3,INDEX1))
!     WRITE(*,*) 'JATOM,INDEX1,INDEX,I',JATOM,INDEX1,INDEX,I       
!!      WRITE(*,*) X1,Y1,Z1                                        
!!!     WRITE(*,*) 'TRY', X,Y,Z                                        
                        IF(X1.LT.TOLER.AND.Y1.LT.TOLER.AND.Z1.LT.TOLER) &
                             THEN        
!         WRITE(*,*) 'FOUND MATCH FOR JATOM=',JATOM,' MULT=',M
                           GOTO 60
                        END IF
!....check positions for centered lattices
            if(lattic(1:1).eq.'B') then
              x1=mod(x1+0.500000001d0,one)
              y1=mod(y1+0.500000001d0,one)
              z1=mod(z1+0.500000001d0,one)
              IF(X1.LT.TOLER.AND.Y1.LT.TOLER.AND.Z1.LT.TOLER) goto 60
            endif                                    
            if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CXY') then
              x1=mod(x1+0.500000001d0,one)
              y1=mod(y1+0.500000001d0,one)
              IF(X1.LT.TOLER.AND.Y1.LT.TOLER.AND.Z1.LT.TOLER) goto 60
              x1=mod(x1+0.5d0,one)
              y1=mod(y1+0.5d0,one)
            endif                 
            if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CXZ') then
              x1=mod(x1+0.500000001d0,one)
              z1=mod(z1+0.500000001d0,one)
              IF(X1.LT.TOLER.AND.Y1.LT.TOLER.AND.Z1.LT.TOLER) goto 60
              x1=mod(x1+0.5d0,one)
              z1=mod(z1+0.5d0,one)
            endif                 
            if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CYZ') then
              y1=mod(y1+0.500000001d0,one)
              z1=mod(z1+0.500000001d0,one)
              IF(X1.LT.TOLER.AND.Y1.LT.TOLER.AND.Z1.LT.TOLER) goto 60
              y1=mod(y1+0.5d0,one)
              z1=mod(z1+0.5d0,one)
            endif                 

 25                  CONTINUE
                     WRITE(*,*) 'NO SYMMETRY OPERATION FOUND IN ROTDEF'
                     WRITE(*,*) JATOM,INDEX
!             WRITE(*,*) 'FOR SITE ',X,Y,Z 
                     WRITE(*,*) (POS(I1,JATOM),I1=1,3)
                     WRITE(*,*) (POS(I1,INDEX),I1=1,3)
                     STOP 'ROTDEF'
                     
!     
 60 		NCOUNT=NCOUNT+1
		IOP(INDEX)=I	
 30               CONTINUE
 20            CONTINUE
               IF(NCOUNT.NE.INDEX) THEN
                  WRITE(6,1000) NCOUNT
                  STOP ' ROTDEF: NCOUNT NE INDEX'
               ENDIF                             
               RETURN
 1000          FORMAT(///,3X &
                    ,'ERROR IN ROTDEF: ROTIJ NOT DEFINED FOR ALL ' &
                    ,'ATOMS OF BASIS',/,20X,'NCOUNT=',I2)
               END                                    
