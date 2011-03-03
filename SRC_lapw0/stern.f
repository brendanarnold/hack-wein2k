      SUBROUTINE STERN (JJA,NST,STG,TAUP)                               
!                                                                       
!     STERN GENERATES THE STAR OF REC LATTICE VECTOR KZZ(1,JJA)         
!     THE STAR VECTORS ARE STORED IN STG, THE STAR-SIZE IN  NST         
!     IZ CONTAINS THE SYMMETRY-MATRICES                                 
!                                                                       
!                                                                       
      use struct, only : NSYM,inversion
      use symetr
      IMPLICIT REAL*8 (A-H,O-Z)
!
      INCLUDE 'param.inc'
!
      COMPLEX*16          TAUP,IMAG                                     
      INTEGER          G,STG,INDEX(NSYM)                                
!                                                                       
      DIMENSION        STG(3,NSYM),TAUP(NSYM),G(3)                      
!      DATA IMAG        /(0.D0,1.D0)/                                    
      PARAMETER (TPI= 6.283185307179586231995926937088D0)
!-------------------------------------------------------------------    
!                                                                       
!      TPI=2.D0*ACOS(-1.D0)                                              
      G(1)=KZZ(1,JJA)                                                   
      G(2)=KZZ(2,JJA)                                                   
      G(3)=KZZ(3,JJA)                                                   
      NST=0                                                             
!                                                                       
!.....START LOOP OVER ALL SYMMETRY OPERATIONS                           
      DO 1 I=1,IORD                                                     
         TK=0.                                                          
         DO 2 J=1,3                                                     
         TK=TK+TAU(J,I)*G(J)*TPI                                     
           K=0                                                         
           DO 3 L=1,3                                                  
   3       K=IZ(J,L,I)*G(L)+K                                          
!   3       K=IZ(l,j,I)*G(L)+K                                          
           STG(J,I)=K                                                     
!        TK=TK-TAU(J,I)*STG(J,I)*TPI                                     
 2      continue
!      if (jja.lt.5 ) then
!      write (6,768) jaa,i,(stg(j,i),j=1,3), tk
! 768  format(2i3,3i5,2f10.3) 
!      end if                                           

         IF(NST.EQ.0) GOTO 7                                            
!                                                                       
!........PROOF, IF THE VECTOR STG(J,I) IS A NEW STARMEMBER OR NOT       
         DO 4 M=1,NST                                                   
            DO 5 J=1,3                                                  
               IF(STG(J,M).NE.STG(J,I)) GOTO 4                          
   5        CONTINUE                                                    
!           STG(J,I) IS NOT A NEW STARMEMBER, IT IS EQUIV TO STG(J,M).  
!           BUT THE TAUPHASE OF STG(J,I) CAN BE NEW.  THEREFORE THE     
!           ALREADY DEFINED PHASE TAUP(M) IS AVERAGED WITH THE PHASE    
!           OF STG(J,M)                                                 
!            TAUP(M)=TAUP(M)+EXP(TK*IMAG)                                
            TAUP(M)=TAUP(M)+DCMPLX(cos(TK),sin(TK))
            INDEX(M)=INDEX(M)+1                                         
            GOTO 1                                                      
   4     CONTINUE                                                       
!                                                                       
!........NEW VECTOR FOUND ]]]                                           
   7     NST=NST+1                                                      
         DO 6 J=1,3                                                     
   6     STG(J,NST)=STG(J,I)                                            
!         TAUP(NST)=EXP(TK*IMAG)                                         
         TAUP(NST)=DCMPLX(cos(TK),sin(TK))
         INDEX(NST)=1                                                   
   1  CONTINUE                                                          
!                                                                       
      DO 10 I=1,NST
        if(inversion)then
                TAUP(I)=DBLE(TAUP(I))/dble(INDEX(I))   
        else           
                TAUP(I)=TAUP(I)/dble(INDEX(I))                                          
        endif
!      if (jja.lt.10 ) then
!      write (6,767) jaa,(stg(j,i),j=1,3), taup(i)
! 767  format(i3,3i5,2f10.3) 
!      end if                                           
  10  continue
      RETURN                                                            
      END                                                               
	subroutine haveinversion
!	Does the structure file have inversion?
	USE struct
        USE symetr
	inversion=.false.
!        write(6,*)'Searching over ', iord,' operations'
  	DO j=1,iord
!     	Check if inversion is present
      	if((IZ(1,1,J).eq.-1).and.(IZ(2,2,J).EQ.-1).and.(IZ(3,3,J).EQ.-1))then
!     	Diagonal elements OK, check others
         if((IZ(2,1,J).EQ.0).and.(IZ(3,1,J).EQ.0).and.(IZ(1,2,J).EQ.0).and. &
          (IZ(3,2,J).EQ.0).and.(IZ(1,3,J).EQ.0).and.(IZ(2,3,J).EQ.0))then
!     	Off-diagonals OK, insist on zero translations (may not be always correct)
          if((abs(tau(1,J)).lt.1d-5).and.(abs(tau(2,J)).lt.1d-5).and.(abs(tau(3,J)).lt.1d-5))then
!
!     	We have inversion
               inversion=.true.
          endif
         endif
      	endif
  	ENDDO
!        write(6,*)'Found inversion as',inversion
  	return
  	end

