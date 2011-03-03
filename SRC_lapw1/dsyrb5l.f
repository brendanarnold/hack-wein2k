!
       SUBROUTINE DSYRB5L( BI, LENO, WIDTHO, W, Y, V, WORK, &
                           LWORK, A, LDA, RESTA, TAU ) 
!
! Parameters:
!
       implicit            none
       INTEGER             BI, LENO, WIDTHO, LWORK, LDA, RESTA
!       DOUBLE PRECISION    W( LDA, WIDTHO ), 
!     $                     Y( LDA, WIDTHO ), 
!     $                     V( LDA, WIDTHO ),
       DOUBLE PRECISION    W( RESTA, WIDTHO ),  &
                           Y( RESTA, WIDTHO ),  &
                           V( RESTA, WIDTHO ), &
                           WORK( LWORK ),  &
                           A( LDA, RESTA ),  &
                           TAU( WIDTHO )
!
! Local variables:
!
       INTEGER             KO, NI, J, SI, RI, LENI, WIDTHI,  &
                           KI, K, I
!
! Local constants:
!
       DOUBLE PRECISION    ZERO, ONE
       PARAMETER           ( ZERO = 0.0D0, ONE = 1.0D0 )
!
!

!           KO is the number of already computed columns of the current outer block
            KO = 0
       
!           NI is the number of inner blocks ( =ceiling((RO-SO+1)/BI) ) 
            NI = (WIDTHO - 1)/BI + 1
	    DO J = 1, NI

!             SI is the first column of the current inner block
	      SI = BI * ( J - 1 ) + 1
!             RI is the last column of the current inner block
	      RI = MIN( WIDTHO, BI * J )

!             length of the current inner block
	      LENI = LENO - SI + 1
!             width of the current inner block
              WIDTHI = RI - SI + 1

              IF( KO .GT. 0 ) THEN
!               update current inner block with existing parts of 
!               current outer block (rank-2k update)
!               A(inner block) = A(inner block) + V Y' + Y V'
!              write(6,*) 'TAU vor dblr2k: ', (tau(k), k=si,si+widthi-1)
!	      write(6,*) (a(i,si),i=si, si+4)
                CALL DBLR2K( LENI + BI, WIDTHI, KO,  &
                             V( SI, 1 ), RESTA, &
                             Y( SI, 1 ), RESTA, &
                             A( SI, SI ), LDA )
              END IF

!             QR factorization of current inner block
!               compute Householder vectors and apply them from left
!               to rest of inner block
!              write(6,*) 'TAU vor dgeqrl: ', (tau(k), k=si,si+widthi-1)
!	      write(6,*) (a(i,si),i=si+bi, si+bi+4)
              CALL DGEQRL( 'R', LENI, WIDTHI, A( SI + BI, SI ), &
                    LDA, ZERO, TAU( SI ), WORK )
!              write(6,*) 'ZERO = ', ZERO

!              write(6,*) 'TAU nach dgeqrl: ', (tau(k), k=si,si+widthi-1)
!             clean unused parts of W and Y and V
              DO K = 1, WIDTHI
                CALL DCOPY( SI + BI, ZERO, 0, W( 1, KO + K), 1)
                CALL DCOPY( SI + BI, ZERO, 0, Y( 1, KO + K), 1)
                CALL DCOPY( SI + BI, ZERO, 0, V( 1, KO + K), 1)
              END DO

!             compute W and Y to get Q=I+WY', using the Householder
!             vectors computed in DGEQRL 
!             (Y holds a copy of the Householder vectors stored in A)
!              write(6,*) 'TAU vor dgewyg: ', (tau(k), k=si,si+widthi-1)
              CALL DGEWYG( 'L', LENI, WIDTHI, A( SI + BI, SI ), &
                           LDA, TAU( SI ),  &
                           W( SI + BI, KO + 1 ), RESTA,  &
                           Y( SI + BI, KO + 1 ), RESTA, KI, WORK )
!              write(6,*) 'TAU nach dgewyg: ', (tau(k), k=si,si+widthi-1)
!	      write(6,*) (a(i,si),i=si+bi, si+bi+4)
              IF ( KI .NE. WIDTHI ) THEN
	        write(6,*) 'KI = ', KI, '  WIDTHI = ', WIDTHI
!                write(6,*) 'Problem in dsyrb5l: degenerated matrix   '
                STOP 'Problem in dsyrb5l: degenerated matrix   '
              END IF

!             find V, so that Q'AQ becomes A+VY'+YV'
	      CALL DSYWYV( 'L', LENI, KO, KI,  &
                           A( SI + BI, SI + BI ), LDA,  &
                           W( SI + BI, KO + 1 ), RESTA,  &
                           Y( SI + BI, KO + 1 ), RESTA, &
                           Y( SI + BI, 1 ), RESTA,  &
                           V( SI + BI, 1 ), RESTA, &
                           V( SI + BI, KO + 1 ), RESTA,  &
                           WORK, MAX(K,KI)*KI )
!              write(6,*) 'TAU nach dsywyv: ', (tau(k), k=si,si+widthi-1)
!	      write(6,*) (a(i,si+bi),i=si+bi, si+bi+4)
                

!             number of (real) columns of Y, W, and V
	      KO = KO + KI

   90       END DO

!
!           Update rest of the matrix
!

!           Rank-2k update of symmetric part
            CALL DSYR2K( 'L', 'N', RESTA - WIDTHO, KO, &
                         ONE, V( WIDTHO + 1, 1), RESTA,  &
                              Y( WIDTHO + 1, 1), RESTA, &
              ONE, A( WIDTHO + 1, WIDTHO + 1), LDA )
!              write(6,*) 'TAU nach dsyr2k: ', (tau(k), k=si,si+widthi-1)
!	      write(6,*) (a(i,widtho+1),i=widtho+1, widtho+4)
!           One-sided update with last inner block
            IF( ( BI - WIDTHI ) .GT. 0) THEN
              CALL DGEWY( 'L', RESTA - WIDTHO + WIDTHI - BI,  &
                          BI - WIDTHI, WIDTHI,  &
                          A( WIDTHO + BI - WIDTHI + 1, WIDTHO + 1), LDA, &
                          W( SI + BI, KO - KI + 1), RESTA, &
                          Y( SI + BI, KO - KI + 1), RESTA, &
                          WORK )
!              write(6,*) 'TAU nach dgewy: ', (tau(k), k=si,si+widthi-1)
!	      write(6,*) (a(i,widtho+1),i=widtho+bi-widthi+1, 
!     +                    widtho+bi-widthi+4)
            END IF

       END 

