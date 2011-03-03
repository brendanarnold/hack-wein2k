      SUBROUTINE REDUCE_ALLOC_R8_2D(N1OLD,N2OLD, &
                                     N1NEW,N2NEW)
!..REDUCES UP TO 4 DIMENSIONAL REAL*8 ARRAYS FROM OLD TO NEW ALLOCATIONS
!
      use variable_fields,  a => pos 
!      REAL*8 allocatable  :: A(:,:)
!      IF(NDIM.EQ.1) THEN
      real*8,  allocatable ::  HELP(:,:)
      ALLOCATE ( HELP(N1NEW,n2new) )
      DO I=1,N1NEW
      DO j=1,N2NEW
        HELP(I,j)=A(I,j)
      ENDDO
      ENDDO
      DEALLOCATE ( A )
      ALLOCATE ( A(N1NEW,n2new) )
      DO I=1,N1NEW
      DO j=1,N2NEW
        A(I,j)=HELP(I,j)
      ENDDO
      ENDDO
      DEALLOCATE ( HELP )      
      END

 
      
