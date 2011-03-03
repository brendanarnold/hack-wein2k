      SUBROUTINE REAN0 (NKK,KZZ,IFFT1,IFFT2,IFFT3,U)                                             
!                                                                       
!     SETS UP IN THE ARRAY U THE STEP FUNCTION FOR ALL K VECTORS        
!     USED IN REAN1 (REANALYSING THE FOURIER COEFFICIENTS)              
!                                                                       
!                                                                       
      use struct
      use parallel

      IMPLICIT REAL*8 (A-H,O-Z)
!
      INCLUDE 'param.inc'
!
      COMPLEX*16      U,IMAG,PHS,P                                   
!                                                                       
      DIMENSION       U(ifft1,ifft2,ifft3)
      DIMENSION       KZZ(3,NKK)

      integer   ifirst(nat)
      complex*16, allocatable:: ubuf(:,:,:)
      
#ifdef Parallel
  include 'mpif.h'
#endif
                                                            
      DATA ZILCH/1.D-10/,IMAG/(0.D0,1.D0)/                  
!--------------------------------------------------------------------   
!                                          
#ifdef Parallel
!  call MPI_BARRIER(comm_coul,ierr)
#endif

         u0=0.0d0

      index=1                      
      DO jatom=1,NAT
         ifirst(jatom)=index
         DO MU=1,MULT(JATOM)
            INDEX=INDEX+1
         enddo
      enddo
      
      PI=ACOS(-1.D0)
      PI2=2.D0*PI
!.....DEFINE STEPFUNCTION U FOR K=0      
! dirty fix for parallelization (U0 is summed over NPE)                       
      if(inversion)then
         do i3=1,ifft3
            ikz=i3-1
            if(i3.gt.ifft3/2) ikz=ikz-ifft3
            do i2=1,ifft2
               iky=i2-1
               if(i2.gt.ifft2/2) iky=iky-ifft2
               do i1=1,ifft1
                  ikx=i1-1
                  if(i1.gt.ifft1/2) ikx=ikx-ifft1
                  SK1=ikx
                  SK2=iky
                  SK3=ikz
                  SSK1=SK1*BR1(1,1)+SK2*BR1(1,2)+SK3*BR1(1,3)
                  SSK2=SK1*BR1(2,1)+SK2*BR1(2,2)+SK3*BR1(2,3)
                  SSK3=SK1*BR1(3,1)+SK2*BR1(3,2)+SK3*BR1(3,3)
                  S=SQRT(SSK1*SSK1+SSK2*SSK2+SSK3*SSK3)
                  IF(S.LT.ZILCH) THEN
!.....DEFINE STEPFUNCTION U FOR K=0
                     U(I1,I2,I3)=U0
!.....DEFINE STEPFUNCTION FOR K#0
                  ELSE
                     UTMP=0.D0
                     DO jatom=nstart(myid),nstop(myid)
                        PHSR=0.
                        index=ifirst(jatom)
                        DO MU=1,MULT(JATOM)
                           ARG1=(POS(1,INDEX)*SK1  +  &
                                POS(2,INDEX)*SK2   +  &
                                POS(3,INDEX)*SK3)*PI2
!                        ARG1=(ARG1+ARG2+ARG3)*PI2
                           PHSR=PHSR+Cos(Arg1)
                           INDEX=INDEX+1
                        enddo
                        Q=S*RMT(JATOM)
                        BESR=(Q*COS(Q)-SIN(Q))/(Q*Q*Q)
                        UTMP=UTMP+VSPH(JATOM)*BESR*PHSR
                     enddo
                     U(I1,I2,I3)=3.0D0*UTMP
                  ENDIF
               end do
            end do
         end do
      else
         do i3=1,ifft3
            ikz=i3-1
            if(i3.gt.ifft3/2) ikz=ikz-ifft3
            do i2=1,ifft2
               iky=i2-1
               if(i2.gt.ifft2/2) iky=iky-ifft2
               do i1=1,ifft1
                  ikx=i1-1
                  if(i1.gt.ifft1/2) ikx=ikx-ifft1
                  SK1=ikx
                  SK2=iky
                  SK3=ikz
                  SSK1=SK1*BR1(1,1)+SK2*BR1(1,2)+SK3*BR1(1,3)
                  SSK2=SK1*BR1(2,1)+SK2*BR1(2,2)+SK3*BR1(2,3)
                  SSK3=SK1*BR1(3,1)+SK2*BR1(3,2)+SK3*BR1(3,3)
                  S=SQRT(SSK1*SSK1+SSK2*SSK2+SSK3*SSK3)
                  IF(S.LT.ZILCH) THEN                                               
!                                                                       
!.....DEFINE STEPFUNCTION U FOR K=0                                     
                     U(I1,I2,I3)=U0                                                 
!                                                                       
!.....DEFINE STEPFUNCTION FOR K#0                                       
                  ELSE                                                              
                     U(I1,I2,I3)=0.D0
                     DO jatom=nstart(myid),nstop(myid)
                        PHS=(0.D0,0.D0) 
                        index=ifirst(jatom)
                        DO MU=1,MULT(JATOM)
                           ARG1=(POS(1,INDEX)*SK1 + &
                                POS(2,INDEX)*SK2 + &
                                POS(3,INDEX)*SK3)*PI2
!                        ARG1=(ARG1+ARG2+ARG3)*PI2
                           PHS=PHS+DCMPLX(Cos(Arg1),sin(Arg1))
                           INDEX=INDEX+1                                               
                        enddo
                        Q=S*RMT(JATOM)
                        BESR=(Q*COS(Q)-SIN(Q))/(Q*Q*Q)                                  
                        U(I1,I2,I3)=U(I1,I2,I3)+VSPH(JATOM)*BESR*PHS
                     enddo
                     U(I1,I2,I3)=3.0D0*U(I1,I2,I3)                                     
                  ENDIF
               end do
            end do
         end do
      endif

#ifdef Parallel      
      allocate(ubuf(ifft1,ifft2,ifft3))
      ubuf=u
      call mpi_allreduce(ubuf,u,ifft1*ifft2*ifft3, MPI_DOUBLE_COMPLEX, MPI_SUM, comm_coul, ierr)
      deallocate(ubuf)
#endif

      U0=VOL
      DO JATOM=1,nat
         U0=U0-VSPH(JATOM)*MULT(JATOM)                                  
      END DO
      u(1,1,1)=u0

      RETURN                                                            
      END                                                               
