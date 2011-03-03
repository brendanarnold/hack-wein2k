      subroutine nose(xmin1,x0,xplus1,xnosmi,xnos0,xnospl,n,fwert, &
                      fnose,temp,delta,weight,mult)
      implicit double precision (a-h,o-z)
      include 'param.inc'
      COMMON /STRUK/  AA,BB,CC,ALPHA(3),PIA(3),VOL
      DIMENSION  xmin1(n),x0(n),xplus1(n),fwert(n)
      DIMENSION  weight(*),mult(*)
!
!     initialize nose 
!
      write(6,*)'input nose:temp,delta,fnose',temp,delta,fnose
      write (6,*) (xmin1(k),x0(k),weight(k),k=1,n)
      fnose=fnose*1.0d12/(4.1341373358254d16)
      boltzk=3.166678911d-6
      API=4.D0*DATAN(1.D0)
      free=-6.0d0
      do 10 i=1,n,3
         help=dble(i+2)/3.0d0
         iatom=nint(help)
         do 20 j=1,mult(iatom)
            free=free+3.0d0
  20     CONTINUE
  10  CONTINUE
      if (free.lt.0.1d0) free=1.0d0
      if (xnosmi.gt.10000.0d0) then
      do 110 i=1,n,3
         help=dble(i+2)/3.0d0
         iatom=nint(help)
!            vmax=1.0d0/0.70d0*dsqrt(boltzk*free*temp/weight(iatom))
       vmax=dsqrt(1.0d0/2.0d0*boltzk*temp/weight(iatom))
            xmin1(i)=x0(i)+gauss(dummy)*vmax*delta
            xmin1(i+1)=x0(i+1)+gauss(dummy)*vmax*delta
            xmin1(i+2)=x0(i+2)+gauss(dummy)*vmax*delta
            write(6,*)'vmax,random positions',iatom,vmax, &
                       xmin1(i),xmin1(i+1),xmin1(i+2)
 110  CONTINUE
      endif
      call ininos(xmin1,x0,delta,ekin,n,weight,mult)
      temp1=2.0d0*ekin/boltzk/free
      write(6,*)'old ekin,temp,xnos0,xnosmi,free', &
                  ekin,temp1,xnos0,xnosmi,free
      QNOSE=2.0d0/((2.0d0*API*FNOSE)**2)*2.0d0*FREE*boltzk*TEMP
      if (xnosmi.gt.10000.0d0) then
         VNOS=delta/qnose*(ekin-FREE*boltzk*TEMP/2.0d0)
         facnos=1.0d0/(1.0d0+VNOS*delta/2.0d0)
         xnospl=(delta**2)/qnose* &
                (2.0d0*ekin*(facnos**2)-FREE*boltzk*TEMP)
         ekin=ekin*(facnos**2)
         anner=(VNOS*delta)/2.0d0
            write(6,*) 'nose:i,vnos,facnos',i,vnos,facnos
      else
         VNOS=(xnos0-xnosmi)/delta+ &
              delta/qnose*(ekin-FREE*boltzk*TEMP/2.0d0)
         facnos=1.0d0/(1.0d0+VNOS*delta/2.0d0)
!
!     determine anner
!
         do 30 i=1,100
            xnospl=2.0d0*xnos0-xnosmi+(delta**2)/qnose* &
                   (2.0d0*ekin*(facnos**2)-FREE*boltzk*TEMP)
            VNOS=(xnospl-xnosmi)/2.0d0/delta
            facnos=1.0d0/(1.0d0+VNOS*delta/2.0d0)
            write(6,*) 'nose:i,vnos,facnos',i,vnos,facnos, &
            (2.0d0*ekin*(facnos**2)-FREE*boltzk*TEMP),(delta**2)/qnose
   30    CONTINUE
         ekin=ekin*(facnos**2)
         VNOS=(xnospl-xnosmi)/2.0d0/delta
         anner=(VNOS*delta)/2.0d0
         write(6,*)'nose:ekin,vnos,free,anner',ekin,vnos,free,anner
      endif
!
!     force and gradient in Ry/bohr
!     determine new position
!
      iatom=1
      ihelp=0
      do 40 i=1,n
         ihelp=ihelp+1
         xplus1(i)=2.0d0*x0(i)/(1.0d0+anner)+ &
                   xmin1(i)*(1.0d0-2.0d0/(1.0d0+anner))+ &
                   fwert(i)/2.0d0* &
                   (delta**2)/weight(iatom)/(1.0d0+anner)
!         write(6,'(10x,''(nwtmin)  i,damp,speed,force =''
!     &      ,i3,f6.2,2f10.4)') i,damp,speed,fwert(i)
         if (ihelp.eq.3) then
            ihelp=0
            iatom=iatom+1
         endif
 40   CONTINUE
      call ininos(x0,xplus1,delta,ekin,n,weight,mult)
      temp=2.0d0*ekin/boltzk/free
      write(6,*)'  new position'
      write(6,*)(xplus1(i),i=1,n)
      write(6,*)'fwert',(fwert(i),i=1,n)
      write(6,*)'new ekin,temp,xnospl',ekin,temp,xnospl
550   format( 20x,3f12.4 )
      return
      end
!    *******************************************************************
!    ** THIS FORTRAN CODE IS INTENDED TO ILLUSTRATE POINTS MADE IN
!    ** THE TEXT. TO OUR KNOWLEDGE IT WORKS CORRECTLY. HOWEVER IT IS
!    ** THE RESPONSIBILITY OF THE USER TO TEST IT, IF IT IS USED IN A
!    ** RESEARCH APPLICATION.
!    *******************************************************************
!    *******************************************************************
!    ** FICHE F.24.  INITIAL VELOCITY DISTRIBUTION
!    *******************************************************************
!    *******************************************************************
!    ** CENTRE OF MASS AND ANGULAR VELOCITIES FOR LINEAR MOLECULES


        REAL*8 FUNCTION GAUSS ( DUMMY )

!    *******************************************************************
!    ** RANDOM VARIATE FROM THE STANDARD NORMAL DISTRIBUTION.
!    **
!    ** THE DISTRIBUTION IS GAUSSIAN WITH ZERO MEAN AND UNIT VARIANCE.**
!    **
!    ** REFERENCE:
!    **
!    ** KNUTH D, THE ART OF COMPUTER PROGRAMMING, (2ND EDITION
!    **    ADDISON-WESLEY), 1978
!    **
!    ** ROUTINE REFERENCED:
!    **
!    ** REAL FUNCTION RANF ( DUMMY )
!    **    RETURNS A UNIFORM RANDOM VARIATE ON THE RANGE ZERO TO ONE
!    *******************************************************************

        REAL*8        A1, A3, A5, A7, A9
        PARAMETER ( A1 = 3.949846138d0, A3 = 0.252408784d0 )
        PARAMETER ( A5 = 0.076542912d0, A7 = 0.008355968d0 )
        PARAMETER ( A9 = 0.029899776d0                   )

        REAL*8      SUM, R, R2
        REAL*8      RANF, DUMMY
        INTEGER     I

!    *******************************************************************

        SUM = 0.0d0

        DO 10 I = 1, 12

           SUM = SUM + RANF ( DUMMY )

10      CONTINUE

        R  = ( SUM - 6.0d0 ) / 4.0d0
        R2 = R * R

        GAUSS = (((( A9 * R2 + A7 ) * R2 + A5 ) * R2 + A3 ) * R2 +A1 ) &
                * R

        write(6,*) 'gauss',gauss
        RETURN
        END
        REAL*8 FUNCTION RANF ( DUMMY )

!    *******************************************************************
!    ** RETURNS A UNIFORM RANDOM VARIATE IN THE RANGE 0 TO 1.
!    **
!    **                 ***************
!    **                 **  WARNING  **
!    **                 ***************
!    **
!    ** GOOD RANDOM NUMBER GENERATORS ARE MACHINE SPECIFIC.
!    ** PLEASE USE THE ONE RECOMMENDED FOR YOUR MACHINE.
!    *******************************************************************

        INTEGER     L, C, M
        PARAMETER ( L = 1029, C = 221591, M = 1048576 )

        INTEGER     SEED
        REAL*8      DUMMY
        SAVE        SEED
        DATA        SEED / 0 /

!    *******************************************************************

        SEED = MOD ( SEED * L + C, M )
        RANF = DBLE ( SEED ) / M

        RETURN
        END



