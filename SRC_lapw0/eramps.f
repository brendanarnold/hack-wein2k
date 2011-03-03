!       Code to generate various types of potential ramps
!       Controlled by two real*8 parameters E and lambda
!       E is the strength, coded as refeld
!
!       Mode 11:        Triangular Ramp of E/2 to -E/2
!       Mode 12:        Triangular Ramp of 0 to -E
!       Mode 1:         =0 if |z|<lambda, otherwise E*cos((z0-lambda)*pi2/(0.5-lambda))
!       Mode 2:         =0 if |Z|<lambda, otherwise E*(1.D0-cos((z0-lambda)*pi2/(0.5-lambda)))
!       Mode 3:         =0 if |z|<lambda, otherwise E*(1-cos((z0-lambda)*pi/(0.5-lambda)))*0.5D0
!       Mode 4:         =E if |z|<lambda, otherwise E*(0.5D0-Z0)/(0.5d0-lambda)
!       Mode 5:         =0 if |z|<lambda, otherwise E*(1.D0-(0.5D0-Z0)/(0.5d0-lambda))
!       Mode 6,7:       Variants of forms from Lozovoi et al, J. Chem Phys, 115, 1661, 2001
!       Mode 8:         =E
!       Mode 9:         =0 if |z|<lambda, otherwise E
!       Mode 10:        =E if |z|<lambda, otherwise E*((1-cos((z0-lambda)*pit/cut))*0.5D0)
!       Mode 0:         Triangular Ramp, analytic form everywhere
!
!       Modes 20 added, use limited Fourier series
!       Modes 40 added, use full Fourier coefficients form not analytic form
!
!       1-D only; generalization is simple
!
        subroutine tester
!
!       This subroutine is just here for testing purposes
!       Comment out the first line, uncomment the program line
!        program tester
!       Not used
        implicit real*8 (a-h,o-z)
        real*8 lambda
        lambda=0.35
10      write(6,2)'Enter mode, lambda'
2       format(a,'_',$)
        read(5,*,err=999,end=999)mode,lambda
        ifourier=0
        do iz=0,10001,100
                z=iz/10000.d0
                b=backget(z,mode,lambda,ifourier)
                write(6,1)z,b
1               format(2F12.6)
        enddo
        goto 10
        return
999     continue
        end
!
        real*8 function backstep(z,lambda)
!       Step function
        implicit real*8 (a-h,o-z)
        real*8 lambda
        pi=acos(-1.0D0)
        z0=z
        if( Z0 .gt. 0.5D0) z0=1.D0-z0
        Z0=abs(z0)
        if(z0 .lt. lambda)then
                backstep=0D0
        else
                backstep=1.0D0
        endif
        return
        end
        real*8 function backstep2(z,lambda)
!       Step, with a soft edge
        implicit real*8 (a-h,o-z)
        real*8 lambda
        pit=acos(-1.0D0)*1.5D0
        cut=0.5-lambda
        top=0.5-cut/3.D0
        z0=z
        if( Z0 .gt. 0.5D0) z0=1.D0-z0
        Z0=abs(z0)
        if(z0 .le. lambda)then
                backstep2=0
        else if(z0 .ge. top)then
                backstep2=1.D0 
        else
                backstep2=(1-cos((z0-lambda)*pit/cut))*0.5D0
        endif
        return
        end

        real*8 function backldm3(z,lambda)
!       Barrier function
        implicit real*8 (a-h,o-z)
        real*8 lambda
        pi=acos(-1.0D0)
        z0=z
        if( Z0 .gt. 0.5D0) z0=1.D0-z0
        Z0=abs(z0)
        if(z0 .lt. lambda)then
                backldm3=0
        else
                backldm3=(1-cos((z0-lambda)*pi/(0.5-lambda)))*0.5D0
        endif
        return
        end
        real*8 function backldm(z,lambda)
!       Different type of barrier function
        implicit real*8 (a-h,o-z)
        real*8 lambda
        pi2=acos(-1.0D0)*0.5D0
        z0=z
        if( Z0 .gt. 0.5D0) z0=1.d0-z0
        Z0=abs(z0)
        if(z0 .lt. lambda)then
                backldm=1.D0
        else
                backldm=cos((z0-lambda)*pi2/(0.5-lambda))
        endif
        return
        end

        real*8 function backldm2(z,lambda)
!       Different type of barrier function
        implicit real*8 (a-h,o-z)
        real*8 lambda
        pi2=acos(-1.0D0)*0.5D0
        z0=z
        if( Z0 .gt. 0.5D0) z0=1.D0-z0
        Z0=abs(z0)
        if(z0 .lt. lambda)then
                backldm2=0.D0
        else
                backldm2=1.D0-cos((z0-lambda)*pi2/(0.5-lambda))
        endif
        return
        end

        real*8 function backramp0(z,lambda)
!       Zig-zag potential
        implicit real*8 (a-h,o-z)
        real*8 lambda
        z0=z
        if( Z0 .gt. 0.5D0) z0=1.D0-z0
        Z0=abs(z0)
        backramp0=(0.5D0-Z0)/0.5d0-0.5
        return
        end

        real*8 function backramp0z(z,lambda)
!       Zig-zag potential
        implicit real*8 (a-h,o-z)
        real*8 lambda
        z0=z
        if( Z0 .gt. 0.5D0) z0=1.D0-z0
        Z0=abs(z0)
        backramp0z=(0.5D0-Z0)/0.5d0-1.0
        return
        end

        real*8 function backramp(z,lambda)
!       Unity with a linear ramp
        implicit real*8 (a-h,o-z)
        real*8 lambda
        z0=z
        if( Z0 .gt. 0.5D0) z0=1.D0-z0
        Z0=abs(z0)
        if(z0 .lt. lambda)then
                backramp=1.D0
        else
                backramp=(0.5D0-Z0)/(0.5d0-lambda)
        endif
        return
        end
!
        real*8 function backramp2(z,lambda)
!       Zero with a linear ramp
        implicit real*8 (a-h,o-z)
        real*8 lambda
        z0=z
        if( Z0 .gt. 0.5D0) z0=1.D0-z0
        Z0=abs(z0)
        if(z0 .lt. lambda)then
                backramp2=0.D0
        else
                backramp2=1.D0-(0.5D0-Z0)/(0.5d0-lambda)
        endif
        return
        end

        real*8 function backlozovoi2(z,lzz)
!       Generate background potential following Lozovoi et al
        implicit real*8 (a-h,o-z)
        real*8 lambda,LZ,LZZ
        parameter (Lambda=0.5D0)
        z0=0
        Lz=LZZ*2.D0
        base=-(z0*z0-LZ*lambda+LZ*LZ*0.25)/LZ
!       Could use sqrt(2)*(exp(-z*z)-exp(-0.25))
!       Shift to between -0.5 & 0.5
        Z0=Z
        if(Z0 .gt. 0.5D0) Z0=1.D0-z0
!       Form
        if(z0 .le. -LZ*0.5D0) then
                backcor = lambda+Z0
        else if(z0 .ge. LZ*0.5D0) then
                backcor = (lambda-Z0)
        else
                backcor = -(z0*z0-LZ*lambda+LZ*LZ*0.25)/LZ
        endif
        backlozovoi2=backcor
        return
        end

	real*8 function backlozovoi(z,lzz)
!	Generate background potential following Lozovoi et al
!       Normalized to a maxium of 1
	implicit real*8 (a-h,o-z)
	real*8 lambda,LZ,lzz
	parameter (lambda=0.5D0)
        LZ=LZZ*2.0D0
        z0=0
        base=-(z0*z0-LZ*lambda+LZ*LZ*0.25)/LZ
!       Could use sqrt(2)*(exp(-z*z)-exp(-0.25))
!	Shift to between -0.5 & 0.5
	Z0=Z
	if(Z0 .gt. 0.5D0) Z0=1.d0-Z0
!	Form
	if(z0 .le. -LZ*0.5D0) then
		backcor = lambda+Z0
	else if(z0 .ge. LZ*0.5D0) then
		backcor = (lambda-Z0)/base
	else
		backcor = -(z0*z0-LZ*lambda+LZ*LZ*0.25)/LZ/base
	endif
        backlozovoi=backcor
	return
	end
!
        real*8 function backget(z,modein,lambda,ifourier)
!       Returns value depending upon mode
        use efeld, only : kefeld, fefeldpw,refeld,jefeld
        use symetr
        implicit real*8 (a-h,o-z)
        real*8 lambda
               mode=modein
               if(mode .ge. 40)then
                       if(ifourier .eq. 1)then
                          pi2=2.D0*acos(-1.d0)
                          b=fefeldpw(1)
                          const=z*pi2
                          do J=2,kefeld
                            jj=jefeld(j)
                            L=kzz(3,jj)
                            b=b+fefeldpw(J)*cos(const*L)
                          enddo
                          backget=b/refeld
                          return
                        else
                          mode=mode-40
                        endif
                else if(mode .ge. 20)then
                        mode=mode-20
                endif
                if(mode .eq. 1)then
                        backget=backldm(z,lambda)
                else if(mode .eq. 2)then
                        backget=backldm2(z,lambda)
                else if(mode .eq. 3)then
                        backget=backldm3(z,lambda)
                else if(mode .eq. 4)then
                        backget=backramp(z,lambda)
                else if(mode .eq. 5)then
                        backget=backramp2(z,lambda)
                else if(mode .eq. 6)then
                        backget=backlozovoi(z,lambda)
                else if(mode .eq. 7)then
                        backget=backlozovoi2(z,lambda)
                else if(mode .eq. 8)then
                        backget=1.D0
                else if(mode .eq. 9)then
                        backget=backstep(z,lambda)
                else if(mode .eq. 10)then
                        backget=backstep2(z,lambda)
                else if(mode .eq. 0)then
                        backget=backramp0(z,lambda)
                else if(mode .eq. 11)then
                        backget=backramp0(z,lambda)
                else if(mode .eq. 12)then
                        backget=backramp0z(z,lambda)
                else
                        CALL OUTERR('LAPW0','Unknown E-field mode')
                        stop 'Unknown E-field mode'
                endif
        return
        end

        subroutine makeback(fefeldpw,lambda,N,wt,mode,myid)
!       Create array of values using an FFT
!       pw(*), for PW's
!       N is the total number of K-vectors (largest is N-1)
        use symetr
        use efeld, only : kebig, jefeld, refeld,iefeld
        implicit real*8 (a-h,o-z)
        real*8 lambda, fefeldpw(*)
        complex*16,allocatable ::  W(:)
        real*8,allocatable     ::  DWORK(:)
        logical debug
        debug=.false.
        if(iefeld .lt. 0)debug=.true.
        fefeldpw(1:N+1)=0
!       Analytical case
        if((mode .eq. 0).or.(mode.eq.40).or.(mode.eq.20))then
            const=wt*2.D0/(acos(-1.D0)**2)
            if(debug)write(6,*)'Efield Fourier Coefficients'
            if(debug)write(6,1)0,1,0.0D0
1           format('KZ=',I4,' Degeneracy ',i2,' Value ',F12.6)
            do j=2,N
               jj=jefeld(j)
               L=kzz(3,jj)
               if(2*(L/2).ne.L)fefeldpw(j)=const*inst(jj)/dble(L*L)
               if(debug)write(6,1)L,inst(jj),fefeldpw(j)
            enddo
            if(debug)then
                write(6,*)'Efield values along z'
                N2=2*kebig+2
                t=1.D0/dble(N2)
                do J=1,N2
                z=dble(J-1)*t
                t1=backget(z,mode,lambda,ifourier)*refeld
                write(6,80)j,z,t1
80              format(i4,3F12.6)
                enddo
            endif
            return
        endif
!
!       Size for just a finite Fourier series
        N2=2*kebig+2
!
!       Use a larger size to approach an infinite Fourier series
        if(mode .lt. 19)N2=max(4*(kebig+1),2048*2)
        write(6,100)N2
        if(myid .eq. 0)write(21,100)N2
100     format(':EFIELD       Fourier Size ',I4)
!
        allocate (W(N2), DWORK(4*N2+20))
        t=1.D0/dble(N2)
!       Establish in w an array of values
!       Ifourier=0 so we don't sum terms if requested for mode>40
        ifourier=0
        if(debug)write(6,*)'Efield values along z'
        do J=1,N2
                z=dble(J-1)*t
                w(j)=backget(z,mode,lambda,ifourier)
                if(debug)write(6,80)j,z,dble(w(j))*refeld
        enddo
!       Do a 1-D fft
        CALL CFFTI(N2,DWORK)
        CALL CFFTF(N2,W,DWORK)
!
!       Copy, taking care of scaling
        w2=wt/dble(N2)
        fefeldpw(1)=dble(w(1))*w2
        if(debug)write(6,1)0,1,fefeldpw(1)
!
!       Let inst handle the degeneracy for us
        if(debug)write(6,*)'Efield Fourier Coefficients'
        do J=2,N
                jj=jefeld(j)
                L=kzz(3,jj)+1
                if(L .lt. 1)L=N2+L
                fefeldpw(j)=DBLE(w(L))*w2*inst(jj)
                if(debug)write(6,1)kzz(3,jj),inst(jj),fefeldpw(j)
        enddo
        deallocate (W,DWORK)
        end
!
        real*8 function backgetder(z)
!       Returns gradient depending upon mode
        use efeld
        use symetr
        use struct
        implicit real*8 (a-h,o-z)
!       Do linear case exactly
        imode=abs(iefeld)/1000
        if((imode .eq. 0).or.(imode.eq.11) .or. (imode.eq.12))then
          if((abs(z) .lt. 1D-5) .or. (abs(z-0.5) .lt. 1D-5))then
                backgetder=0
          else
                backgetder=-refeld*2.D0/a(3)
                if(z .gt. 0.5)backgetder=-backgetder
          endif
          return
        else if (imode .eq. 8)then
                backgetder=0
                return
        endif
        pi2=2.D0*acos(-1.d0)
        b=0
        const=z*pi2
        do J=2,kefeld
            jj=jefeld(j)
            L=kzz(3,jj)
            b=b-fefeldpw(J)*sin(const*L)*L
        enddo
        backgetder=b*pi2/a(3)
!        write(*,*)a(3)
        return
	end

