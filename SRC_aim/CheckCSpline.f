      subroutine CheckCSpline(ErrorSum,ErrorAbs,weit)
!
!     This code does an error analysis of the Bader surface in case.surface
!     The method is:
!       a) Fit cubic splines to the theta, phi values
!       b) From these, estimate two vectors in the surface numerically
!       c) Generate the normal to the surface from this
!       d) Calculate the dot product with the gradient, and print data to the
!          file (currently) Surface_Errors
!       e) Accumulate the error metrics
!
!               E1 = Sum (n.g)^2*wt/Sum wt
!               E2 = Sum |n.g|*wt/Sum wt
!               E3 = Sum |n.g|*wt*|g|/Sum wt*|g|
!               E4 = Sum |n.g|*wt*rho/Sum wt*rho
!
!         where n is the normal vector, g the gradient on the Bader surface, wt
!         a quadrature weight and rho the density. Most useful are E3 and E4
!         since these seem to be reasonable metrics. Values converted into
!         degrees in the form E_reported = acos(E)*180/pi are printed
!
!       The code also gives a crude estimate of the error, more to caution the
!       user than as a true number
!
      use sphe
      use zsphe
      implicit real*8 (a-h,o-z)
!
      real*8,allocatable  :: rsurface(:,:),rtheta(:),rphi(:),bcoef(:,:),WT(:,:)
      real*8,allocatable  :: TX(:),TY(:),work(:),DTH(:,:),DPH(:,:)
      real*8 D1(3),D2(3),D3(3)
      real*8 grho(3),rho,v(3),hrho(3,3),VC(2)
      real*8 ThetaInt(2),PhiInt(2),RT(2),RP(2),NError
      integer IC(2)
      logical weit
!
!       Initialize things
        ErrorSum=0.D0
        ErrorAbs=0.D0
        GSsum=0.0D0
        GAsum=0.0D0
        GGTOT=0.0D0
        GTOT=0.0D0
        RSUM=0.0D0
        RErr=0.0D0
!       Read the data
        rewind(unit=21)
        read(21,*) index,shift
        read(21,*) nth,themin,themax
        read(21,*) nph,phimin,phimax
!       Fix numerical accuracy errors in angles
        tol=1d-5
        tpi=2.D0*acos(-1.D0)
        do i=-8,8
          if(i.eq.0)then
                if(abs(themin) .lt. tol)themin=0.D0
                if(abs(themax) .lt. tol)themax=0.D0
                if(abs(phimin) .lt. tol)phimin=0.D0
                if(abs(phimax) .lt. tol)phimax=0.D0
          else
                t=tpi/dble(i)
                if(abs(themin-t) .lt. tol)themin=t
                if(abs(themax-t) .lt. tol)themax=t
                if(abs(phimin-t) .lt. tol)phimin=t
                if(abs(phimax-t) .lt. tol)phimax=t
          endif
        enddo
!       Workspace
        allocate (rphi(nph),rtheta(nth),rsurface(nph,nth))
        allocate (DTH(nph,nth), DPH(nph,nth),WT(nph,nth))
        NWK=(max(nph,nth)+1)*2
        allocate ( work(NWK) )
!       Read the data, storing values
        do i=1,nth
          do j=1,nph
            if(weit)then
                read(21,*) theta,phi,rmax,WT(j,i)
            else
                read(21,*) theta,phi,rmax
                WT(j,i)=1.D0
            endif
            rsurface(j,i)=rmax
            rtheta(i)=theta
            rphi(j)=phi
          end do
        end do
!
	IC(1)=0
	IC(2)=0
	switch=0
	INCFD=1
!       Setup Phi interpolation table
	do j=1,nth
	   call DPCHIC(IC,VC,SWITCH,NPH,rphi,rsurface(1,j),DPH(1,j), &
           INCFD,Work,NWK,IERR)
        enddo
!       Setup Theta interpolation table
!       Note: it would be reasonable to use cos(theta), but this does not seem to matter
        INCFD=nph
        do j=1,nph
	   call DPCHIC(IC,VC,SWITCH,NTH,rtheta,rsurface(j,1),DTH(j,1), &
           INCFD,Work,NWK,IERR)
        enddo
!
!       Error analysis, skipping end points
        NError=0
!       Output, should be changed
!        open(unit=77,file='Surface_Errors')
        do j=2,nth-1 ! nth-1,2,-1
             theta=rtheta(j)
             ct=cos(theta)
             st=sin(theta)
!            Two locations, on either side in theta
             TH1=Theta+1D-4
             TH2=Theta-1D-4
             CT1=cos(TH1)
             CT2=cos(TH2)
             ST1=sin(TH1)
             ST2=sin(TH2)
             ThetaInt(1)=TH1 !CT1
             ThetaInt(2)=TH2 !CT2
!
             do k=2,nph-1
                phi=rphi(k)
                cf=cos(phi)
                sf=sin(phi)
                rmax=rsurface(k,j)
!
!               Get gradient at the point
                v(1)=pos(1,index)+rmax*st*cf
                v(2)=pos(2,index)+rmax*st*sf
                v(3)=pos(3,index)+rmax*ct
                call vgh_rho(v,rho,grho,hrho,.true.,.true.,.false.,&
                rr,iat,ipos,ist,ibest,ipbest)
!
                PH1=PHI+1D-4
                PH2=PHI-1D-4
                PhiInt(1)=PH1
                PhiInt(2)=PH2
!	Get values along phi
	call DPCHFE(NPH,Rphi,rsurface(1,j),DPH(1,j),1,.true.,2,PhiInt,RP,IERR)
!       First vector in surface
                D1(1)=RP(1)*st*cos(ph1)-RP(2)*st*cos(ph2)
                D1(2)=RP(1)*st*sin(ph1)-RP(2)*st*sin(ph2)
                D1(3)=RP(1)*CT-RP(2)*CT

!       Get values along theta
	call DPCHFE(NTH,Rtheta,rsurface(k,1),DTH(k,1),NPH,.true.,2,ThetaInt,RT,IERR)
!       Second vector
                D2(1)=RT(1)*ST1*cf-RT(2)*ST2*cf
                D2(2)=RT(1)*ST1*sf-RT(2)*ST2*sf
                D2(3)=RT(1)*CT1-RT(2)*CT2
!
!       Rescale so they have unit length (for safety)
                DD1=vnorm(D1)
                DD2=vnorm(D2)
                G=vnorm(grho)
                do jj=1,3
                        D1(jj)=D1(jj)/DD1
                        D2(jj)=D2(jj)/DD2
                        grho(jj)=grho(jj)/G
                enddo
!       Cross product
                D3(1)=D1(2)*D2(3)-D1(3)*D2(2)
                D3(2)=D1(3)*D2(1)-D1(1)*D2(3)
                D3(3)=D1(1)*D2(2)-D1(2)*D2(1)
!       Renormalize cross product
                DD1=vnorm(D3)
                do jj=1,3
                        D3(jj)=D3(jj)/DD1
                enddo
!       Dot with gradient, and accumulate error metrics
                DOTP = D3(1)*Grho(1)+D3(2)*grho(2)+D3(3)*grho(3)
                GGSUM=GGSUM+G*G*WT(k,j)
                GSUM=GSUM+G*wt(k,j)
                RSUM=RSUM+rho*wt(k,j)
                GSsum=GSsum+DOTP*DOTP*G*G*wt(k,j)
                GAsum=GAsum+abs(dotp)*G*wt(k,j)
                Rerr=Rerr+abs(dotp)*rho*wt(k,j)
                ErrorSum=ErrorSum+DOTP*DOTP*wt(k,j)
                ErrorAbs=ErrorAbs+abs(dotp)*wt(k,j)
                NError=NError+wt(k,j)
!                write(*,77)theta,phi,DOTP,D3,grho
                write(77,77)theta,phi,rmax,acos(DOTP)*180.D0/acos(-1.D0),D3,grho,rho,G
77              format(3f10.6,'  Angle ',F10.5,'  Normal',3f9.5,'  Grad',3F9.5,' Rho ', &
                D12.5,' |G| ',D12.5)
           enddo
        enddo
      GSsum=sqrt(GSsum/GGSUM)
      GAsum=GAsum/Gsum
      Rerr = Rerr/Rsum
      GSsum=acos(GSsum)*180.D0/acos(-1.D0)-90.D0
      GAsum=acos(GAsum)*180.D0/acos(-1.D0)-90.D0
      Rerr =acos(Rerr )*180.D0/acos(-1.D0)-90.D0
      ErrorSum=sqrt(ErrorSum/NError)
      ErrorSum=acos(ErrorSum)*180.D0/acos(-1.D0)-90.D0
      ErrorAbs=ErrorAbs/NError
      ErrorAbs=acos(ErrorAbs)*180.D0/acos(-1.D0)-90.D0
      write(6,303) 'Error analysis, values in degrees'
      write(6,303) 'RMS Angular Error       ',abs(ErrorSum),' Sum of moduli ',abs(ErrorAbs)
      write(77,303)'RMS Angular Error       ',abs(ErrorSum),' Sum of moduli ',abs(ErrorAbs)
      write(6,303) 'Gradient weighted Error ',abs(Gasum),   ' Rho weighted  ',abs(Rerr)
      write(77,303)'Gradient weighted Error ',abs(Gasum),   ' Rho weighted  ',abs(Rerr)
!
!     Roughly correct, based upon tests with TiO2
      errest=(abs(Gasum)*0.003+abs(Rerr)*0.0045)*0.5
      write(6,303) 'Estimated Charge Error  ',errest, ' Electrons (only a rough estimate!)'
      close(unit=77)
!     Tidy up
      deallocate (rphi,rtheta,rsurface)
      deallocate (DTH, DPH,WT, work)
!
      return
 303  format(a,F13.8,a,f12.6)
      end
	
