      real*8 function Energy(translat,distances,br,nuse,pred,N2)
!     Pairwise Hessian initializer
!     Version 1.0, Feb 2006
!     J. Rondinelli, B. Deng and L. D. Marks
!
      implicit real*8 (a-h,o-z)
      include 'params.inc'
      dimension br(3,3),nuse(natmax)
      integer translat(4,neigmax,48*natmax)
      dimension pall(3,48*natmax),pred(3,natmax)
      dimension pwork(3,natmax),pose(3,48*natmax)
      dimension distances(neigmax,natmax*48)
!      save debug,icall
!      logical debug
!      data debug/.true./,icall/0/
!      icall=icall+1
!      write(78,*)icall
!
!       Evaluate a simple pair-wise energy
        E=0.D0
!       Get equivalent expansions
        call genequivs(pred,pose)
!       First atom of expanded set
        II=0
        do j=1,NAT
!        Loop over equivalents for each site
         EAT=0.D0
         DO NE=1,nineq(J)
          II=II+1
          do K=1,NAT
             pwork(1:3,K)=pred(1:3,K)
          enddo
          pwork(1:3,J)=pose(1:3,II)
          call expandset(pwork,pall)
           M=nuse(j)
           dmin=distances(1,j)
!          First atom, convert to cartesians
           x11=pwork(1,j)*BR(1,1)+pwork(2,j)*BR(2,1)+pwork(3,j)*BR(3,1)
           x12=pwork(1,j)*BR(1,2)+pwork(2,j)*BR(2,2)+pwork(3,j)*BR(3,2)
           x13=pwork(1,j)*BR(1,3)+pwork(2,j)*BR(2,3)+pwork(3,j)*BR(3,3)
!          Loop over neighbors
           do n=1,M
!           i1=translat(1,n,j)
!           i2=translat(2,n,j)
!           i3=translat(3,n,j)
!           k =translat(4,n,j)
!           x21=(pall(1,k)-i1)*BR(1,1)+(pall(2,k)-i2)*BR(2,1)+(pall(3,k)-i3)*BR(3,1)
!           x22=(pall(1,k)-i1)*BR(1,2)+(pall(2,k)-i2)*BR(2,2)+(pall(3,k)-i3)*BR(3,2)
!           x23=(pall(1,k)-i1)*BR(1,3)+(pall(2,k)-i2)*BR(2,3)+(pall(3,k)-i3)*BR(3,3)
!           p11=x21-x11
!           p12=x22-x12
!           p13=x23-x13
              k =translat(4,n,j)
              t1=pall(1,k)-translat(1,n,j)
              t2=pall(2,k)-translat(2,n,j)
              t3=pall(3,k)-translat(3,n,j)
!          Second atom, convert to cartesians
!          Differences
           p11=t1*BR(1,1)+t2*BR(2,1)+t3*BR(3,1)-x11
           p12=t1*BR(1,2)+t2*BR(2,2)+t3*BR(3,2)-x12
           p13=t1*BR(1,3)+t2*BR(2,3)+t3*BR(3,3)-x13
!          Add energy term
            if(Simple)then
                p11=p11-atrans(1,n,j)
                p12=p12-atrans(2,n,j)
                p13=p13-atrans(3,n,j)
                EAT=EAT+(P11*P11+P12*P12+P13*P13)*atrans(4,n,j)
           else
!              Distance
               d0=distances(n,j)
               d=sqrt(p11*p11+p12*p12+p13*p13)/d0-1.D0
               EAT=EAT+d*d*atrans(4,n,j)
           endif
           enddo
         enddo
         E=E+EAT*WMULT(J)
        enddo
        Energy=E
        return
        end
!
      real*8 function Energyb(translat,distances,br,nuse,pred,N2)
!     Pairwise Hessian initializer
!     Version 1.0, Feb 2006
!     J. Rondinelli, B. Deng and L. D. Marks
!     Version that fixes atoms
!
      implicit real*8 (a-h,o-z)
      include 'params.inc'
      dimension br(3,3),nuse(natmax)
      integer translat(4,neigmax,48*natmax)
      dimension pall(3,48*natmax),pred(3,natmax)
      dimension pwork(3,natmax),pose(3,48*natmax)
      dimension distances(neigmax,natmax*48)
!      save debug,icall
!      logical debug
!      data debug/.true./,icall/0/
!      icall=icall+1
!      write(78,*)icall
!
!       Evaluate a simple pair-wise energy
        E=0.D0
!       Get equivalent expansions
        call expandset(pred,pall)
!       First atom of expanded set
        II=0
        do j=1,NAT
!        Loop over equivalents for each site
           EAT=0.D0
           M=nuse(j)
           dmin=distances(1,j)
!          First atom, convert to cartesians
           x11=pred(1,j)*BR(1,1)+pred(2,j)*BR(2,1)+pred(3,j)*BR(3,1)
           x12=pred(1,j)*BR(1,2)+pred(2,j)*BR(2,2)+pred(3,j)*BR(3,2)
           x13=pred(1,j)*BR(1,3)+pred(2,j)*BR(2,3)+pred(3,j)*BR(3,3)
!          Loop over neighbors
           do n=1,M
              k =translat(4,n,j)
              t1=pall(1,k)-translat(1,n,j)
              t2=pall(2,k)-translat(2,n,j)
              t3=pall(3,k)-translat(3,n,j)
!             Second atom, convert to cartesians
!             Differences
              p11=t1*BR(1,1)+t2*BR(2,1)+t3*BR(3,1)-x11
              p12=t1*BR(1,2)+t2*BR(2,2)+t3*BR(3,2)-x12
              p13=t1*BR(1,3)+t2*BR(2,3)+t3*BR(3,3)-x13
!             Add energy term
              if(Simple)then
                p11=p11-atrans(1,n,j)
                p12=p12-atrans(2,n,j)
                p13=p13-atrans(3,n,j)
                EAT=EAT+(P11*P11+P12*P12+P13*P13)*atrans(4,n,j)
              else
!              Distance
               d0=distances(n,j)
               d=sqrt(p11*p11+p12*p12+p13*p13)/d0-1.D0
               EAT=EAT+d*d*atrans(4,n,j)
              endif
           enddo
         E=E+EAT*WMULT(J)
        enddo
        Energyb=E
        return
        end
!
