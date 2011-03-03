      subroutine distan(nkk,kzz,roknew,rokold,rhok, &
           nat,jspin,jri,rnot,dx,mult,r,rho, &
           clmnew,clmold,iff1,iff2,iff3,qpw,qel,cfft,ust,qtot)
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      character*8 atype
      character*2 atype1
      complex*16 cfft(IFF1*IFF2*IFF3),ust((IFF1+1)*(IFF2+1)*(IFF3+1))
      complex*16 roknew(nkk,2),rokold(nkk,2),rhok(nkk,2),cvalue
      dimension kzz(3,nkk),jri(nat),rnot(nat),dx(nat),mult(nat)
      dimension r(nrad),clmnew(NRAD,NCOM,NAT,2),rho(NRAD,2)
      dimension qpw(2),qel(nat,2),clmold(NRAD,NCOM,NAT,2)
!
      qtot=0.d0
!      do 1 ispin=1,jspin
!         do 2 k=1,nkk
! 2          rhok(k,ispin)=(roknew(k,ispin)-rokold(k,ispin))
!c........the proper way would be: take difference as above,
!c........fft to real space, square, fft to rec.space
!         CALL SETFFX(NKK,IFF1,IFF2,IFF3,RHOK(1,ispin),CFFT,kzz)
!         CALL REAN3(NKK,KZZ,CVALUE,IFF1,IFF2,IFF3,CFFT,UST,UST0)
!         qpw(ispin)=sqrt(CVALUE*UST0)
!         qtot=qtot+qpw(ispin)
!         if(ispin.eq.1) atype='SPIN-UP '
!         if(ispin.eq.2) atype='SPIN-DN '
!         if(jspin.eq.1) atype='TOTAL   '
!         write(6,2041) atype,qpw(ispin)
!         write(21,2041)atype,qpw(ispin)
! 1    continue
!      if(jspin.eq.2) then
!          atype='TOTAL   '
!          write(6,2041) atype,qpw(1)+qpw(2)
!          write(21,2041)atype,qpw(1)+qpw(2)
!      end if
      denom=0.d0
      do 10 ispin=1,jspin
         do 10 jatom=1,nat
         do 20 j=1,jri(jatom)
            R(J)=RNOT(JATOM)*EXP( DX(JATOM)*(J-1) )                        
 20         rho(j,ispin)=(clmnew(j,1,jatom,ispin)- &
                          clmold(j,1,jatom,ispin))**2
         CALL CHARGE(R,DX(JATOM),rho(1,ispin),1,JRI(JATOM), &
         QEL(jatom,ispin)) 
         qel(jatom,ispin)=sqrt(qel(jatom,ispin))  
         qtot=qtot+qel(jatom,ispin)*mult(jatom)
         denom=denom+mult(jatom)
         if(ispin.eq.1) atype='SPIN-UP '
         if(ispin.eq.1) atype1='UP'
         if(ispin.eq.2) atype='SPIN-DN '
         if(ispin.eq.2) atype1='DN'
         if(jspin.eq.1) atype='TOTAL   '
         if(jspin.eq.1) atype1='TO'
         WRITE(6,2040)  atype1,JATOM,atype,JATOM,QEL(jatom,ispin)
 10      WRITE(21,2040) atype1,JATOM,atype,JATOM,QEL(jatom,ispin)
 2040 FORMAT(':D',a2,i3.3,':',1X, &
      A8,'DIFFERENCE CHARGE**2 IN SPHERE',I3,' = ',F10.7)            
 2041 FORMAT(3X,A8,'DIFFERENCE INTERSTITIAL CHARGE**2= ',2F10.7)              
 2042 FORMAT(/,':DIS  :',1x,' CHARGE DISTANCE ',f15.7)
      if(jspin.eq.2) then
         atype='TOTAL   '
         atype1='TO'
         DO 671 JATOM=1,NAT                                                
            WRITE(6,2040) atype1,JATOM,atype,JATOM, &
                          QEL(jatom,1)+QEL(jatom,2) 
 671     WRITE(21,2040)atype1,JATOM,atype,JATOM, &
                          QEL(jatom,1)+QEL(jatom,2)    
      end if       
      write(6,2042) qtot
      qtot=qtot/denom*jspin
      write(21,2042) qtot
      return
      end
