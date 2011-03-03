      subroutine integr(nkk,kzz,rhok,nat,mult,jspin,jri,rnot, &
      dx,r,clmnew,iff1,iff2,iff3,qpw,qel,cfft,ust,margn)
!.....integrate charge density
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      character*79 margn
      character*8 atype
      character*2 atype1
      complex*16 cfft(IFF1*IFF2*IFF3),ust((IFF1+1)*(IFF2+1)*(IFF3+1))
      complex*16 rhok(nkk,2),cvalue
      dimension kzz(3,nkk),jri(nat),rnot(nat),dx(nat)
      dimension r(nrad),clmnew(NRAD,NCOM,NAT,2)
      dimension qpw(2),qel(nat,2)
      integer mult(nat)  
!
!     calc interstitial charge
      do 10 ispin=1,jspin
         CALL SETFFX(NKK,IFF1,IFF2,IFF3,RHOK(1,ispin),CFFT,kzz)
         CALL REAN3(NKK,KZZ,CVALUE,IFF1,IFF2,IFF3,CFFT,UST,UST0)
         qpw(ispin)=CVALUE*UST0
         if(ispin.eq.1) atype='SPIN-UP '
         if(ispin.eq.1) atype1='UP'
         if(ispin.eq.2) atype='SPIN-DN '
         if(ispin.eq.2) atype1='DN'
         if(jspin.eq.1) atype='TOTAL   '
         if(jspin.eq.1) atype1='TO'
         write(6,2040) margn,atype1,atype,qpw(ispin)
         write(21,2040)margn,atype1,atype,qpw(ispin)
 2040 FORMAT(':',a1,a2,'  :',1X,A8,'INTERSTITIAL CHARGE=',F11.7)              
 10   continue
      if(jspin.eq.2) then
          atype='TOTAL   '
          atype1='TO'                                
          write(6,2040) margn,atype1,atype,qpw(1)+qpw(2)
          write(21,2040)margn,atype1,atype,qpw(1)+qpw(2)
            qmtot=qpw(1)-qpw(2)
      end if
!
!.....INTEGRATE CHARGES IN SPHERES                                      
!     
         atype='SPIN-UP '
         atype1='UP'
         if(jspin.eq.1) atype='TOTAL   '
         if(jspin.eq.1) atype1='TO'
      DO 570 JATOM=1,NAT                                                
         DO 575 J=1,JRI(JATOM)                                          
575      R(J)=RNOT(JATOM)*EXP( DX(JATOM)*(J-1) )                        
      CALL CHARGE(R,DX(JATOM),CLMNEW(1,1,JATOM,1),1,JRI(JATOM), &
         QEL(jatom,1))    
      WRITE(6,2041)  margn,atype1,jatom,atype,JATOM,QEL(jatom,1)                     
570   WRITE(21,2041) margn,atype1,jatom,atype,JATOM,QEL(jatom,1)
 2041 FORMAT(':',a1,a2,i3.3,':',1X, &
             A8,'CHARGE IN SPHERE',I3,' = ',2F10.7)            
!                                                                       
      if(jspin.eq.2) then
!                                                                       
      DO 670 JATOM=1,NAT                                                
         DO 675 J=1,JRI(JATOM)                                          
675      R(J)=RNOT(JATOM)*EXP( DX(JATOM)*(J-1) )                        
      CALL CHARGE(R,DX(JATOM),CLMNEW(1,1,JATOM,2),1,JRI(JATOM), &
         QEL(jatom,2))   
         atype='SPIN-DN '
         atype1='DN'
      WRITE(6,2041)  margn,atype1,jatom,atype,JATOM,QEL(jatom,2)
670   WRITE(21,2041) margn,atype1,jatom,atype,JATOM,QEL(jatom,2)
      DO 671 JATOM=1,NAT                                                
         atype='TOTAL   '
         atype1='TO'
      WRITE(6,2041) margn,atype1,jatom,atype,JATOM, &
      QEL(jatom,1)+QEL(jatom,2)
 671  WRITE(21,2041)margn,atype1,jatom,atype,JATOM, &
      QEL(jatom,1)+QEL(jatom,2)
!
      if(margn(1:1).eq.'C') then
        write(6,2042) qmtot
        write(21,2042) qmtot
        do jatom=1,nat
        qmtot=qmtot+(QEL(jatom,1)-QEL(jatom,2))*mult(jatom)
        write(6,2043) jatom,jatom,QEL(jatom,1)-QEL(jatom,2)
        write(21,2043) jatom,jatom,QEL(jatom,1)-QEL(jatom,2)
        enddo
        write(6,2044) qmtot
        write(21,2044) qmtot
      endif
 2042 FORMAT(/,'       MAGNETIC MOMENTS OF MIXED CHARGE DENSITY',/ &
       ':MMINT:',1X, &
       'MAGNETIC MOMENT IN INTERSTITIAL = ',2F10.5)            
 2043 FORMAT(':MMI',i3.3,':',1X, &
       'MAGNETIC MOMENT IN SPHERE ',i3,'    = ',2F10.5)            
 2044 FORMAT(':MMTOT:',1X, &
       'TOTAL MAGNETIC MOMENT IN CELL   = ',2F10.5)            
!
      end if
      return
      end
