      subroutine parse(config,norb0,iex,zk,z,nuc,dval,TITRE, &
      TTIRE,aplot)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      character*1 config(20)
      character*1 aplot(30,2)                                     
      logical allplot
      COMMON DEN(30,2),DQ1(30),DFL(30),NQN(30),NQL(30),NK(30),NMAX(30),  &
      ZEL(30,2),NORB                                                    
      DIMENSION TITRE(30),TTIRE(9)                                        
!
         nqn(1)=1
         nk(1)=-1
         nqn(2)=2
         nk(2)=-1
         nqn(3)=2
         nk(3)= 1
         nqn(4)=2
         nk(4)=-2
         nqn(5)=3
         nk(5)=-1
         nqn(6)=3
         nk(6)= 1
         nqn(7)=3
         nk(7)=-2
         nqn(8)=3
         nk(8)= 2
         nqn(9)=3
         nk(9)=-3
         nqn(10)=4
         nk(10)=-1
         nqn(11)=4
         nk(11)= 1
         nqn(12)=4
         nk(12)=-2
         nqn(13)=4
         nk(13)= 2
         nqn(14)=4
         nk(14)=-3
         nqn(15)=5
         nk(15)=-1
         nqn(16)=5
         nk(16)= 1
         nqn(17)=5
         nk(17)=-2
         nqn(18)=4
         nk(18)= 3
         nqn(19)=4
         nk(19)=-4
         nqn(20)=5
         nk(20)= 2
         nqn(21)=5
         nk(21)=-3
         nqn(22)=6
         nk(22)=-1
         nqn(23)=6
         nk(23)= 1
         nqn(24)=6
         nk(24)=-2
!
      norb0=0
      i=1
      allplot=.false.
!.....start parsing config
 1    if(config(i).eq.'R') then
         norb0=24
      else if(config(i).eq.'X') then
         norb0=17
      else if(config(i).eq.'K') then
         norb0=12
      else if(config(i).eq.'A') then
         norb0=7
      else if(config(i).eq.'N') then
         norb0=4
      else if(config(i).eq.'H') then
         norb0=1
      else if(config(i).eq.'P') then
         allplot=.true.
      else if(ichar(config(i)).gt.47.and.ichar(config(i)).lt.58) then
         if(ichar(config(i+1)).gt.47.and.ichar(config(i+1)).lt.58) then
            norb=(ichar(config(i))-48)*10+ichar(config(i+1))-48+norb0
         else
            norb=ichar(config(i))-48+norb0
         endif
         do 3 i1=i+2,20
         if(ichar(config(i1)).gt.47.and.ichar(config(i1)).lt.58) then
!            if(ichar(config(i1+1)).gt.47.and.
!     *      ichar(config(i1+1)).lt.58) then
!cc               iex=(ichar(config(i1))-48)*10+ichar(config(i1+1))-48
!            else
!cc               iex=ichar(config(i1))-48
!            endif
            goto 2
         endif
 3       continue
      goto 2
      endif 
      i=i+1
      if(i.gt.18) stop 'invalid atomic configuration'
      goto 1
!.....end of parsing config
 2    continue 
      DO 24 I=1,norb0                                                   
      DO 24 ISPIN=1,2                                              
      zel(i,ispin)=iabs(nk(i))

      if (allplot) then 
         aplot(i,ispin)='P'
      else
         aplot(i,ispin)=' '
      endif

!      READ(5,1070) NQN(I),NK(I),ZEL(I,ISPIN),aplot(i,ispin)     
! DEN ENERGIE DE L ORBITALE EN UNITE ATOMIQUE ET NEGATIVE               
! NQN NOMBRE QUANTIQUE PRINCIPAL   NK NOMBRE QUANTIQUE KAPPA            
! ZEL OCCUPATION DE L ORBITALE                                          
 1070 FORMAT(I1,1X,i2,1x,f5.3,a1)                                   
! **********************************************************************
      ZK=ZK+ZEL(I,ISPIN)                                                
 18   DEN(I,ISPIN)=-Z*Z/(4.*NQN(I)*NQN(I))                              
 19   NQL(I)=IABS(NK(I))                                                
      IF (NK(I).LT.0) NQL(I)=NQL(I)-1                                   
      IF (NUC.GT.0) GO TO 21                                            
      DFL(I)=NK(I)*NK(I)                                                
      DFL(I)=SQRT(DFL(I)-DVAL)                                          
      GO TO 22                                                          
 21   DFL(I)=IABS(NK(I))                                                
 22   L=2*IABS(NK(I))                                                   
      NEL=ZEL(I,ISPIN)                                                  
      IF (NQL(I).LE.NQN(I).AND.NEL.LE.L.AND.NQN(I).GT.0.AND.NQL(I).LE.4) GO TO 23                                                      
 2002 FORMAT (' ERROR IN INPUT ',E15.8,I1,I2,F5.2)                      
      WRITE (6,2002) DEN(I,ISPIN),NQN(I),NQL(I),J,ZEL(I,ISPIN)          
      stop '598 in parse'                                                  
 23   J=NQL(I)+IABS(NK(I))                                              
      TITRE(I)=TTIRE(J)                                                 
 2005 FORMAT (7X,I1,A2,F16.3,1PE23.7,5x,a1)                                   
      WRITE (6,2005) NQN(I),TITRE(I),ZEL(I,ISPIN),DEN(I,ISPIN), &
       aplot(i,ispin)          
 24   continue
!
      return
      end
