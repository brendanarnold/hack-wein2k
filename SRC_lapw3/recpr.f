      SUBROUTINE RECPR(NKK,INDMAX,KXMAX,KYMAX,KZMAX,GMAX)           
      use kgrid 
      use reallocate 
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      LOGICAL KDELTA,ORTHO
!                                                                       
      CHARACTER*4 LATTIC                                
      CHARACTER*10 TITLE(8)                                             
      COMMON /CHAR/   TITLE,LATTIC                   
      COMMON /ORTH/   ORTHO                                            
      COMMON /GENER/ BR1(3,3),BR2(3,3)                                  
      COMMON /STRUK/ AA,BB,CC,alpha(3),PIA(3),VOL,ndif
      COMPLEX*16          TAUP
      COMMON /FOUR/ IM(3),NST,ISTM(3,48),TAUP(48)                       
      COMMON /SYM2/ IORD,IZ(3,3,48),ITAU,TAU(3,48)                      
      DIMENSION AM(3),M(3)                                 
      save l,test,nk1,ll
!C    L=16                                                              
      L=32                                                              
!      L=08                                                             
      TEST=1.d-05                                                       
      nkd=11111
      NK=nkd                                                          
      allocate (absk(nkd+1), kzz(3,nkd+1),inst(nkd))
      LL=2*L+1                                                          
!                                                                       
!..... READ IN SYMMETRY OPERATIONS AND NONPRIMITIVE TRANSLATIONS        
!                                                                       
      READ(20,100) IORD                                                 
      ITAU=0                                                            
      DO 15 I=1,IORD                                                    
      READ(20,101) ((IZ(I1,I2,I),I1=1,3),TAU(I2,I),I2=1,3)              
      DO 16 IND=1,3                                                     
   16 IF(TAU(IND,I).NE.0.0d0) GOTO 7777                                   
      GOTO 15                                                           
 7777 ITAU=ITAU+1                                                       
   15 CONTINUE                                                          
!                                      
      nfk=11111                                 
      allocate ( TAUK(NFK),KREC(3,NFK))                            
      ENTRY RECPR2(NKK,INDMAX,KXMAX,KYMAX,KZMAX,GMAX)                 
      NK1=NK+1                                                          
      J=1                                                               
      NK=0                                                              
      LL1=KXMAX+1                                                       
      LL2=KYMAX+1                                                       
      LL3=KZMAX+1                                                       
      ILL1=-ll1+2                                                       
      ILL2=-ll2+2                                                       
      ILL3=-ll3+2                                                       
!      
      DO 1 I1=ILL1,LL1                                                 
      M(1)=I1-1                                                         
      DO 1 I2=ILL2,LL2                                                   
      M(2)=I2-1                                                         
      DO 1 I3=ILL3,LL3                                                
      M(3)=I3-1                                                         
!                                                                       
      ABSM=0.0d0                                                          
      DO 22 IND1=1,3                                                    
      AM(IND1)=0.0d0                                                      
      DO 23 IND2=1,3                                                    
  23  AM(IND1)=AM(IND1)+M(IND2)*BR2(IND1,IND2)                          
      ABSM=ABSM+AM(IND1)**2                                             
      AHELP=AM(IND1)/PIA(IND1)                                          
      IM(IND1)=AHELP+SIGN(0.01d0,AHELP)                                 
      IF(ORTHO.or.lattic(1:3).eq.'CXZ') GOTO 22              
      IM(1)=M(1)                                                        
      IM(2)=M(2)                                                        
      IM(3)=M(3)
   22 CONTINUE                                                          
!  transformation into primitiv monoclinic basis
      IF(.not.ORTHO.and.lattic(1:3).eq.'CXZ') then
        im(1)=m(1)+m(3)
        im(2)=m(2)
        im(3)=m(1)-m(3)
      endif 
!                                                                       
       ABSM=SQRT(ABSM)                                                  
      if(absm.gt.gmax) goto 1
!
      IF(NK.EQ.0) GOTO 7                                                
      DO 2 J=1,NK                                                       
      ABST=ABSM-ABSK(J)                                                 
      IF(ABS(ABST).LT.TEST)GOTO 5                                       
      IF(ABST.LT.0.d0) GOTO 4                                            
   2  CONTINUE                                                          
! .... NEW VECTOR GT. OLD MAX                                           
      J=NK                                                              
      IF(J.GE.NK1) then
          nkd=nkd*3
          nk1=(nk1-1)*3+1
          call doreallocate(absk,nkd+1)
          call doreallocate(inst,nkd)
          call doreallocate(kzz,3,nkd+1)
          endif
      J=J+1                                                             
      GOTO 7                                                            
  5   CALL STERN                                                        
! .... NEW VECTOR EQU. OLD ONE                                          
      DO 3 I=J,NK                                                       
      ABST=ABSM-ABSK(I)                                                 
      IF(ABS(ABST).GT.2.d0*TEST) GOTO 44                                       
      IF(KDELTA(KZZ(1,I))) GOTO 1                                       
  3   CONTINUE                                                          
      i=nk                          
 44   j=i
  4   NJ=NK-J+1                                                         
! .... NEW VECTOR .NE. OLD ONE                                          
      DO 6 JJ=1,NJ                                                      
      NJJ=NK-JJ+1                                                       
      ABSK(NJJ+1)=ABSK(NJJ)                                             
      DO 6 N=1,3                                                        
  6   KZZ(N,NJJ+1)=KZZ(N,NJJ)                                           
   7  ABSK(J)=ABSM                                                      
! .... PUT NEW VECTOR IN LIST                                           
      DO 8 N=1,3                                                        
   8  KZZ(N,J)=IM(N)                                                    
!      KZZ(3,J)=IABS(KZZ(3,J))                                          
      NK=NK+1                                                           
!      IF(NK.GE.NK1) NK=NK1-1                                            
      IF(NK.GE.NK1) then                                                
          nkd=nkd*3
          nk1=(nk1-1)*3+1
          call doreallocate(absk,nkd+1)
          call doreallocate(inst,nkd)
          call doreallocate(kzz,3,nkd+1)
          endif
   1  CONTINUE                                                          
!
      IND=0                                                             
      DO 30 I=1,NKK
      DO 31 J=1,3                                                       
   31 IM(J)=KZZ(J,I)                                                    
      CALL STERN                                                        
      INST(I)=NST                                                       
      IF(IND+NST.GT.NFK) THEN
        nfk=3*nfk
        call doreallocate(tauk,nfk)
        call doreallocate(krec,3,nfk) 
!      WRITE(66,102) IND,I                                               
! 102  FORMAT('0 MORE THAN NFK K REQUESTED, TRUNCATED TO',I5,            &
!      '   NKK+1=',I5,//)                                                
!      GOTO 33                                                           
      END IF                                                            
      DO 32 JJ=1,NST                                                    
      IND=IND+1                                                         
      TAUK(IND)=TAUP(JJ)
      DO  JJJ=1,3                                                     
        KREC(JJJ,IND)=ISTM(JJJ,JJ)  
      enddo
 32   continue
   30 CONTINUE                                                          
!
 33   NKK=I-1                                                           
      INDMAX=IND                                                        
      RETURN                                                            
  100 FORMAT(I4)                                                        
  101 FORMAT(3(3I2,F10.5/))                                              
      END                                                               
