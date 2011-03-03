!BOP
! !ROUTINE: Insld
! !INTERFACE:
      SUBROUTINE INSLD(RELA,npxx,jatom,nprim,lprim,kappa)                                            
! !USES:
      use program_control,only : verbosity
	  use dimension_constants,only : nrad
      use initialstate1
      use struct, izi => iz
! !INPUT/OUTPUT PARAMETERS:
!   rela     :   switch for relativistic or non-relativistic calculation
!   npxx     :   radial mesh has jrj(jatom)+npxx points
!   jatom    :   atom for which core state is calculated
!   nprim    :   main quantum number of core state 
!   lprim    :   orbital quantum number of core state
!   kappa    :   relativistic quantum number of core state
! !DESCRIPTION:
!   Read the spherical potential from file 18 and initializes a few parameters
!   for the calling routine hfsd.
! !REVISION HISTORY:
!     Taken from SRC_lcore.
!     Heavily modified November 2004 (Kevin Jorissen)
!EOP
	  implicit none
!  INPUT/OUTPUT
      integer jatom,nprim,lprim,kappa,npxx
	  logical rela
!  LOCALS
      CHARACTER *10 BAR
	  integer nqn,nk,nql,nmax
	  real*8 nel
	  real*8 den,dexv,dq1,dfl,dexe,dcop,test,teste,testy,testv,ion,ion1
	  real*8 dvc,dsal,dr1,dval,val
	  integer iz1,nuc1,j,l,niter,i,k
      COMMON DEN,DQ1,DFL,NQN,NQL,NK,NMAX,NEL                                                         
      COMMON/PS2/DEXV,DEXE,DCOP,TEST,TESTE,TESTY,TESTV,ION,NITER,TITRE
      COMMON /CHAR/ BAR
	  character ttire(9),titre
      DATA TTIRE/4HS   ,4HP*  ,4HP   ,4HD*  ,4HD   ,4HF*  ,4HF   ,4HG*  ,4HG   /                                                          

1     IF (NSTOP.ne.0) then               
        DVC=dble(137.0373)                                                      
        IF(.NOT.RELA) DVC=1.E30                                           
        DSAL=DVC+DVC                                                      
        IZ1=0                                                             
        ION1=0                                                            
        NUC1=-1                                                           
        IF (NSTOP.ne.1) goto 598
	  endif
      DPAS=dble(0.05)                                                         
      DR1=dble(0.01)
      TESTE=5.E-06
      TESTY=1.E-05
      TESTV=1.E-05
      TEST=dble(1.E-07)
      NSTOP=30
      DEXV=dble(2)/dble(3)
      DEXE=dble(1)
      DCOP=dble(0.3)

! **********************************************************************

      np=jrj(jatom)
	  dpas=dh(jatom)
	  dr1=rO(jatom)
      z=zz(jatom)
	  bar=name(jatom)

      I=NP
      J=30
      L=0
      NUC=0
      NES=J                                          
      NITER=1                                               

! Z NUMERO ATOMIQUE    ION=Z-NOMBRE D ELECTRONS                       
! I NCMBRE DE POINTS POUR L INTEGRATION  NRAD PAR DEFAUT                 
! J NCMBRE D ESSAIS POUR AJUSTER L ENERGIE  15 PAR DEFAUT               
! NITER NOMBRE D ITERATIONS  50 PAR DEFAUT                                  
! L=0 OPTION STANDARD POUR LE BLOC DES POINTS ET LES PRECISIONS         
! NOYAU DE DIMENSIONS FINIES SI NUC POSITIF                             
! **********************************************************************
                                                                       
! DPAS PAS EXPONENTIEL   DR1 DEFINIT LE PREMIER POINT =DR1/Z           
! TEST PRECISION SUR L ENERGIE DANS RESLS                               
! TESTE CRITERE DE SELF CONSISTENCE POUR LES ENERGIES MONDELECTRONIQUES 
! TESTY CRITERE DE SELF CONSISTENCE POUR LES FONCTIONS D ONDE           
! TESTV CRITERE DE SELF CONSISTENCE POUR LE POTENTIEL                   
! **********************************************************************
      if(verbosity.ge.2) then
        WRITE (6,1990) BAR,NITER,TESTE,TESTY,TESTV                        
        WRITE (6,2008) NP,DR1,Z,DPAS                                     
        WRITE (6,2009) TEST,NES
      endif
                                                            
      DVAL=Z*Z/(DVC*DVC)                                                
                                                                       
! DONNEES POUR LES ORBITALES                                            
! **********************************************************************

	  nqn=nprim
	  nk=kappa
	  nel=dble(1)

! DEN ENERGIE DE L ORBITALE EN UNITE ATOMIQUE ET NEGATIVE               
! NQN NOMBRE QUANTIQUE PRINCIPAL
! NK NOMBRE QUANTIQUE KAPPA            
! NEL OCCUPATION DE L ORBITALE                                          
! **********************************************************************
                                                    
      DEN=-Z*Z/(4.*NQN*NQN)                                    
      NQL=IABS(NK)                                                
      IF (NK.LT.0) NQL=NQL-1 
      if(nql.ne.lprim) stop ':WARNING lprim != nql in insld'
      DFL=DSQRT(NK*NK-DVAL)                                          
      L=2*IABS(NK)                                                   
      J=NQL+IABS(NK)
      TITRE=TTIRE(J)
      WRITE (6,2004) NQN,TITRE,NEL,DEN

      ION=Z-NEL
      IF (0.NE.NUC1) GO TO 29                                         
      IF (Z.EQ.IZ1.AND.ION.EQ.ION1) GO TO 101                          
      IF(Z.EQ.IZ1) GO TO 35                                            
 29   DR(1)=DR1
      do  I=2,NRAD+41
        DR(I)=DR(1)*DEXP((I-1)*DPAS)
	  enddo
                                                                       
! POTENTIEL DE DEPART                                                   
! **********************************************************************
 35   VAL=-ION-1

                                                    
      READ(18,'(x,//)')
      do k=1,jatom
        READ(18,'(x,///)')
        READ(18,'(/)')
        READ(18,2021) (DV(I),I=1,NP)
        READ(18,'(/////)')
      enddo
2021  FORMAT(3X,4E19.12)                                                 

      NPXX=40                                                           
      I=2*(np/2)+1  
      if(i.gt.np) npxx=41                                                     
      DO I=1,NP+NPXX                                                 
        IF(I.GT.NP)DV(I)=DV(NP)                                           
        DV(I)= DV(I)/2./DR(I)
      enddo                                             
      NP=NP+NPXX                                                        
! **********************************************************************
                                                             
 101  IZ1=Z                                                            
      ION1=ION                                                          
      NUC1=NUC                                                          
      NMAX=NP                                                        
      L=1                                                               
      J=NQN-NQL          
      IF((J-2*(J/2)).EQ.0) L=-L                                         
      DQ1=L*NK/IABS(NK)                                        

      RETURN                                                            



! KJ I think you only come here in case of trouble :
 598  if(verbosity.ge.2) WRITE (6,2018)                                                    
      NSTOP=1                                                           
      GO TO 1


 1990 FORMAT(A10,/,5X,'NOMBRE D ITERATIONS',I4,/,5X,'PRECISION SUR LES ENERGIES',&
      1PE9.2,/,23X,'FONCTIONS D ONDE',1PE9.2,/,23X &
      ,'POTENTIEL',1PE9.2)                                            
 2004 FORMAT (2X,'ORBITALE',3x,i1,a2,5X,'OCCUPATION',3x,f8.4,5X,'ENERGIE DE DEPART',3x,1PE23.7)
 2008 FORMAT (' L INTEGRATION EST FAITE EN ',I3,' POINTS LE PREMIER EST EGAL A ',e10.4,'/',f7.2,' ET LE PAS A ',F7.4)                      
 2009 FORMAT (' DANS LE SOUS PROGRAMME RESLD LA PRECISION RELATIVE A OBTENIR SUR L ENERGIE EST ', &
      1PE9.2,' ET LE NOMBRE D ESSAIS ',I3,/)                         
 2018 FORMAT ('  CAS SUIVANT')                                          


      END                                                               
