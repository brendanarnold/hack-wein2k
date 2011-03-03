!BOP
! !ROUTINE: Somm
! !INTERFACE:
    subroutine somm(dr,dp,dq,dpas,da,m,np)                           
! !INPUT/OUTPUT PARAmETERS:
!   dr     :   radial mesh used for integration
!   dpas   :   exponential step of the radial mesh
!   dp,dq  :   large and small components of the radial wave function
!              of the core state, expressed on mesh dr
!   np     :   integrate up to r=dr(np)
!   m      :   specify metric as r^m
!   da     :  I don't know
! !DESCRIPTION:
!   INTEGRATION PAR LA mETHODE DE SImPSON DE (dp+dq)*dr**m DE 0 A R=dr(NP)
!   POUR R VOISIN DE ZERO (dp+dq)=CTE*R**da        
! !REVISION HISTORY:
!   Originally taken from SRC_lcore/
!   Updated November 2004 (Kevin Jorissen)
!   Cleaned June 2005 (Kevin Jorissen)
!EOP

      implicit none
!  IN/OUTPUT
      real*8 dr,dp,dq,dpas,da
      dimension dr(*),dp(*),dq(*)                                       
	  integer m,np
!  LOCAL VARIABLES
      integer mm,i
	  real*8 db,dl,dc,d1



      mm=m+1                                                            
      d1=da+mm                                                          
      da=dble(0)                                                             
      db=dble(0)                                                             
      do i=1,np                                                      
        dl=dr(i)**mm                                                      
        if (.not.(i.eq.1.or.i.eq.np)) then
          dl=dl+dl                                                          
          if ((i-2*(i/2)).eq.0) dl=dl+dl
	    endif                                    
        dc=dp(i)*dl
		if(dc.lt.0) then
		  db=db+dc
		elseif(dc.gt.0) then
		  da=da+dc
		endif                                                       
        dc=dq(i)*dl
		if (dc.lt.0) then
		  db=db+dc
		elseif (dc.gt.0) then
		  da=da+dc
		endif                                                       
      enddo                                                          
      da=dpas*(da+db)/dble(3)                                                
      dc=dexp(dpas)-dble(1)                                                   
      db=d1*(d1+dble(1))*dc*dexp((d1-dble(1))*dpas)                                
      db=dr(1)*(dr(2)**m)/db                                            
      dc=(dr(1)**mm)*(dble(1)+dble(1)/(dc*(d1+dble(1))))/d1                            
      da=da+dc*(dp(1)+dq(1))-db*(dp(2)+dq(2))                           

      return                                                            
      end                                                               
