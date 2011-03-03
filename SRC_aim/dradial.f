      subroutine dradial(clm,drho,r,ir,jatom)
      use rad
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      dimension clm(nrad)
!      common /rad/ rm(nrad,nato),jri(nato)
      r1=rm(ir,jatom)                                                   
      r2=rm(ir+1,jatom)                                                 
      dr=r2-r1                                                          
      ddr=r-r1                                                          
      c1=clm(ir)                                                        
      c2=clm(ir+1)                                                      
      dc=c2-c1                                                          
      rho=c1+(ddr/dr)*dc                                                
      drho=-2.D0*rho/(r*r*r)+dc/dr/(r*r)      
      return                                                            
      end                                                               
