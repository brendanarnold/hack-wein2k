      SUBROUTINE WRTSCF(nat,iter,jspin,minmod,temp2,q)
!
!     read scf file and sums up total energy and forces
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
!
      dimension      temp2(maxit)
      dimension      q(nat,2,maxit)
      CHARACTER*79   MARGN                                             
      CHARACTER*4    minmod
!
!.....GOTO BOTTOM OF TAPE23=SCFMINI                                 
!
      ISCF=0                                                            
 110  READ(23,700,END=120) MARGN
      IF(MARGN(2:4).EQ.'ITE') then
	 READ(MARGN,715) ISCF
	 index=1
      endif
      IF(MARGN(2:6).EQ.'CTO  ') goto 110
      IF((MARGN(2:4).EQ.'CTO').and.(jspin.eq.1)) then
	 READ(MARGN,760) q(index,1,iscf)
	 index=index+1
      endif
      IF(MARGN(2:6).EQ.'CUP  ') goto 110
      IF(MARGN(2:4).EQ.'CUP') then
	 jspin=2
	 READ(MARGN,760) q(index,1,iscf)
	 if (index.eq.nat) index=0
	 index=index+1
      endif
      IF(MARGN(2:6).EQ.'CDN  ') goto 110
      IF(MARGN(2:4).EQ.'CDN') then
	 READ(MARGN,760) q(index,2,iscf)
	 index=index+1
      endif
      GO TO 110                                                          
 120  CONTINUE  
      Backspace(23)                                                        
!
!.....READ TAPE22=SCFDATA UNTIL BOTTOM                                  
!                                                                       
      ISCF=0                                                            
  10  READ(22,700,END=20) MARGN
      IF(MARGN(2:4).EQ.'ITE')GO TO 30                                          
      GO TO 10                                                          
  30  READ(MARGN,710) ISCF 
      GOTO 10                                                           
  20  CONTINUE                                                          
!
!..... READ RIGHT SCF SECTION
!
      REWIND 22                                                         
 60   READ(22,700,END=40) MARGN
      IF (MARGN(2:4).NE.'ITE') GO TO 60                                        
      READ(MARGN,710) ISCFT 
!      write (*,*)'iter,iscf,iscft',iter,iscf,iscft
      IF (ISCFT.NE.ISCF) GO TO 60                                        
      WRITE(24,705) iter,iter
      READ(22,700,END=40) MARGN
      index=1
 50   READ(22,700,END=40) MARGN
      write(24,700) MARGN
      IF(MARGN(2:6).EQ.'CTO  ') goto 50
      IF((MARGN(2:4).EQ.'CTO').and.(jspin.eq.1)) then
	 READ(MARGN,760) q(index,1,iter)
	 index=index+1
      endif
      IF(MARGN(2:6).EQ.'CUP  ') goto 50
      IF(MARGN(2:4).EQ.'CUP') then
	 jspin=2
	 READ(MARGN,760) q(index,1,iter)
!	 write(6,*) 'cup',index,iter,q(index,1,iter)
	 if (index.eq.nat) index=0
	 index=index+1
      endif
      IF(MARGN(2:6).EQ.'CDN  ') goto 50
      IF(MARGN(2:4).EQ.'CDN') then
	 READ(MARGN,760) q(index,2,iter)
	 index=index+1
      endif
      GO TO 50                                                          
 40   CONTINUE                                                          
      if (minmod.eq.'NOSE'.or.minmod.eq.'MOLD') &
             write (24,750) temp2(iter)
      RETURN
!
!
!
 700  FORMAT(A79)                                                    
 705  FORMAT(//,12X,9(1H-),/,':ITE',i3.3,':',I3,'. ITERATION',/,12X, &
             9(1H-),/)
 710  FORMAT(4X,i3)                                              
 715  FORMAT(4X,i3)                                              
 720  FORMAT(30X,f20.6)                                              
 730  FORMAT(20X,f20.6)                                             
 740  FORMAT(15X,4f15.3)                                             
 750  FORMAT(':TEM',2x,':',1x,'TEMPERATURE (K): ',f10.5)
 760  FORMAT(38X,f10.7)                                             
!
      END 
