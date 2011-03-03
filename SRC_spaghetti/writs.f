      subroutine writs(string,size,strnam)
      use irr_param
      IMPLICIT REAL*8 (A-H,O-Z)
!      INCLUDE 'param.inc'
      character*12     string
      character*6      strnam
      character*40     sfont(7), font
!      common /prtbc/   ibnd(NEVL),iprto(3),iprtf,ifont 
      common /prtft/   sfont
!
!-------------------------------------------------------------------
!
!.....no text
      if(ifont.eq.0) then
        return
      endif
!
!.....select fonts
      if(strnam(1:4).eq.'text') then
         font=sfont(1)
      else if(strnam(1:4).eq.'symb') then
         font=sfont(2)
      else if(strnam(1:4).eq.'ferm') then
         font=sfont(3)
      else if(strnam(1:4).eq.'labe') then
         font=sfont(4)
      else if(strnam(1:4).eq.'head') then
         font=sfont(5)
      else
         stop 'incorrect strnam' 
      endif

!.....(in order to not writing out k-point number)
      if (string(3:5).eq.'   ') then
        write(11,100) font,size*20.d0, string(1:2)
      return
      endif
!
      if (string(1:5).eq.'GAMMA') then  
         font=sfont(2)
         write(11,110) font,size*20.d0,'G'
      return
      endif
      if (string(1:5).eq.'DELTA') then 
         font=sfont(2)
         write(11,110) font,size*20.d0,'D'
      return
      endif
      if (string(1:6).eq.'LAMBDA') then 
         font=sfont(2)
         write(11,110) font,size*20.d0,'L'
      return
      endif
      if (string(1:5).eq.'SIGMA') then 
         font=sfont(2)
         write(11,110) font,size*20.d0,'S'
      return
      endif
!
      do i=12,1,-1
      if (string(i:i).ne.' ') then
         write(11,120) font,size*20.d0,string(1:i)
      return
      endif
      enddo  
      write(11,130) font,size*20.d0,string
!
      return
 100  format(A40,' findfont',f10.5,' scalefont setfont',/, &
             '(',a2,') show')
 110  format(A40,' findfont',f10.5,' scalefont setfont',/, &
             '(',a1,') show')
 120  format(A40,' findfont',f10.5,' scalefont setfont',/, &
             '(',a,') show')
 130  format(A40,' findfont',f10.5,' scalefont setfont',/, &
             '(',a12,') show')
!
      end

