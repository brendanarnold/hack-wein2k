function watsonw(r,rw,ew)

  IMPLICIT REAL*8 (A-H,O-Z)
  
  if (r.lt.rw) then
     wa=-ew/rw
  else
     wa=-ew/r
  endif
  
  watsonw=wa 
  
end function watsonw
