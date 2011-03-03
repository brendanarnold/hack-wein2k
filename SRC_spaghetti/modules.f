module irr_param
      dimension      skir(3),iprto(3)
      real*8,pointer ::        eneik(:,:)                          !NEVL,NKP)
      integer,pointer ::       ngrp(:,:,:), ngde(:,:,:)            !(4,NEVL,NKP)
      integer,pointer ::       ndeg(:,:)                           !NEVL,NKP)
      integer,pointer ::       neik(:),ikg(:),ntab4(:)             !NKP)
      integer,pointer ::       iscol(:),iscol2(:),iste(:),ibnd(:)  !NEVL)
      integer,pointer ::       lkg(:,:,:)                          !48,48,NKP)
      integer,pointer ::       nlkg(:,:)                           !48,NKP)
      character*12,pointer ::  grpnam(:)                           !NKP) 
      character*12,pointer ::  cnam(:,:)                           !48,NKP)
      character*11    col
      character*3    tjunk
      integer iprtf,ifont 
      real*8 dlwid,dawid,dfwid 
end module irr_param
