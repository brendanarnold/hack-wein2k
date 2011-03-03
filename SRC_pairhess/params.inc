!      Maximum number of atoms
       parameter (natmax=300)
!      Maximum number of neighbours
       parameter (neigmax=200)
!      Exponential decay part
       logical useother,Simple,isfree(3,natmax),depends(3,natmax)
!      Common for symmetry operations
       integer icode(5,natmax*48)
       dimension iz(3,3,48),tau(3,48)
       common /symms/icode,iz,tau,iord,nat,index
       common /cell/alat(3),scl,cutoff
       common /plogicals/Simple,isfree,depends
       dimension wmult(natmax),mult(natmax)
       common /which/wmult,mult
       integer indeq(5,natmax,48),nineq(natmax)
       common /aequiv/atrans(4,natmax,neigmax),indeq,nineq
