!BOP
! !ROUTINE: DumpWorkspace
! !INTERFACE:
      subroutine DumpWorkspace(n)
! !USES:
      use constants
	  use dimension_constants
	  use units
	  use work_atpar
	  use energygrid
      use radial_functions
	  use potential
	  use input
	  use program_control
	  use qvectors
	  use initialstate1
	  use initialstate2
      use leftovers
	  use densityofstates
      use struct
	  use cross_dos
	  use rotation_matrices
! !INPUT/OUTPUT PARAMETERS:
!     n   :   integer specifying a particular module
! !DESCRIPTION:
!    This routine is not really a part of the ELNES package; it is only a debugging tool.
!    Since modules are not visible in my Visual Fortran Debugger, I have made a routine
!    that dumps modules to file 97.
!    It can be called simply as Call DumpWorkspace(n), where n specifies a particular module.
!    See the debug.f file for the list of modules.
!    Warning : may not have been fully updated.
! !REVISION HISTORY:
!   Created November 2004 (Kevin Jorissen)
!EOP
	  implicit none
	  integer n
	  integer j,k,l

    open(97,file='debug.txt',form='formatted',position='append')
    if (n.eq.1) then
	 write(97,*) 'module constants:'
	 write(97,*) 'pi= ',pi
	elseif(n.eq.2) then
	 write(97,*) 'module dimension_constants:'
	 write(97,*) 'nrad= ',nrad
	 write(97,*) 'lmax= ',lmax
	 write(97,*) 'nofcross= ',nofcross
	elseif(n.eq.3) then
       write(97,*) 'module units:'
	 write(97,*) 'a0= ',a0
	elseif(n.eq.4) then
       write(97,*) 'module work_atpar:'
	 write(97,*) 'a= ',(a(j),j=1,nrad)
	 write(97,*) 'a1= ',((a1(j,k),j=1,nrad),k=1,iemax)
	 write(97,*) 'ae= ',(ae(j),j=1,nrad)
	 write(97,*) 'b= ',(a(j),j=1,nrad)
	 write(97,*) 'b1= ',((a1(j,k),j=1,nrad),k=1,iemax)
	 write(97,*) 'be= ',(ae(j),j=1,nrad)
      elseif(n.eq.5) then
	 write(97,*) 'module energygrid:'
	 write(97,*) 'ef= ',ef
	 write(97,*) 'jemin= ',jemin
	 write(97,*) 'jemax= ',jemax
	 write(97,*) 'iemax= ',iemax
	 write(97,*) 'ene= ',(ene(j),j=0,iemax)
	elseif(n.eq.6) then
     write(97,*) 'module radial_functions'
	 write(97,*) 'uell= ',(((uell(j,k,l),j=1,nrad),k=1,iemax),l=0,lmax)
	 write(97,*) 'upell= ',(((upell(j,k,l),j=1,nrad),k=1,iemax),l=0,lmax)
	 write(97,*) 'uelluell= ',(((uelluell(j,k,l),j=1,iemax),k=0,lmax),l=0,lmax)
	elseif(n.eq.7) then
	 write(97,*) 'module potential'
	 write(97,*) 'vr= ',((vr(j,k),k=1,nat),j=1,nrad)
	elseif(n.eq.8) then
	 write(97,*) 'module input:' 
	 write(97,*) 'natom= ',natom
	 write(97,*) 'neqat= ',neqat
	 write(97,*) 'neqatdeb= ',neqatdeb
	 write(97,*) 'nc= ',nc
	 write(97,*) 'lc= ',lc
	 write(97,*) 'emin= ',emin
	 write(97,*) 'de= ',de
	 write(97,*) 'emax= ',emax
	 write(97,*) 'split= ',split
	 write(97,*) 'deltae= ',deltae
	 write(97,*) 'energy= ',energy
	 write(97,*) 'k0len= ',k0len
     write(97,*) 'thetax= ',thetax
	 write(97,*) 'thetay= ',thetay
	 write(97,*) 'aalpha= ',aalpha
	 write(97,*) 'abeta= ',abeta
	 write(97,*) 'agamma= ',agamma
	 write(97,*) 'specaperture= ',specaperture
	 write(97,*) 'convergangle= ',convergangle
	 write(97,*) 'thpart= ',thpart
	 write(97,*) 'th0= ',th0
	 write(97,*) 'selrule= ',selrule
	 write(97,*) 'title= ',title
	 write(97,*) 'corewffile= ',corewffile
	 write(97,*) 'finalwffile= ',finalwffile
	 write(97,*) 'xqtlfile= ',xqtlfile
	 write(97,*) 'nr= ',nr
	 write(97,*) 'nt= ',nt
	 write(97,*) 'npos= ',npos
	 write(97,*) 'branchingratio= ',branchingratio
	 write(97,*) 'w1= ',w1
	 write(97,*) 'w2= ',w2
	elseif(n.eq.9) then
	 write(97,*) 'module program_control'
	 write(97,*) 'verbosity= ',verbosity
	 write(97,*) 'modus= ',modus
	 write(97,*) 'qmodus= ',qmodus
	 write(97,*) 'mrot= ',mrot
	 write(97,*) 'mdos= ',mdos
	 write(97,*) 'calcsplit= ',calcsplit
	 write(97,*) 'file2core= ',file2core
	 write(97,*) 'file2final= ',file2final
	 write(97,*) 'xqtlalaNovak= ',xqtlalaNovak
	elseif(n.eq.10) then
	 write(97,*) 'module qvectors:'
	 write(97,*) 'qv= ',(((qv(k,j,l),k=1,3),j=1,mult(natom)),l=1,npos)
	 write(97,*) 'qlenv= ',((qlenv(k,j),k=1,mult(natom)),j=1,npos)
	 write(97,*) 'weightv= ',weightv
	 write(97,*) 'thxv= ',thxv
	 write(97,*) 'thyv= ',thyv
	elseif(n.eq.11) then
	 write(97,*) 'module initialstate1:'
	 write(97,*) 'dpas= ',dpas
	 write(97,*) 'z= ',z
	 write(97,*) 'tets= ',tets
	 write(97,*) 'nstop= ',nstop
	 write(97,*) 'nes= ',nes
	 write(97,*) 'np= ',np
	 write(97,*) 'nuc= ',nuc
	 write(97,*) 'dv= ',(dv(j),j=1,nrad+41)
	 write(97,*) 'dr= ',(dr(j),j=1,nrad+41)
	 write(97,*) 'dp= ',(dp(j),j=1,nrad+41)
	 write(97,*) 'dq= ',(dq(j),j=1,nrad+41)
	elseif(n.eq.12) then
	 write(97,*) 'module initialstate2:'
	 write(97,*) 'dpdp= ',(dpdp(j),j=1,nrad+41)
	 write(97,*) 'dqdq= ',(dqdq(j),j=1,nrad+41)
	elseif(n.eq.13) then
	 write(97,*) 'module leftovers:'
	 write(97,*) 'gammaa0= ',gammaa0
	 write(97,*) 'cfein1= ',cfein1
	 write(97,*) 'cfein2= ',cfein2
	elseif(n.eq.14) then
	 write(97,*) 'module densityofstates:'
	 write(97,*) 'indexlm= ',((indexlm(j,k),k=-lmax,lmax),j=0,lmax)
	 write(97,*) 'dos= ',((dos(j,k),k=0,lmax),j=1,iemax)
	elseif(n.eq.15) then
     write(97,*) 'module spectra: removed from list'
	elseif(n.eq.16) then
	 write(97,*) 'module data_inilpw'
	 write(97,*) 'nat= ',nat
	 write(97,*) 'nats= ',nats
	 write(97,*) 'rel= ',rel
	 write(97,*) 'alat= ',alat(1:3)
	 write(97,*) 'alpha= ',alpha(1:3)
	 write(97,*) 'zz= ',zz(1:nat)
	 write(97,*) 'iatnr= ',iatnr(1:nat)
	 write(97,*) 'mult= ',mult(1:nat)
	 write(97,*) 'isplit= ',isplit(1:nat)
	 write(97,*) 'rmt= ',rmt(1:nat)
	 write(97,*) 'pos= ',((pos(j,k),j=1,3),k=1,nats)
	 write(97,*) 'jrj= ',(jrj(j),j=1,nat)
	 write(97,*) 'dh= ',(dh(j),j=1,nat)
	 write(97,*) 'rO= ',(rO(j),j=1,nat)
	elseif(n.eq.17) then
	 write(97,*) 'module cross_dos:'
	 write(97,*) 'to be done'
	elseif(n.eq.18) then
	 write(97,*) 'rotation_matrices:'
       write(97,*) 'rotij= ',(((rotij(l,k,j),k=1,3),l=1,3),j=1,nats)
	 write(97,*) 'rotloc= ',(((rotloc(l,k,j),k=1,3),l=1,3),j=1,nats)
	 write(97,*) 'generalm= ',(((generalm(l,k,j),k=1,3),l=1,3),j=1,nats)
	else
	 write(97,*) 'module ',n,' not described yet.'
	endif


      write(97,*);write(97,*)
      close(97)
	return
	end
