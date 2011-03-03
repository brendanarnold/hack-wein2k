#!/bin/csh -f
unalias rm

set running = ".running.$$.`hostname`.`date +%d%m%H%M%S`"
echo $$ >$running
onintr clear
alias error	'echo ">>> ($name) \!* -> exit"; goto error'

set name	= $0		#full name of script-file
set bin		= $name:h	#default directory for WIEN-executables
if !(-d $bin) set bin = .
set name	= $name:t

set def 	= def		#def-file
set tmp		= ( )		#tmp-files, to be deleted before exit
set log		= :log		#log-file
set t		= time		#
set updn			# spinpolarization switch
set dnup	= 'dn'		# spinpolarization switch
set settol=0.00001              # tol parameter for sgroup
set sc				# semicore-switch
set so				# spinorbit-switch
set sodum			# spinorbit-switch
set cmplx
set para                        # parallel execution 
set dopara                        # parallel execution 
set sel                         # standard vector file in lapw7
set band1                       # regular run in lapw1
set eece                       
unset orb                      # LDA+U in lapw1
unset nmat_only                 # determines matrixsize in lapw1
unset iter                      # iterative diagonalization
unset nohns                     # without HNS in lapw1
unset qtl                       # regular run in lapw2
unset sigma                      # regular run in lstart
unset band                       # regular run in lapw2
unset fermi                     # regular run in lapw2
unset efg                       # regular run in lapw2
unset command			#command
unset file			#file-head
unset deffile_only		#create just def file
unset fftstop			#for dstart	
set scratch =                   #eg: /scr1/pblaha/   for vectors, help-files,

if ( $?SCRATCH ) then
  set scratch=`echo $SCRATCH  | sed -e 's/\/$//'`/ # we are afraid
				# different settings in different 
				# computing centers
			        #use global variable for scratch if set
endif

set argv1=($argv)

while ($#argv)
  switch ($1)
  case -f:
    shift; set file = $1
    shift; breaksw
  case -h:
  case -H: 
    set help
    shift; breaksw
  case -[t|T]:
    set t
    shift; breaksw
  case -up:
    set updn = 'up'
    set dnup = 'dn'
    shift; breaksw
  case -dn:
    set updn = 'dn'
    set dnup = 'up'
    shift; breaksw
  case -du:
    set du = 'du'
    set updn = 'up'
    shift; breaksw
  case -sc:
    set sc = 's'
    shift; breaksw
  case -d:
    set deffile_only
    shift; breaksw
  case -fft:
    set fftstop
    shift; breaksw
  case -c:
    set cmplx = c
    shift; breaksw
  case -[p|P]
    set para = para
    shift; breaksw
  case -it:
    set iter
    shift; breaksw	
  case -eece:
    set eece=eece
    shift; breaksw	
  case -orb:
    set orb
    shift; breaksw	
  case -nmat_only:
    set nmat_only
    shift; breaksw	
  case -nohns:
    set nohns
    shift; breaksw	
  case -so:
    set so = 'so'
    set sodum = 'dum'
    set cmplx = c
    shift; breaksw
  case -band:
    set band
    set band1='_band'
    shift; breaksw
  case -qtl:
    set qtl
    shift; breaksw
  case -sigma:
    set sigma
    shift; breaksw
  case -efg:
    set efg
    shift; breaksw
  case -sel:
    set sel = f
    shift; breaksw
  case -fermi:
    set fermi
    shift; breaksw  
  case -settol:
    shift;set settol = $1
    shift; breaksw  
  default:
    set optiontest = `echo $1|cut -c1`
    if ( "$optiontest" == '-' ) then
        echo "error in your arguments:  $1   is not a valid option"
        exit(3)
    endif
    set command = $1
    shift; breaksw
  endsw
end

if !($?file) then 
set file    = `pwd`
set file    = $file:t		#tail of file-names
endif


if ( $?help ) goto help

if !($?command) then
  error no command
endif

if ($command == lapw1c) then
  set cmplx = c
  set command = lapw1
else if ($command == lapw1it) then
  set iter
  set command = lapw1
else if ($command == lapw1itc) then
  set cmplx = c
  set iter
  set command = lapw1
else if ($command == lapw2c) then
  set cmplx = c
  set command = lapw2
else if ($command == lapwdmc) then
  set cmplx = c
  set command = lapwdm
else if ($command == dmatc) then
  set cmplx = c
  set command = dmat
else if ($command == filtvecc) then
  set cmplx = c
  set command = filtvec
else if ($command == lapw7c) then
  set cmplx = c
  set command = lapw7
else if ($command == lapw5c) then
  set cmplx = c
  set command = lapw5
else if ($command == spinorbitc) then
  set cmplx = c
  set command = spinorbit
else if ($command == opticc) then
  set cmplx = c
  set command = optic
endif

echo "`date`> ($name) $command $argv1[2-]" >> $log

set def	= $updn$command$sc.def
#touch $def

switch ($command)


case lcore:
set exe = lcore
cat << theend > $def
 5,'$file.inc',         'old',    'formatted',0
 6,'$file.outputc$updn','unknown','formatted',0
 8,'$file.vsp$updn',    'old',    'formatted',0
 9,'$file.clmcor$updn', 'unknown','formatted',0
19,'$file.vns$updn',    'unknown','formatted',0
20,'$file.struct',      'old',    'formatted',0         
21,'$file.scfc$updn',   'unknown','formatted',0
28,'$file.vrespcor$updn',   'unknown','formatted',0
29,'$file.corewf$updn',   'unknown','formatted',0
theend
breaksw

case dstart:
set exe = dstart
   if ($updn == '') then
cat << theend > $def
 6,'$file.outputd',  'unknown','formatted',0
13,'$file.in0_std','unknown',    'formatted',0
14,'$file.in0','old',    'formatted',0
15,'$file.in2$cmplx','old',    'formatted',0
17,'$file.in1$cmplx','old',    'formatted',0
20,'$file.struct',   'old',    'formatted',0
51,'$file.clmsum',   'unknown','formatted',0
81,'$file.rsp',	     'old',    'formatted',0
theend
   else
cat << theend > $def
 6,'$file.outputd$updn','unknown','formatted',0
13,'$file.in0_std','unknown',    'formatted',0
14,'$file.in0','old',    'formatted',0
15,'$file.in2$cmplx',   'old',    'formatted',0
17,'$file.in1$cmplx',   'old',    'formatted',0
16,'$file.test',        'unknown','formatted',0
20,'$file.struct',      'old',    'formatted',0
51,'$file.clm$updn',    'unknown','formatted',0
81,'$file.rsp$updn',    'old',    'formatted',0
theend
   endif
if ($?fftstop) then
echo " 4,'$file.fftstop',  'unknown','formatted',0">>$def
endif
breaksw

case kgen:
set exe = kgen
cat << theend > $def
 8,'$file.klist',     'unknown','formatted',0    
15,'$file.kgen',      'unknown','formatted',0      
66,'$file.outputkgen','unknown','formatted',0
theend
if ($so == so ) then
echo "20,'$file.ksym',      'old','formatted',0" >>$def
else
echo "20,'$file.struct',    'old','formatted',0" >>$def   
endif
breaksw

case lapw0: 
set exe = lapw0$para
set dopara = 1
set tmp = ($file.poissn $file.r2v $file.vcoul)
if (-e $file.clmup &&  ! -z $file.clmup) then
    set updn = 'up'
endif
set endup = 'up'
if ($eece == 'eece') set endup = 'valupeece'
set enddn = 'dn'
if ($eece == 'eece') set enddn = 'valdneece'

cat << theend > $def
 3,'$file.rhopw',   'unknown','formatted',0
 4,'$file.inm',     'unknown','formatted',0
 5,'$file.in0$eece',     'old',    'formatted',0
 6,'$file.output0', 'unknown','formatted',0
 7,'$file.vorbup',  'unknown','formatted',0
 8,'$file.clmsum',  'old',    'formatted',0
 9,'$file.vtotal',  'unknown','formatted',0
10,'$file.vcoul',   'unknown','formatted',0
11,'$file.r2v',     'unknown','formatted',0
12,'$file.clm$endup',   'unknown','formatted',0
13,'$file.clm$enddn',   'unknown','formatted',0
16,'$file.vsp$updn','unknown','formatted',0
17,'$file.vspdn',   'unknown','formatted',0
18,'$file.vns$updn','unknown','formatted',0
19,'$file.vnsdn',   'unknown','formatted',0
20,'$file.struct',  'old',    'formatted',0         
21,'$file.scf0',    'unknown','formatted',0
28,'$file.vrespsum','unknown','formatted',0
29,'$file.vrespup', 'unknown','formatted',0
30,'$file.vrespdn', 'unknown','formatted',0
50,'$file.eeceup', 'unknown','formatted',0
51,'$file.eecedn', 'unknown','formatted',0
theend
breaksw


case lapw1:
#if ($?band) then
#  set old="`head -1 $file.in1$cmplx |cut -c1-5`"
#  cp $file.in1$cmplx .oldin1
#  sed "s/UNIT:./UNIT:5/" .oldin1 >$file.in1$cmplx
#endif
set exe = $command$cmplx$para
set dopara = 1
set nmat_string
if ($?nmat_only) set nmat_string='_nmat_only'
cat << theend > $def
 4,'$file.klist$sc${band1}',          'unknown','formatted',0
 5,'$file.in1$cmplx$sc',   'old',    'formatted',0
 6,'$file.output1$sc$updn$nmat_string','unknown','formatted',0
10,'${scratch}$file.vector$sc$updn$nmat_string', 'unknown','unformatted',9000
11,'$file.energy$sc$updn$nmat_string', 'unknown','formatted',0
18,'$file.vsp$updn',       'old',    'formatted',0
19,'$file.vns$updn',       'unknown','formatted',0
20,'$file.struct',         'old',    'formatted',0
21,'$file.scf1$sc$updn$nmat_string',   'unknown','formatted',0
55,'$file.vec',            'unknown','formatted',0
71,'$file.nsh$sc$updn$nmat_string',    'unknown','formatted',0
theend
if ($?orb) echo " 7,'$file.vorb$updn'     ,'unknown','formatted',0" >> $def
if ($?nmat_only) echo "72,'$file.nmat_only'     ,'unknown','formatted',0" >> $def
if ($?iter) echo "98,'${scratch}$file.vector$sc$updn.old', 'unknown','unformatted',9000" >> $def
if ($?nohns) echo "97,'${scratch}$file.nohns', 'unknown','formatted',9000" >> $def
breaksw


case lapw2:
if ($?band) then
  set old="`head -1 $file.in2$cmplx |cut -c1-5`"
  cp $file.in2$cmplx .oldin2
  sed "3s/^...../ROOT /" .oldin2 >.oldin2t
  sed "1s/^...../QTL  /" .oldin2t >$file.in2$cmplx
  rm .oldin2t
  unset qtl
endif
if ($?qtl) then
  set old="`head -1 $file.in2$cmplx |cut -c1-5`"
  cp $file.in2$cmplx .oldin2
  sed "1s/^...../QTL  /" .oldin2 >$file.in2$cmplx
endif
if ($?efg) then
  set old="`head -1 $file.in2$cmplx |cut -c1-5`"
  cp $file.in2$cmplx .oldin2
  sed "1s/^...../EFG  /" .oldin2 >$file.in2$cmplx
endif
if ($?fermi) then
  set old="`head -1 $file.in2$cmplx |cut -c1-5`"
  cp $file.in2$cmplx .oldin2
  sed "1s/^...../FERMI/" .oldin2 >$file.in2$cmplx
endif
set exe = $command$cmplx$para
set dopara = 1
set ending = 'val'
if ($sc == 's') set ending = 'sc'
if ($sodum != 'dum') then
  set sodum=$dnup
endif
cat << theend > $def
 2,'$file.nsh$sc$updn',    'unknown','formatted',0
 3,'$file.in1$cmplx$sc',   'unknown','formatted',0
 4,'$file.inso',           'unknown','formatted',0
 5,'$file.in2$cmplx$sc$eece',   'old',    'formatted',0
 6,'$file.output2$sc$updn$eece','unknown','formatted',0
 8,'$file.clm$ending$updn$eece','unknown','formatted',0
10,'${scratch}$file.vector$sc$so$updn', 'unknown','unformatted',9000
11,'$file.weight$sc$updn',    'unknown','formatted',0
13,'$file.recprlist',      'unknown','unformatted',9000
14,'$file.kgen$sc',        'unknown','formatted',0
15,'$file.tmp$updn',       'unknown','formatted',0
16,'$file.qtl$updn',       'unknown','formatted',0
17,'$file.weightaver$so$updn','unknown','formatted',0
18,'$file.vsp$updn',       'old',    'formatted',0
19,'$file.vns$updn',       'unknown','formatted',0
20,'$file.struct',         'old',    'formatted',0
21,'$file.scf2$sc$updn',   'unknown','formatted',0
22,'$file.rotlm',   'unknown',    'formatted',0
23,'$file.radwf',   'unknown',    'formatted',0
24,'$file.almblm',   'unknown',    'formatted',0
26,'$file.weigh$updn',   'unknown','unformatted',0
27,'$file.weigh$dnup',   'unknown','unformatted',0
28,'$file.vresp$ending$updn',   'unknown','formatted',0
29,'$file.energy$sc$sodum','unknown','formatted',0
30,'$file.energy$sc$so$updn', 'unknown','formatted',0
31,'${scratch}$file.help$updn', 'unknown','formatted',0
theend
if ($so == so ) echo "12,'$file.norm$sc$so$updn',    'unknown','formatted',0" >>$def
breaksw

case lapwso:
set dnup1
if ( $updn == up ) set dnup1=dn
set exe = lapwso$para
set dopara = 1
set def = lapwso.def
cat << theend > $def
4 ,'$file.in1$cmplx$sc',   'old',    'formatted',0   
5 ,'$file.inso', 'old',    'formatted',0
6 ,'$file.outputso',   'unknown','formatted',0
8 ,'$file.scfso',       'unknown','formatted',0
9 ,'${scratch}$file.vector$dnup1',    'old',    'unformatted',9000
10 ,'${scratch}$file.vectorup',    'unknown',    'unformatted',9000
18,'$file.vsp$dnup1',  'old','formatted',0
19,'$file.vspup',  'unknown','formatted',0
20 ,'$file.struct',    'old',    'formatted',0
22,'$file.vns$dnup1',  'old','formatted',0
23,'$file.vnsup',  'unknown','formatted',0
41,'${scratch}$file.vectorsodn',  'unknown','unformatted',9000
42,'${scratch}$file.vectorso$updn',  'unknown','unformatted',9000
44,'$file.vect1',  'unknown','unformatted',9000
45,'$file.normsodn',  'unknown','formatted',0
46,'$file.normsoup',  'unknown','formatted',0
51,'$file.energysodn',  'unknown','formatted',9000
52,'$file.energyso$updn',  'unknown','formatted',9000
53,'$file.energydum',  'unknown','formatted',9000
54,'$file.energy$dnup1',  'old','formatted',9000
55 ,'$file.energyup',    'unknown',    'formatted',9000
theend
if($?orb) echo "11,'$file.vorbdn',  'unknown','formatted',0" >> $def
if($?orb) echo "12,'$file.vorbup',  'unknown','formatted',0" >> $def
if($?orb) echo "13,'$file.vorbdu',  'unknown','formatted',0" >> $def
breaksw

case lapwdm:
set exe = $command$cmplx$para
set dopara = 1
cat << theend > $def
 4,'$file.inso',   'unknown',    'formatted',0
 5,'$file.indm$cmplx',   'old',    'formatted',0
 6,'$file.outputdm$updn','unknown','formatted',0
 7,'$file.dmat$updn','unknown','formatted',0
 8,'$file.dmat$dnup','unknown','formatted',0
 9,'${scratch}$file.vector$sc$so$updn', 'unknown','unformatted',9000
10,'${scratch}$file.vector$sc$so$dnup', 'unknown','unformatted',9000
14,'$file.kgen$sc',        'unknown','formatted',0
18,'$file.vsp$updn',       'old',    'formatted',0
19,'$file.vsp$dnup',       'unknown','formatted',0
20,'$file.struct',         'old',    'formatted',0
21,'$file.scfdm$updn',   'unknown','formatted',0
26,'$file.weigh$updn',   'unknown','unformatted',0
50,'$file.energy$sc$so$updn', 'unknown','formatted',9000
51,'$file.energy$sc$so$dnup', 'unknown','formatted',9000
theend
if ($?so) echo "11,'$file.dmatud'     ,'unknown','formatted',0" >> $def
breaksw

case orb:
set pip
set exe = orb
set updn_old=$updn{_old}
if ( $para == para ) then
    set pip = '_1'
    set para
    endif  
cat << theend > $def
 5,'$file.inorb',         'old',    'formatted',0
 6,'$file.outputorb$updn',   'unknown','formatted',0
 9,'$file.dmat$dnup',      'unknown','formatted',0
10,'$file.dmat$updn',      'unknown','formatted',0
20,'$file.struct',        'old',    'formatted',0
31,'$file.br1orb$updn',    'unknown','unformatted',0
32,'$file.br2orb$updn',    'unknown','unformatted',0
theend
if ($?du) then
cat << theend >> $def
11,'$file.dmatud'     ,'unknown','formatted',0
12,'$file.vorbdu',     'unknown','formatted',0
13,'$file.vorbdu_old',     'unknown','formatted',0
14,'$file.energyup$pip',      'unknown','formatted',0
18,'$file.vspup',       'unknown','formatted',0
21,'$file.scforbdu',      'unknown','formatted',0
theend
else
cat << theend >> $def
12,'$file.vorb$updn',     'unknown','formatted',0
13,'$file.vorb$updn_old',     'unknown','formatted',0
14,'$file.energy$updn$pip',      'unknown','formatted',0
18,'$file.vsp$updn',       'unknown','formatted',0
21,'$file.scforb$updn',      'unknown','formatted',0
50,'$file.eece$updn',    'unknown','formatted',0
theend
endif

breaksw

case averx:
set exe = averx
set def = averx.def
cat << theend > $def
5 ,'$file.inaverx', 'old',    'formatted',0
6 ,'$file.outputaverx',   'unknown','formatted',0
8 ,'$file.scfaverx',       'unknown','formatted',0
9 ,'${scratch}$file.vectordn',    'old',    'unformatted',9000
10 ,'${scratch}$file.vectorup',    'old',    'unformatted',9000
16,'$file.weightaversoup',  'old','formatted',0
18,'$file.vspdn',  'old','formatted',0
19,'$file.vspup',  'old','formatted',0
20 ,'$file.struct',    'old',    'formatted',0
26,'$file.weightaverdn',  'unknown','formatted',0
27,'$file.weightaverup',  'unknown','formatted',0
41,'${scratch}$file.vectorsodn',  'unknown','unformatted',9000
42,'${scratch}$file.vectorsoup',  'unknown','unformatted',9000
43,'${scratch}$file.vectordum',  'unknown','unformatted',9000
44,'$file.vect1',  'unknown','unformatted',9000
theend
breaksw

case optimize:
set exe = optimize
  touch optimize.job
  chmod +x optimize.job
cat << theend > $def
20,'${file}_initial.struct',      'unknown',    'formatted',0         
16,'optimize.job','unknown','formatted',0
17,'$file.struct',      'old',    'formatted',0         
theend
breaksw

case eosfit:
set exe = eosfit
cat << theend > $def
55,'$file.vol',        'old',    'formatted',0         
66,'$file.outputeos',     'unknown','formatted',0         
9,'$file.eosfit',     'unknown','formatted',0         
11,'$file.eosfitb',     'unknown','formatted',0         
theend
breaksw

case eosfit6:
set exe = eosfit6
cat << theend > $def
10,'$file.ene',        'old',    'formatted',0         
11,'$file.latparam',        'old',    'formatted',0         
12,'$file.enefit',     'unknown','formatted',0         
66,'$file.outputeos6',     'unknown','formatted',0         
theend
breaksw

case optic:
set exe = $command$cmplx$para
set dopara = 1
cat << theend > $def
4, '${scratch}$file.mommat$updn' ,  'UNKNOWN',    'FORMATTED',  0
5, '$file.inop'     ,      'OLD',    'FORMATTED',  0
6, '$file.outputop$updn' ,  'UNKNOWN',    'FORMATTED',  0
3, '${scratch}$file.symmat$updn'   ,  'UNKNOWN',    'FORMATTED',  0
9, '${scratch}$file.mat_diag$updn' ,  'UNKNOWN',    'FORMATTED',  0
10,'${scratch}$file.vector$sc$so$updn'   ,      'OLD',  'UNFORMATTED',  0
11,'${scratch}$file.vector$sc$so$dnup'  ,   'UNKNOWN'  , 'UNFORMATTED' ,0
18,'$file.vsp$updn'      ,      'OLD',    'FORMATTED',  0
19,'$file.vsp$dnup'     ,    'UNKNOWN' , 'FORMATTED',  0
20,'$file.struct'   ,      'OLD',    'FORMATTED',  0
28,'$file.inso'    ,     'UNKNOWN',    'FORMATTED' , 0
24,'${scratch}$file.mme$updn'       ,  'UNKNOWN',   'FORMATTED', 0
25,'$file.symop'     ,  'UNKNOWN',   'FORMATTED', 0
theend

breaksw

case nlo_core:
set exe = $command
cat << theend > $def
3,'$file.mme$updn',  'OLD',   'FORMATTED', 0
4,'$file.klist'   ,  'OLD','FORMATTED',  0
10,'$file.kgen'   ,  'OLD','FORMATTED',  0
17,'$file.weight$updn','OLD','FORMATTED',  0
20,'$file.meta$updn' ,'UNKNOWN',    'FORMATTED',  0
21,'$file.innlo_core'     ,'OLD',    'FORMATTED',  0
25,'$file.symop'     ,'OLD',   'FORMATTED', 0
theend

breaksw

case nlo_tet:
set exe = $command
cat << theend > $def
4,'$file.klist'   ,  'OLD','FORMATTED',  0
10,'$file.kgen'   ,  'OLD','FORMATTED',  0
17,'$file.innlo_tet',    'OLD','FORMATTED',  0
18,'$file.nlo$updn' ,'UNKNOWN',    'FORMATTED',  0
20,'$file.meta$updn','OLD',    'FORMATTED',  0
theend

breaksw


case joint:
set exe = $command
cat << theend > $def
 3,'${scratch}$file.symmat$updn' ,  'OLD','FORMATTED',  0
 4,'$file.weight$updn' ,  'OLD','FORMATTED',  0
 5,'$file.injoint'   ,  'OLD','FORMATTED',  0
 6,'$file.outputjoint$updn',  'UNKNOWN','FORMATTED',  0
 7,'$file.joint$updn'  ,  'UNKNOWN','FORMATTED',  0
 8,'$file.sigma_intra$updn'  ,  'UNKNOWN','FORMATTED',  0
9, '${scratch}$file.mat_diag$updn' ,  'UNKNOWN',    'FORMATTED',  0
11,'$file.intra$updn'  ,  'UNKNOWN','FORMATTED',  0
14,'$file.kgen'   ,  'OLD','FORMATTED',  0
20,'$file.struct' ,  'OLD',    'FORMATTED',  0
theend

breaksw

case kram:
set exe = $command
cat << theend > $def
 5,'$file.inkram'   ,  'OLD','FORMATTED',  0
10,'$file.joint$updn'  ,  'OLD','FORMATTED',  0
11,'$file.intra$updn'  ,  'UNKNOWN','FORMATTED',  0
12,'$file.epsilon$updn'  ,  'UNKNOWN','FORMATTED',  0
13,'$file.sigmak$updn'  ,  'UNKNOWN','FORMATTED',  0
14,'$file.absorp$updn'  ,  'UNKNOWN','FORMATTED',  0
15,'$file.eloss$updn'  ,  'UNKNOWN','FORMATTED',  0
16,'$file.refraction$updn',  'UNKNOWN','FORMATTED',0
17,'$file.reflectivity$updn','UNKNOWN','FORMATTED',0
77,'$file.sumrules$updn'  ,  'UNKNOWN','FORMATTED',  0
theend

breaksw

case sumpara:
set exe = sumpara
set deffile_only
cat << theend > $def
 6,'$file.outputsum',   'unknown','formatted',0
 8,'$file.scfdm$updn',    'unknown','formatted',0
17,'$file.clmval$updn$eece',    'unknown','formatted',0
20,'$file.struct',    'old',    'formatted',0
21,'$file.scf2$updn',      'unknown','formatted',0
22,'$file.scf2p',       'unknown','formatted',0
theend
if ($?du) then
cat << theend >> $def
10,'$file.dmatup'     ,'unknown','formatted',0
11,'$file.dmatdn'     ,'unknown','formatted',0
12,'$file.dmatud'     ,'unknown','formatted',0
theend
else
cat << theend >> $def
10,'$file.dmat$updn',    'unknown','formatted',0
theend
endif

breaksw


case sumpara_vresp:
set exe = sumpara
set deffile_only
cat << theend > $def
 6,'$file.outputsum',   'unknown','formatted',0
17,'$file.vrespval$updn',    'unknown','formatted',0
20,'$file.struct',    'old',    'formatted',0
21,'$file.scf2$updn',      'unknown','formatted',0
22,'$file.scf2p',       'unknown','formatted',0
theend

breaksw

case lstart:
set exe = lstart  
#touch $file.inst_sigma
if ($?sigma) then
#  if(  -z $file.inst_sigma ) then
    sed "s/ N/ P/" $file.inst >$file.inst_sigma
#  endif
endif
cat << theend > $def
 6,'$file.outputst',      'unknown','formatted'
 7,'$file.tmpden',        'unknown','formatted'
10,'$file.tmp',	          'unknown','unformatted'
11,'$file.vsp_st',           'unknown','formatted'
12,'$file.vspdn_st',         'unknown','formatted'
13,'$file.sigma',         'unknown','formatted'
14,'$file.in0_st',        'unknown','formatted'
15,'$file.in1_st',        'unknown','formatted'
16,'$file.inc_st',        'unknown','formatted'
17,'$file.in2_ls',        'unknown','formatted'
18,'$file.inm_st',        'unknown','formatted'
19,'$file.inm_restart_st','unknown','formatted'
20,'$file.struct',        'old',    'formatted'
81,'$file.rspup',         'unknown','formatted'
82,'$file.rspdn',         'unknown','formatted'
83,'$file.rsp',           'unknown','formatted'
94,'$file.rsigma',         'unknown','formatted'
theend
if ($?sigma) then
  echo " 5,'$file.inst_sigma',          'old',    'formatted'" >>$def
else
  echo " 5,'$file.inst',          'old',    'formatted'" >>$def
endif
breaksw
 
case mini:
set exe = mini
if !(-e $file.clmvalup || -z $file.clmvalup) then
cat << theend > $def
5 ,'$file.inM',       'old','formatted',0
6 ,'$file.outputM',   'unknown','formatted',0
20,'$file.struct',    'old',    'formatted',0
10,'$file.clmsum',    'unknown',    'formatted',0
13,'$file.clmhist',   'unknown',    'formatted',0
15,'$file.finM',      'old',    'formatted',0
16,'$file.tmpM',      'unknown','formatted',0
17,'$file.tmpM1',     'unknown','formatted',0
21,'$file.struct1',   'unknown','formatted',0
22,'$file.scf','unknown','formatted',0
23,'$file.scf_mini','unknown','formatted',0
24,'$file.scf_mini1','unknown','formatted',0
25,'$file.constraint','unknown','formatted',0
51,'$file.clmsum_inter',    'unknown',    'formatted',0
theend
        else
cat << theend > $def
5 ,'$file.inM',       'old','formatted',0
6 ,'$file.outputM',   'unknown','formatted',0
20,'$file.struct',    'old',    'formatted',0
10,'$file.clmup',    'old',    'formatted',0
60,'$file.clmdn',    'old',    'formatted',0
13,'$file.clmhist',      'unknown',    'formatted',0
15,'$file.finM',      'old',    'formatted',0
16,'$file.tmpM',      'unknown','formatted',0
17,'$file.tmpM1',      'unknown','formatted',0
21,'$file.struct1',   'unknown','formatted',0
22,'$file.scf','unknown','formatted',0
23,'$file.scf_mini','unknown','formatted',0
24,'$file.scf_mini1','unknown','formatted',0
11,'$file.clmsum_inter',    'unknown',    'formatted',0
51,'$file.clmup_inter',    'unknown',    'formatted',0
52,'$file.clmdn_inter',    'unknown',    'formatted',0
theend
endif
breaksw


case mixer:
set exe = mixer
	if !(-e $file.clmvalup || -z $file.clmvalup) then
cat << theend > $def
 5,'$file.inm',       'old',    'formatted',0
 6,'$file.outputm',   'unknown','formatted',0
 7,'$file.inc',       'old',    'formatted',0
10,'$file.clmsum_old','unknown','formatted',0
17,'$file.clmval',    'unknown','formatted',0
18,'$file.clmsc',     'unknown','formatted',0
19,'$file.clmcor',    'unknown','formatted',0
20,'$file.struct',    'old',    'formatted',0
21,'$file.scfm',      'unknown','formatted',0
22,'$file.scf',       'unknown','formatted',0
31,'$file.broyd1',    'unknown','unformatted',164896
32,'$file.broyd2',    'unknown','unformatted',164896
51,'$file.clmsum',    'unknown','formatted',0
71,'$file.dmatup','unknown','formatted',0
72,'$file.dmatdn','unknown','formatted',0
73,'$file.dmatud','unknown','formatted',0
74,'$file.dmatup_old','unknown','formatted',0
75,'$file.dmatdn_old','unknown','formatted',0
76,'$file.dmatud_old','unknown','formatted',0
theend
	else
cat << theend > $def
 5,'$file.inm',      'old',    'formatted',0
 6,'$file.outputm',  'unknown','formatted',0
 7,'$file.inc',       'old',    'formatted',0
10,'$file.clmup_old','unknown','formatted',0
60,'$file.clmdn_old','unknown','formatted',0
17,'$file.clmvalup', 'unknown','formatted',0
18,'$file.clmscup',  'unknown','formatted',0
19,'$file.clmcorup', 'unknown','formatted',0
20,'$file.struct',   'old',    'formatted',0
21,'$file.scfm',     'unknown','formatted',0
22,'$file.scf',      'unknown','formatted',0
31,'$file.broyd1',   'unknown','unformatted',164896
32,'$file.broyd2',   'unknown','unformatted',164896
47,'$file.clmvaldn', 'unknown','formatted',0
48,'$file.clmscdn',  'unknown','formatted',0
49,'$file.clmcordn', 'unknown','formatted',0
51,'$file.clmup',    'unknown','formatted',0
52,'$file.clmdn',    'unknown','formatted',0
11,'$file.clmsum',   'unknown','formatted',0
71,'$file.dmatup','unknown','formatted',0
72,'$file.dmatdn','unknown','formatted',0
73,'$file.dmatud','unknown','formatted',0
74,'$file.dmatup_old','unknown','formatted',0
75,'$file.dmatdn_old','unknown','formatted',0
76,'$file.dmatud_old','unknown','formatted',0
theend
	endif
   breaksw
#77,'$file.dmatup','unknown','formatted',0
#78,'$file.dmatdn','unknown','formatted',0
#79,'$file.dmatdu','unknown','formatted',0


case mixer_vresp:
set exe = mixer
	if !(-e $file.clmvalup || -z $file.clmvalup) then
cat << theend > $def
 5,'$file.inm_vresp',       'old',    'formatted',0
 6,'$file.outputm',   'unknown','formatted',0
 7,'$file.inc',       'old',    'formatted',0
10,'$file.vrespsum_old','unknown','formatted',0
17,'$file.vrespval',    'unknown','formatted',0
18,'$file.vrespsc',     'unknown','formatted',0
19,'$file.vrespcor',    'unknown','formatted',0
20,'$file.struct',    'old',    'formatted',0
31,'$file.vrespbroyd1',    'unknown','unformatted',164896
32,'$file.vrespbroyd2',    'unknown','unformatted',164896
51,'$file.vrespsum',    'unknown','formatted',0
theend
	else
cat << theend > $def
 5,'$file.inm_vresp',      'old',    'formatted',0
 6,'$file.outputm',  'unknown','formatted',0
 7,'$file.inc',       'old',    'formatted',0
10,'$file.vrespup_old','unknown','formatted',0
60,'$file.vrespdn_old','unknown','formatted',0
17,'$file.vrespvalup', 'unknown','formatted',0
18,'$file.vrespscup',  'unknown','formatted',0
19,'$file.vrespcorup', 'unknown','formatted',0
20,'$file.struct',   'old',    'formatted',0
31,'$file.vrespbroyd1',   'unknown','unformatted',164896
32,'$file.vrespbroyd2',   'unknown','unformatted',164896
47,'$file.vrespvaldn', 'unknown','formatted',0
48,'$file.vrespscdn',  'unknown','formatted',0
49,'$file.vrespcordn', 'unknown','formatted',0
51,'$file.vrespup',    'unknown','formatted',0
52,'$file.vrespdn',    'unknown','formatted',0
11,'$file.vrespsum',   'unknown','formatted',0
theend
	endif
   breaksw


case nn:
set exe = nn
cat << theend > $def
66,'$file.outputnn','unknown','formatted',0
67,'$file.nnshells','unknown','formatted',0
20,'$file.struct',  'old',    'formatted',0
21,'$file.struct_nn','unknown','formatted',0
theend
breaksw

case dipan:
set exe = dipan
cat << theend > $def
5 ,'$file.indipan',    'old','formatted',0
6 ,'$file.outputdipan','unknown','formatted',0
20,'$file.struct',     'old',    'formatted',0
theend
breaksw

case tetra:
set exe = tetra  
cat << theend > $def
4 ,'$file.qtl$updn',    'old','formatted',0
5 ,'$file.int',         'old','formatted',0
6 ,'$file.outputt$updn','unknown','formatted',0
3 ,'$file.kgen',        'old',    'formatted',0
7 ,'$file.dos1$updn',   'unknown','formatted',0
57,'$file.dos1ev$updn', 'unknown','formatted',0
theend
breaksw

case lapw5:
set exe = lapw5$cmplx  
cat << theend > $def
5 ,'$file.in5$cmplx', 'old',    'formatted',0
6 ,'$file.output5',   'unknown','formatted',0
8 ,'$file.struct',    'old',    'formatted',0
9 ,'$file.clmval$updn',    'old',    'formatted',0
10,'$file.tmp',       'unknown','unformatted',0
11,'$file.clmval$dnup',  'unknown','formatted',0
12,'$file.sigma',     'unknown','formatted',0
20,'$file.rho_onedim','unknown','formatted',0
21,'$file.rho',       'unknown','formatted',0
theend
breaksw

case filtvec:
set exe = filtvec$cmplx
cat << theend > $def
5 ,'$file.inf$cmplx',   'old',    'formatted',0
6 ,'$file.outputf$updn','unknown','formatted',0
10,'${scratch}$file.vector$updn', 'old',    'unformatted',0
12,'$file.vectorf$updn','unknown','unformatted',0
20,'$file.struct',      'old',    'formatted',0
theend
breaksw

case lapw7:
set exe = lapw7$cmplx
set vecloc = $scratch; if ($sel == f) set vecloc =
cat << theend > $def
2 ,'$file.tmp'             'unknown','formatted',0
5 ,'$file.in7$cmplx',      'old',    'formatted',0
6 ,'$file.output7$updn',   'unknown','formatted',0
7 ,'$file.grid',           'unknown','formatted',0
8 ,'$file.struct',         'old',    'formatted',0
10,'${vecloc}$file.vector$sel$updn','old',    'unformatted',0
12,'$file.abc$updn',       'unknown','unformatted',0
18,'$file.vsp$updn',       'old',    'formatted',0
21,'$file.psink',          'unknown','formatted',0
theend
breaksw

case lapw3:
set exe = $command$cmplx  
cat << theend > $def
66,'$file.output3',  'unknown','formatted',0
4 ,'$file.in2$cmplx','old',    'formatted',0
8 ,'$file.clmsum',   'old',    'formatted',0
15,'$file.tmp3',      'unknown','unformatted',0
20,'$file.struct',   'old',    'formatted',0
theend
breaksw

case momentum:
set exe = $command$cmplx  
cat << theend > $def
66,'$file.outputmom',  'unknown','formatted',0
4 ,'$file.in2$cmplx','old',    'formatted',0
8 ,'$file.clmsum',   'old',    'formatted',0
10,'${scratch}$file.vector$sc$so$updn', 'unknown','unformatted',9000
15,'$file.tmp',      'unknown','unformatted',0
20,'$file.struct',   'old',    'formatted',0
23,'$file.radwf',   'old',    'formatted',0
24,'$file.almblm',   'old',    'formatted',0
theend
breaksw

case symmetry:
set exe = symmetry  
cat << theend > $def
6, '$file.outputs','unknown','formatted',0
17,'$file.in2_sy', 'unknown','formatted',0
20,'$file.struct', 'old',    'formatted',0
21,'$file.struct_st', 'unknown',    'formatted',0
theend
breaksw

case symmetso:
set exe = symmetso  
cat << theend > $def
 5,'$file.inso',   'old',    'formatted',0
 6,'$file.outsymso','unknown','formatted',0
25,'$file.vspdn',         'unknown',    'formatted',0
45,'$file.vspdn_so',         'unknown',    'formatted',0
26,'$file.vspup',         'unknown',    'formatted',0
46,'$file.vspup_so',         'unknown',    'formatted',0
27,'$file.vnsdn',         'unknown',    'formatted',0
47,'$file.vnsdn_so',         'unknown',    'formatted',0
28,'$file.vnsup',         'unknown',    'formatted',0
48,'$file.vnsup_so',         'unknown',    'formatted',0
20,'$file.struct_interm',         'unknown',    'formatted',0
21,'$file.struct_so',         'unknown',    'formatted',0
22,'$file.struct'           'old',    'formatted',0
23,'$file.ksym',         'unknown',    'formatted',0
24,'$file.temp',         'unknown',    'formatted',0
29,'$file.in1$cmplx',         'unknown',    'formatted',0
49,'$file.in1${cmplx}_so',         'unknown',    'formatted',0
30,'$file.inc',         'unknown',    'formatted',0
50,'$file.inc_so',         'unknown',    'formatted',0
31,'$file.inorb',         'unknown',    'formatted',0
51,'$file.inorb_so',         'unknown',    'formatted',0
32,'$file.vorbdn',         'unknown',    'formatted',0
52,'$file.vorbdn_so',         'unknown',    'formatted',0
33,'$file.vorbup',         'unknown',    'formatted',0
53,'$file.vorbup_so',         'unknown',    'formatted',0
34,'$file.in2$cmplx',         'unknown',    'formatted',0
54,'$file.in2${cmplx}_so',         'unknown',    'formatted',0
35,'$file.clmsum',         'unknown',    'formatted',0
55,'$file.clmsum_so',         'unknown',    'formatted',0
36,'$file.clmup',         'unknown',    'formatted',0
56,'$file.clmup_so',         'unknown',    'formatted',0
37,'$file.clmdn',         'unknown',    'formatted',0
57,'$file.clmdn_so',         'unknown',    'formatted',0
theend
breaksw


case spaghetti:
set exe = spaghetti  
set typ=1
if ($so == so ) then
    set typ=so
endif
if($para == para) then
    set cnum=0
    if(-f .processes) then
	@ cnum=`grep -ce "^[0-9][0-9]* :" .processes`
#	@ cnum+=1
    endif
    #the original file will only be removed if there are sources to concat
    if(-f "$file.output$typ$updn" && $cnum>0) then
	rm "$file.output$typ$updn"
    endif
    unset irr
    if(-f "$file.irrep$so${updn}_$cnum" && $cnum>0) then
	if(-f "$file.irrep$so${updn}" ) rm "$file.irrep$so${updn}"
        head -2 "$file.irrep$so${updn}_$cnum" >"$file.irrep$so${updn}"
        set irr
    endif
    set count=1
    while($count <= $cnum)
        if($?irr) then
           tail +3 $file.irrep$so${updn}_$count >>$file.irrep$so${updn}
	endif
        cat "$file.output$typ${updn}_$count" >>"$file.output$typ$updn"
	@ count++
    end
endif
cat << theend > $def
5, '$file.insp',         'old',    'formatted',0
6, '$file.outputsp$updn',     'unknown','formatted',0
9, '$file.qtl$updn',          'unknown','formatted',0
10,'$file.spaghetti${updn}_ene','unknown','formatted',0
11,'$file.spaghetti${updn}_ps', 'unknown','formatted',0
20,'$file.struct',       'old',    'formatted',0
30,'$file.irrep$so$updn',          'unknown','formatted',0
40,'$file.bands${updn}.agr',       'unknown','formatted',0
theend
#if ($so == so ) then
#echo "7, '$file.outputso',      'old','formatted',0" >>$def
#else
echo "7, '$file.output$typ$updn',  'old','formatted',0" >>$def
#endif
breaksw


case xspec:
set exe = $command
cat << theend > $def
 5,'$file.inxs',   'old',    'formatted',0
 6,'$file.outputx','unknown','formatted',0
 7,'$file.inc',   'old',    'formatted',0
 8,'$file.int','unknown','formatted',0
 9,'$file.corewfx$updn', 'unknown','formatted',0
18,'$file.vsp$updn',       'old',    'formatted',0
20,'$file.struct',         'old',    'formatted',0
30,'$file.qtl$updn', 'old',    'formatted',0
32,'$file.dos1ev$updn',        'unknown','formatted',0
46,'$file.xspec$updn',            'unknown','formatted',0
47,'$file.txspec$updn',            'unknown','formatted',0
53,'$file.m1$updn',            'unknown','formatted',0
54,'$file.m2$updn',            'unknown','formatted',0
theend
breaksw

case init_xspec:
set exe = $command
cat << theend > $def
 5,'$file.inxs',   'old',    'formatted',0
 8,'$file.int','unknown','formatted',0
30,'$file.qtl$updn', 'old',    'formatted',0
theend
breaksw

case txspec:
set exe = $command
cat << theend > $def
 5,'$file.inxs',   'old',    'formatted',0
 6,'$file.outputx','unknown','formatted',0
 7,'$file.inc',   'old',    'formatted',0
 9,'$file.corewfx$updn',   'unknown','formatted',0
18,'$file.vsp$updn',        'old',    'formatted',0
20,'$file.struct',          'old',    'formatted',0
21,'$file.scfc$updn',       'unknown','formatted',0
30,'$file.qtl$updn',        'old',    'formatted',0
32,'$file.dos1ev$updn',     'old','formatted',0
47,'$file.txspec$updn',     'unknown','formatted',0
53,'$file.m1$updn',         'unknown','formatted',0
54,'$file.m2$updn',         'unknown','formatted',0
theend
breaksw

case lorentz:
set exe = $command
cat << theend > $def
 5,'$file.inxs',             'old',    'formatted',0
46,'$file.xspec$updn',       'unknown','formatted',0
47,'$file.txspec$updn',      'old',    'formatted',0
theend
breaksw

case telnes2:
set exe = $command
cut -f1 -d' ' $file.innes > .temp_making_elnes.def
if (`grep -i xqtl .temp_making_elnes.def` != ) then
  set number = `head -n 2 $file.innes | tail -n 1`
  if ($#number > 1) set number = $number[1]
else
  set number =
endif
cat << theend > $def
3, '$file.kgen','unknown','formatted',0
5, '$file.innes',           'old',    'formatted',0
6, '$file.outputelnes$updn',       'unknown','formatted',0
7, '$file.inc',             'unknown',    'formatted',0
8, '$file.inb','unknown','formatted',0
9, '$file.corewavef$updn',   'unknown','formatted',0
10,'$file.final$updn',              'unknown','formatted',0
18,'$file.vsp$updn',        'unknown',    'formatted',0
20,'$file.struct',          'old',    'formatted',0
30,'$file.qtl$number$updn',        'unknown',    'formatted',0
31,'$file.rotij',           'unknown',    'formatted',0
46,'$file.ortho$updn',         'unknown','formatted',0
47,'$file.elnes$updn',        'unknown','formatted',0
48,'$file.ctr$updn',        'unknown','formatted',0
49,'$file.sdlm$updn',        'unknown','formatted',0
50,'$file.matrix$updn','unknown','formatted',0
56,'$file.cdos$updn','unknown','formatted',0
57,'$file.dos$updn',     'unknown',    'formatted',0
58,'$file.xdos$updn','unknown','formatted',0
59,'$file.sp2$updn','unknown','formatted',0
60,'$file.angular$updn','unknown','formatted',0
theend
breaksw

case broadening:
set exe = broadening
cat << theend > $def
 5,'$file.inb','old','formatted',0
 6,'$file.outputbroadening','unknown','formatted',0
46,'$file.broadspec','unknown','formatted',0
47,'$file.elnes$updn','old','formatted',0
theend
breaksw

case hex2rhomb:
set exe = $command
breaksw

case rhomb_in5:
set exe = $command
breaksw

case veccopy:
if($sc == 's') set ending = 'sc'
set exe = veccopy$cmplx
cat << theend > $def
5, '$file.inveccopy',    'old',    'formatted', 0
6, '$file.outputveccopy','unknown','formatted', 0
20,'$file.struct',       'old',    'formatted', 0
10,'${scratch}$file.vectorup', 'old','unformatted',9000
30,'${scratch}$file.vectordn', 'unknown','unformatted',9000
71,'$file.nshup',    'unknown','formatted',0
81,'$file.nshdn',    'unknown','formatted',0
theend
breaksw

case clmcopy:
set ending = 'val'
set updn = 'up'
if($sc == 's') set ending = 'sc'
set exe = clmcopy
cat << theend > $def
5, '$file.inclmcopy',    'old',    'formatted', 0
6, '$file.outputclmcopy','unknown','formatted', 0
7, '$file.dmatup',       'unknown','formatted', 0
8, '$file.dmatdn',       'unknown','formatted', 0
17,'$file.clm$ending$updn',     'old',    'formatted',  0
18,'$file.clm$ending$dnup',     'unknown','formatted',  0
20,'$file.struct',       'old',    'formatted', 0
21,'$file.scf2$sc$updn',       'unknown',    'formatted',       0
22,'$file.scf2$sc$dnup',       'unknown',    'formatted',       0
theend
breaksw


case RMTCheck:
set exe = RMTCheck$cmplx
   if ($updn == '') then
cat << theend > $def
 6,'$file.outputRMT',  'unknown','formatted',0
 8,'$file.struct',   'old',    'formatted',0
 9,'$file.clmval',   'unknown','formatted',0
theend
   else
cat << theend > $def
 6,'$file.outputRMT$updn','unknown','formatted',0
 8,'$file.struct',      'old',    'formatted',0
 9,'$file.clmval$updn',    'unknown','formatted',0
theend
   endif
breaksw

case struct_afm_check:
set exe = struct_afm_check
cat << theend > $def
5, '$file.inclmcopy',    'old',    'formatted', 0
6, '$file.outputstruct_afm_check','unknown','formatted', 0
20,'$file.struct',       'old',    'formatted', 0
21,'$file.struct_afm_check',       'unknown',    'formatted',       0
theend
breaksw

case clminter:
if ($updn == '' ) set updn = 'sum'
set exe = clminter
cat << theend > $def
6, '$file.outputclminter','unknown','formatted', 0
17,'$file.clm$updn',     'old',    'formatted',  0
18,'$file.clm${updn}_new',     'unknown','formatted',  0
20,'$file.struct',       'old',    'formatted', 0
21,'$file.struct_new',       'old',    'formatted',       0
theend
breaksw

case afminput:
set exe = afminput
cat << theend > $def
14, '$file.outputafminput',  'unknown','formatted', 0
15, '$file.inclmcopy_st',    'unknown','formatted', 0
16,'$file.clmup',            'old',    'formatted', 0
17,'$file.clmdn',            'old',    'formatted', 0
20,'$file.struct',           'old',    'formatted', 0
21,'$file.struct_supergroup','unknown','formatted', 0
theend
breaksw

case irrep:
set exe = irrep$para
set dopara = 1
cat << theend > $def
5, '$file.irrep$so$updn',      'unknown','formatted',0
6, '$file.outputir$so$updn',    'unknown','formatted',0
9, '${scratch}$file.vector$so$dnup',      'unknown',    'unformatted',9000
10,'${scratch}$file.vector$so$updn',      'old',    'unformatted',9000
20,'$file.struct',          'old',    'formatted',0
theend
breaksw

case xyz2struct:
if (-e $file.xyz ) then
xyz2struct <$file.xyz
cp fort.21 $file.struct_xyz
echo "File $file.struct_xyz has been created"
else
 echo "You need a xyz-file: $file.xyz " 
 exit(1)
endif
breaksw

case supercell:
set exe = supercell
cat << theend > $def

theend
breaksw

case plane:
if (-e plane.input) then
 plane <plane.input
else
 echo "You need a file:  plane.input (see top of $WIENROOT/SRC_trig/plane.f)" 
 exit(1)
endif
breaksw

case qtl:
set exe = $command$para
set dopara = 1
cat << theend > $def
 4,'$file.inso',   'unknown',    'formatted',0
 5,'$file.inq',   'old',    'formatted',0
 6,'$file.outputq$updn','unknown','formatted',0
 7,'$file.in1c',   'unknown',    'formatted',0
 9,'${scratch}$file.vector$so$updn', 'unknown','unformatted',9000
10,'${scratch}$file.vector$so$dnup', 'unknown','unformatted',9000
18,'$file.vsp$updn',       'old',    'formatted',0
19,'$file.vsp$dnup',       'unknown','formatted',0
20,'$file.struct',         'old',    'formatted',0
45,'$file.normsoup','unknown','formatted',0
46,'$file.normsodn','unknown','formatted',0
59,'$file.energy$so$dnup', 'old','formatted',0
60,'$file.energy$so$updn', 'unknown','formatted',0
theend
set natom=`head -5 $file.inq | tail -1 | cut -c0-10`
set iatom=1
while ($iatom <= $natom)
set help=`expr 30 + $iatom`
set cf=`expr 50 + $iatom`
     echo "$help,'$file.qtl$updn$iatom','unknown','formatted',0">>$def
     echo "$cf,'$file.cf$iatom','unknown','formatted',0">>$def
   @ iatom ++
end
breaksw

case struct2mol:
set exe = struct2mol
set molf='xtl'
cat << theend > $def
2, '$file.$molf',       'unknown',    'formatted',0
20,'$file.struct','old','formatted',0
theend
breaksw

case orbstr:
set exe = orbstr
cat << theend > $def
6, '$file.outputorbstr', 'unknown', 'formatted',0
8, '$file.inso',       'old',    'formatted',0
20,'$file.struct','old','formatted',0
21,'$file.struct_orbstr','unknown','formatted',0
theend
breaksw

case aim:
set exe = $command$cmplx
cat << theend > $def
5 ,'$file.inaim', 'old',    'formatted',0
6 ,'$file.outputaim',   'unknown','formatted',0
8 ,'$file.struct',    'old',    'formatted',0
9 ,'$file.clmsum',    'old',    'formatted',0
21,'$file.surf',       'unknown','formatted',0
22,'$file.crit',       'unknown','formatted',0
77,'$file.aim_surface_errors',       'unknown','formatted',0
theend
breaksw

case sgroup:
set settol1 = `echo "scale = 8; $settol / 10" | bc`
set exe = $command
$bin/$command -wi $file.struct  -set-TOL=$settol $file.outputsgroup 
$bin/$command -wi $file.struct  -set-TOL=$settol1 $file.outputsgroup1 
set test=`diff $file.outputsgroup $file.outputsgroup1 |wc`
if ( $test[1] != '0' ) then
    echo "Accuracy problem. Please run with different tolerance (x sgroup -settol $settol1)"
endif
set def = "-wi $file.struct -wo $file.struct_sgroup  -set-TOL=$settol"
breaksw

case nlo_KK:
set exe = $command
set def = "-w -i$file.innlo_KK -o$file.KK$updn $file.nlo$updn "
breaksw

case kzsurf:
set exe = $command
cat << theend > $def
6 ,'$file.outputkzsurf',   'unknown','formatted',0
20 ,'$file.struct',    'old',    'formatted',0
10 ,'$file.vector$updn',    'old',    'unformatted',0
theend
breaksw

case pairhess:
set exe = pairhess
cat << theend > $def
 6,'$file.outputpair','unknown','formatted',0
10,'$file.inpair',    'unknown',    'formatted',0
20,'$file.struct',    'old',    'formatted',0
21,'$file.inM_st',    'unknown',    'formatted',0
22,'$file.hess',    'unknown',    'formatted',0
theend
breaksw

case patchsymm:
set exe = patchsymm
cat << theend > $def
 6,'$file.outputpatch','unknown','formatted',0
20,'$file.struct',    'old',    'formatted',0
21,'$file.struct_new',    'unknown',    'formatted',0
theend
breaksw

case clmcreate:
set exe = clmcreate
cat << theend > $def
6, '$file.outputcreate','unknown','formatted', 0
17,'$file.clmsum',     'old',    'formatted',  0
18,'$file.clmsum_new',     'unknown','formatted',  0
20,'$file.struct',       'old',    'formatted', 0
21,'$file.increate',    'unknown','formatted', 0
theend
breaksw

default:
error	command $command does not exist
breaksw

endsw

if ($?deffile_only) then
  rm $running
  exit(0)
endif

if ( $para == para && $dopara) then
   if ( $updn == up ) then
      set exe = ($exe -up)
   else if ($updn == dn ) then
      set exe = ($exe -dn)
   endif
   if ($cmplx == c ) then
      set exe = ($exe -c)
   endif
   if ($so == so) then
      set exe = ($exe -so)
   endif
   if ($eece == eece) then
      set exe = ($exe -eece)
   endif
endif  	

echo $exe >>$running
$t $bin/$exe $def 
if($status != 0) then
  echo "error: command   $bin/$exe $def   failed"
  if ($?qtl || $?band || $?fermi || $?efg) then
    if( -e .oldin2 ) mv .oldin2 $file.in2$cmplx
  endif
  if (-f $running) rm $running
  exit(9)
endif

clear:
#cleanup
if ($?qtl || $?band || $?fermi || $?efg) then
  if( -e .oldin2 ) mv .oldin2 $file.in2$cmplx
endif
if (-f $running) rm $running
exit(0)

error:
rm $running

help:					#help exit 
cat << theend 

USAGE:	$0 PROGRAMNAME [flags] 

PURPOSE:runs WIEN executables:  afminput,aim,broadening,clmcopy,clminter,dipan,
      dstart,eosfit,eosfit6,filtvec,init_xspec,hex2rhomb,irrep,joint,kgen,kram,
      lapw0,lapw1,lapw2,lapw3,lapw5,lapw7,lapwdm,lapwso,lcore,lorentz,lstart,
      mini,mixer,nn,pairhess,plane,qtl,optic,optimize,orb,rhomb_in5,sgroup,
      spaghetti,struct_afm_check,sumpara,supercell,symmetry,symmetso,telnes2,
      tetra,txspec,xspec

FLAGS:
-f FILEHEAD ->	FILEHEAD for path of struct & input-files
-t/-T ->	suppress output of running time
-h/-H ->	help
-d    ->	create only the def-file 
-up   ->	runs up-spin
-dn   ->	runs dn-spin
-du   ->	runs up/dn-crossterm
-sc   ->	runs semicore calculation
-c    ->	complex calculation (no inversion symmetry present)	
-p    ->        run lapw0/1/2/so/dm/optic in parallel (needs .machines file)
-orb  ->	runs lapw1 with LDA+U/OP or B-ext correction
-it   ->	runs lapw1 with iterative diagonalization
-nohns->	runs lapw1 without HNS
-nmat_only->	runs lapw1 and yields only the matrixsize
-qtl  ->        calculates QTL in lapw2
-band ->        for bandstructures: uses *klist_band, sets QTL and ROOT (in2)
-fermi->        calculates lapw2 with FERMI switch
-efg  ->        calculates lapw2 with EFG switch
-so   ->	runs lapw2/optic/spaghetti with def-file for spin-orbit calc.
-fft  ->        runs dstart only up to case.in0_std creation
-sel  ->        use reduced vector file in lapw7
-settol 0.000x -> run sgroup with different tolerance
-sigma->        run lstart with case.inst_sigma (autogenerated) for diff.dens.
theend
if !($?command) then
  set command
endif

switch ($command)

case afminput:
echo 'x afminput'
breaksw
case aim:
echo 'x aim [-c]'
breaksw
case clmcopy:
echo 'x clmcopy'
breaksw
case clminter:
echo 'x clminter |-up/dn]'
breaksw
case dipan:
echo 'x dipan'
breaksw
case dstart:
echo 'x dstart [-up/-dn -c -fft]'
breaksw
case eosfit: 
echo 'x eosfit'
breaksw
case eosfit6 
echo 'x eosfit6'
breaksw
case init_xspec: 
echo 'x init_xspec [-up/-dn]'
breaksw
case irrep: 
echo 'x irrep [-so -up/-dn -p]'
breaksw
case joint: 
echo 'x joint [-up/-dn]'
breaksw
case kgen: 
echo 'x kgen [-so]     (-so switch only when you have a case.ksym file generated by initso_lapw)'
breaksw
case kram: 
echo 'x kram [-up/-dn]'
breaksw
case lapw0: 
echo 'x lapw0 [-p -eece]'
breaksw
case lapw1: 
echo 'x lapw1 [-c -up/-dn -it -p -nohns -orb -band -nmat_only]'
breaksw
case lapw2: 
echo 'x lapw2 [-c -up/-dn -p -so -qtl -fermi -efg -band -eece]'
breaksw
case lapw3: 
echo 'x lapw3 [-c]'
breaksw
case lapw5: 
echo 'x lapw5 [-c -up/-dn]'
breaksw
case filtvec: 
echo 'x filtvec [-c -up/-dn]'
breaksw
case lapw7: 
echo 'x lapw7 [-c -up/-dn -sel]'
breaksw
case lapwdm: 
echo 'x lapwdm [ -up -p -c]'
breaksw
case lapwso: 
echo 'x lapwso [ -up -p -c -orb]'
breaksw
case lcore: 
echo 'x lcore [-up/-dn]'
breaksw
case lorentz: 
echo 'x lorentz [-up/-dn]'
breaksw
case lstart: 
echo 'x lstart [-sigma]'
breaksw
case mini: 
echo 'x mini'
breaksw
case mixer: 
echo 'x mixer'
breaksw
case nn: 
echo 'x nn'
breaksw
case optic: 
echo 'x optic [-c -up/-dn -so -p]'
breaksw
case optimize: 
echo 'x optimize'
breaksw
case orb: 
echo 'x orb [ -up/-dn/-du ]'
breaksw
case pairhess: 
echo 'x pairhess'
breaksw
case plane: 
echo 'x plane  ( needs plane.input, see   head $WIENROOT/SRC_trig/plane.f )'
breaksw
case qtl: 
echo 'x qtl [ -up/-dn -so -p]'
breaksw
case sgroup: 
echo 'x sgroup [-settol 0.000x]  (try 0.001 - 0.0000001)'
breaksw
case spaghetti: 
echo 'x spaghetti [-up/-dn -so -p]'
breaksw
case struct_afm_check:
echo 'x struct_afm_check '
breaksw
case sumpara: 
echo 'x sumpara -d [-up/-dn/-dui -eece]     and than '
echo 'sumpara [up/dn]sumpara.def #_of_processors'
breaksw
case supercell: 
echo 'x supercell'
breaksw
case symmetry: 
echo 'x symmetry'
breaksw
case symmetso: 
echo 'x symmetso [-c]'
breaksw
case tetra: 
echo 'x tetra [-up/-dn]'
breaksw
case txspec: 
echo 'x txspec [-up/-dn]'
breaksw
case xspec: 
echo 'x xspec [-up/-dn]'
breaksw
case telnes2:
echo 'x telnes2 [-up/-dn]'
breaksw

default:
echo 'USE: x -h PROGRAMNAME   for valid flags for a specific program'
breaksw

endsw
if (-e $running) rm $running

exit(1)
