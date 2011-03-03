#!/usr/bin/perl


require "../../libs/w2web.pl";
require "../../libs/struct.pl";
$debug=0;
&GetInput;
&GetSession;
&StructRead;

$show=0;

$submit = "<i>View only mode</i>  --&gt;<a href=\"/util/structstart.pl?SID=$SID\">edit STRUCT file</A><br>";
$bgcolor="$red";
$comment="StructView $DIR";
$readonly="readonly";#if structedit is true, we have to set this to ""
$disabled="disabled";#if structedit is true, we have to set this to ""

if ($structedit) {
	$show=1;
	$bgcolor="$green";
	$readonly="";
	$disabled="";
	$comment="StructEdit $DIR";
	$submit =  <<__STOP__;
<i>You have to click "Save Structure" for changes to take effect!</i>
<br>
<INPUT TYPE=SUBMIT VALUE="Save Structure">
__STOP__
}

$lat = "<SELECT NAME=\"s_lattice\" SIZE=9>";
$s_spacegr=~ s/ /_/;
foreach $v (@lattype) {
	if ( $s_spacegr =~ / */ && $s_lattice =~ /$v/ ) {
		$lat .= "<OPTION VALUE=$v SELECTED>$v\n";
	} elsif ( $s_spacegr =~ /$v/ ) {
		$lat .= "<OPTION VALUE=$v SELECTED>$v\n";
	} else {
		$lat .= "<OPTION VALUE=$v>$v\n";
	}
}
$lat .= "</SELECT>";

$rela = "<SELECT NAME=\"s_rels\">";
$rel = "";
$nrel = "";

if ($s_rels =~ /RELA/) {$rel = "SELECTED" ;}
if ($s_rels =~ /NREL/) {$nrel = "SELECTED" ;}

$rela .= <<__STOP__;
<OPTION VALUE=RELA $rel>relativistic
<OPTION VALUE=NREL $nrel>non relativistic
</SELECT>
__STOP__

$atom = 0;

if ($FORM{'split'}) { 
	($atom, $eq) = split(/,/, $FORM{'split'}, 2);
	$comment = "split atom/position mode: ";
	$s_nsym=0;
	$show=0;
		# serious business ahead: we have to add an atom
		$s_nato++;
		$show=0;
		$s_name[$s_nato]=$s_name[$atom];
		$s_nameadd[$s_nato]=$s_nameadd[$atom];
		$s_zz[$s_nato]=$s_zz[$atom];
		$s_mult[$s_nato]=1;
		$s_x[$s_nato][1]=$s_x[$atom][$eq];
		$s_y[$s_nato][1]=$s_y[$atom][$eq];
		$s_z[$s_nato][1]=$s_z[$atom][$eq];
		$s_cubic[$s_nato]=0;
		$s_npt[$s_nato]=$s_npt[$atom];
		$s_ro[$s_nato]=$s_ro[$atom];
		$s_rmt[$s_nato]=$s_rmt[$atom];
		
		$s_rot11[$s_nato]= $s_rot11[$atom];
		$s_rot12[$s_nato]= $s_rot12[$atom];
		$s_rot13[$s_nato]= $s_rot13[$atom];
		$s_rot21[$s_nato]= $s_rot21[$atom];
		$s_rot22[$s_nato]= $s_rot22[$atom];
		$s_rot23[$s_nato]= $s_rot23[$atom];
		$s_rot31[$s_nato]= $s_rot31[$atom];
		$s_rot32[$s_nato]= $s_rot32[$atom];
		$s_rot33[$s_nato]= $s_rot33[$atom];

		# now  an equivalent position
		$comment .= "position $atom,$eq is now inequivalent atom $s_nato";
		for ($i = $eq; $i < $s_mult[$atom]; $i++) {
			my $next=$i;
			$next++;
			$s_x[$atom][$eq]=$s_x[$atom][$next];
			$s_y[$atom][$eq]=$s_y[$atom][$next];
			$s_z[$atom][$eq]=$s_z[$atom][$next];
		}
		$s_mult[$atom]--;
}

if ($FORM{'add'}) { 
	# we have to add atoms
	($atom, $eq) = split(/,/, $FORM{'add'}, 2);
	$comment = "add atom/position mode";
	$s_nsym=0;
	$show=0;
	if ($atom > $s_nato) {

		# serious business ahead: we have to add an atom
		$show=0;

		$s_name[$atom]="";
		$s_nameadd[$atom]="";
		$s_zz[$atom]="";
		$s_mult[$atom]=1;
		$s_x[$atom][1]="";
		$s_y[$atom][1]="";
		$s_z[$atom][1]="";
		$s_cubic[$atom]=0;
		$s_npt[$atom]=781;
		$s_ro[$atom]=0.0005;
		$s_rmt[$atom]=2.0;
		
		$s_rot11[$atom]=0.0;
		$s_rot12[$atom]=0.0;
		$s_rot13[$atom]=0.0;
		$s_rot21[$atom]=0.0;
		$s_rot22[$atom]=0.0;
		$s_rot23[$atom]=0.0;
		$s_rot31[$atom]=0.0;
		$s_rot32[$atom]=0.0;
		$s_rot33[$atom]=0.0;

		$s_nato++;
	} else {
		$s_x[$atom][$eq]="";
		$s_y[$atom][$eq]="";
		$s_z[$atom][$eq]="";
		$s_mult[$atom]++;
	}
}

if ($FORM{'del'}) { 
	# we have to add atoms
	($atom, $eq) = split(/,/, $FORM{'del'}, 2);
	$s_nsym=0;
	$show=0;

	if ($eq ==0) {
		$comment = "atom $atom removed" ;
		# gheee remove a whole site ...
		for ($i = $atom; $i < $s_nato; $i++) {
			my $ni = $i;
			$ni++;
			for ($j = 1; $j <= $s_mult[$ni]; $j++) {
				$s_x[$i][$j]=$s_x[$ni][$j];
				$s_y[$i][$j]=$s_y[$ni][$j];
				$s_z[$i][$j]=$s_z[$ni][$j];
			}
			# don't forget site specific parameters
			$s_cubic[$i] = $s_cubic[$ni];
			$s_mult[$i] = $s_mult[$ni];
			$s_isplit[$i] = $s_isplit[$ni];
			$s_name[$i] = $s_name[$ni];
			$s_nameadd[$i] = $s_nameadd[$ni];
			$s_zz[$i] = $s_zz[$ni];
			$s_npt[$i] = $s_npt[$ni];
			$s_ro[$i] = $s_ro[$ni];
			$s_rmt[$i] = $s_rmt[$ni];
			$s_rot11[$i] = $s_rot11[$ni];
			$s_rot12[$i] = $s_rot12[$ni];
			$s_rot13[$i] = $s_rot13[$ni];
			$s_rot21[$i] = $s_rot21[$ni];
			$s_rot22[$i] = $s_rot22[$ni];
			$s_rot23[$i] = $s_rot23[$ni];
			$s_rot31[$i] = $s_rot31[$ni];
			$s_rot32[$i] = $s_rot32[$ni];
			$s_rot33[$i] = $s_rot33[$ni];
		}	
		$s_nato--;


	} else {
		# remove an equivalent position
		$comment = "position $atom,$eq removed";
		for ($i = $eq; $i < $s_mult[$atom]; $i++) {
			my $next=$i;
			$next++;
			$s_x[$atom][$eq]=$s_x[$atom][$next];
			$s_y[$atom][$eq]=$s_y[$atom][$next];
			$s_z[$atom][$eq]=$s_z[$atom][$next];
		}
		$s_mult[$atom]--;
	}	
}

$OUT .=  <<__STOP__;

<H2>StructGen<font size=-2><sup>TM</sup></font> </H2>
<table bgcolor=$bgcolor><tr><td>
<FORM ACTION="/util/structsave.pl" METHOD=POST>
$submit 
<br>
<br>
__STOP__

&PassHiddenParms;

$OUT .=  <<__STOP__;
<b>Title:&nbsp; </b> <INPUT NAME="s_title" VALUE="$s_title" SIZE=40 $readonly>
<br>
<INPUT TYPE=HIDDEN NAME="s_rels" VALUE="$s_rels">
<table><tr><td>
<b>Lattice:</b>
<br>
__STOP__
#Treatment of core: 
#$rela
#<br>

if ($generateeq) {
	$OUT .= "Spacegroup: $s_spacegr<br>$lat<br><br>Splitting of equivalent positions not available.<br> To split you must select a lattice type";
#<INPUT NAME=\"s_spacegr\" VALUE=\"$s_spacegr\">";
} else {
$OUT .= "Type: $s_lattice<br>$lat";
}

$bohr="";
$ang="";
$ang="SELECTED" if ($unit =~ /ang/);
$bohr="SELECTED" if ($unit =~ /bohr/);

$OUT .=  <<__STOP__;
</td><td><A HREF="http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-wp-list" TARGET=spgrp>Spacegroups from<br>Bilbao Cryst Server</A></td></tr></table>

<b>Lattice parameters </b>in 
<SELECT NAME="unit">
<OPTION VALUE="ang" $ang>&Aring;
<OPTION VALUE="bohr" $bohr>bohr
</SELECT>
<br>
$indent a=<INPUT NAME="s_a" VALUE="$s_a" SIZE=10 $readonly>
b=<INPUT NAME="s_b" VaLUE="$s_b" SIZE=10 $readonly>
c=<INPUT NAME="s_c" VALUE="$s_c" SIZE=10 $readonly>
<br>$indent &alpha;=<INPUT NAME="s_aa" VALUE="$s_aa" SIZE=10 $readonly>
&beta;=<INPUT NAME="s_bb" VaLUE="$s_bb" SIZE=10 $readonly>
&gamma;=<INPUT NAME="s_cc" VALUE="$s_cc" SIZE=10 $readonly>
<br>
<br>
<b>Inequivalent Atoms: $s_nato</b>
<INPUT TYPE=HIDDEN NAME="s_nato" VALUE="$s_nato">
<INPUT TYPE=HIDDEN NAME="generateeq" VALUE="$generateeq">
__STOP__


	
for ($i = 1; $i <= $s_nato; $i++) {
	$remove= "";
	if($FORM{'split'}) {
	} else {
		if($show && !$generateeq) {
			$remove = <<__STOP__;
<A HREF="structgen.pl?SID=$SID&del=$i,1">remove</A>
<A HREF="structgen.pl?SID=$SID&split=$i,1">split</A> 
__STOP__
		} else {
			$remove = "&lt;-- edit only this position!" if $structedit ;
		}
	}
	$OUT .= <<__STOP__;
<br><b>Atom $i</b>: 
<INPUT TYPE=HIDDEN NAME="s_mult[$i]" VALUE=$s_mult[$i]>
<INPUT SIZE=3 NAME="s_name[$i]" VALUE="$s_name[$i]" MAXLENGTH=2 $readonly>
<INPUT SIZE=8 NAME="s_nameadd[$i]" VALUE="$s_nameadd[$i]" $readonly>
<INPUT TYPE=HIDDEN NAME="s_isplit[$i]" VALUE="$s_isplit[$i]">
Z=<INPUT SIZE=4 NAME="s_zz[$i]" VALUE="$s_zz[$i]" $readonly>
RMT=<INPUT SIZE=5 NAME="s_rmt[$i]" VALUE="$s_rmt[$i]" $readonly>
__STOP__
    if ($show){
	$OUT .= <<__STOP__;
$indent <A HREF="/util/structgen.pl?SID=$SID&del=$i,0">remove atom</A>
__STOP__
	}
	$OUT .= <<__STOP__;

<INPUT TYPE=HIDDEN NAME="s_npt[$i]" VALUE="$s_npt[$i]"> 
<INPUT TYPE=HIDDEN NAME="s_ro[$i]" VALUE="$s_ro[$i]"> 
<INPUT TYPE=HIDDEN NAME="s_cubic[$i]" VALUE="$s_cubic[$i]"> 
__STOP__
#<br>
#$indent NPT=<INPUT SIZE=5 NAME="s_npt[$i]" VALUE="$s_npt[$i]" $readonly>
#R0=<INPUT SIZE=7 NAME="s_ro[$i]" VALUE="$s_ro[$i]" $readonly>
#Non-cubic <INPUT TYPE=CHECKBOX NAME="s_cubic[$i]" $readonly
#	$OUT.="CHECKED" if ($s_cubic[$i]);
#>
	$OUT .= <<__STOP__;
<br>
__STOP__
	if ($s_lattice =~ /R/) {
		$OUT .= "$indent positions must be specified in rhombohedral coordinates!<br>";
	}
	$OUT .= <<__STOP__;
$indent Pos 1:
x=<INPUT NAME="s_x[$i][1]" VALUE="$s_x[$i][1]" SIZE=10 $readonly>
y=<INPUT NAME="s_y[$i][1]" VALUE="$s_y[$i][1]" SIZE=10 $readonly>
z=<INPUT NAME="s_z[$i][1]" VALUE="$s_z[$i][1]" SIZE=10 $readonly>
$remove
<br>
__STOP__

	if ($s_mult[$i] > 1) {
		for ($j = 2; $j <= $s_mult[$i]; $j++) {
			$remove= "";
			$mysplit="split" if $show;
			$mysplit="";
			if($show && !$generateeq) {
				$remove = <<__STOP__;
<A HREF="structgen.pl?SID=$SID&del=$i,$j">remove</A>
__STOP__
				$mysplit = <<__STOP__;
<A HREF="structgen.pl?SID=$SID&split=$i,$j">split</A>
__STOP__
			}
			$OUT .= <<__STOP__;
$indent Pos $j:
x=<INPUT NAME="s_x[$i][$j]" VALUE="$s_x[$i][$j]" SIZE=10 $readonly>
y=<INPUT NAME="s_y[$i][$j]" VALUE="$s_y[$i][$j]" SIZE=10 $readonly>
z=<INPUT NAME="s_z[$i][$j]" VALUE="$s_z[$i][$j]" SIZE=10 $readonly>
$remove
$mysplit
<br>
__STOP__
		}
	}
	# local rotation matix
	$OUT .= <<__STOP__;
<INPUT NAME="s_rot11[$i]" VALUE="$s_rot11[$i]" SIZE=10 TYPE=HIDDEN>
<INPUT NAME="s_rot12[$i]" VALUE="$s_rot12[$i]" SIZE=10 TYPE=HIDDEN>
<INPUT NAME="s_rot13[$i]" VALUE="$s_rot13[$i]" SIZE=10 TYPE=HIDDEN>
<INPUT NAME="s_rot21[$i]" VALUE="$s_rot21[$i]" SIZE=10 TYPE=HIDDEN>
<INPUT NAME="s_rot22[$i]" VALUE="$s_rot22[$i]" SIZE=10 TYPE=HIDDEN>
<INPUT NAME="s_rot23[$i]" VALUE="$s_rot23[$i]" SIZE=10 TYPE=HIDDEN>
<INPUT NAME="s_rot31[$i]" VALUE="$s_rot31[$i]" SIZE=10 TYPE=HIDDEN>
<INPUT NAME="s_rot32[$i]" VALUE="$s_rot32[$i]" SIZE=10 TYPE=HIDDEN>
<INPUT NAME="s_rot33[$i]" VALUE="$s_rot33[$i]" SIZE=10 TYPE=HIDDEN>
__STOP__

	if($show && !$generateeq) {
		$OUT .= <<__STOP__;
$indent <A HREF="/util/structgen.pl?SID=$SID&add=$i,$j">add position</A>
<br>
__STOP__
	}

	
}

if ($show) {
	$OUT .= "<br><A HREF=\"/util/structgen.pl?SID=$SID&add=$i,1\">add an atom</A><br>";
}

# symmetry operations
$OUT .= "<br>";
$OUT .= "Number of symmetry operations: \n";

if($s_nsym > 0) {
	$OUT .= <<__STOP__;
<SELECT NAME="s_nsym">
<OPTION VALUE="$s_nsym">$s_nsym
<OPTION VALUE="0">generate
</SELECT>
__STOP__
} else {
	$OUT .= "generate";
}
$OUT .= "<br>";

for ($i = 1; $i <= $s_nsym; $i++) {
	$OUT .= <<__STOP__;
<INPUT TYPE=HIDDEN NAME="s_s11[$i]" VALUE="$s_s11[$i]">
<INPUT TYPE=HIDDEN NAME="s_s12[$i]" VALUE="$s_s12[$i]">
<INPUT TYPE=HIDDEN NAME="s_s13[$i]" VALUE="$s_s13[$i]">
<INPUT TYPE=HIDDEN NAME="s_t1[$i]" VALUE="$s_t1[$i]">
<INPUT TYPE=HIDDEN NAME="s_s21[$i]" VALUE="$s_s21[$i]">
<INPUT TYPE=HIDDEN NAME="s_s22[$i]" VALUE="$s_s22[$i]">
<INPUT TYPE=HIDDEN NAME="s_s23[$i]" VALUE="$s_s23[$i]">
<INPUT TYPE=HIDDEN NAME="s_t2[$i]" VALUE="$s_t2[$i]">
<INPUT TYPE=HIDDEN NAME="s_s31[$i]" VALUE="$s_s31[$i]">
<INPUT TYPE=HIDDEN NAME="s_s32[$i]" VALUE="$s_s32[$i]">
<INPUT TYPE=HIDDEN NAME="s_s33[$i]" VALUE="$s_s33[$i]">
<INPUT TYPE=HIDDEN NAME="s_t3[$i]" VALUE="$s_t3[$i]">
__STOP__
}
$OUT .= "<br>$submit</FORM>";
$OUT .= "</td></tr></table>";

&PrintPage("Context",$OUT);


