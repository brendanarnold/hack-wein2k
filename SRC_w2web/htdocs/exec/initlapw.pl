#!/usr/bin/perl

require "../../libs/w2web.pl";
$phase=0;
$doinst=0;
$endphase=0;
#$complex="";
#$c="";
#$spinpol="";
$spin="up";
&GetInput;
&GetSession;
$restart_init = "";
$addtext="";
$endtext="";
$nextinteractive=1;

if (!$phase) {
    $OUT .= "Phase not found" if ($debug);
    if ( -e "$DIR/.initphase" ) {
	$oldphase = qx(cat $DIR/.initphase);
	$tmp = $oldphase -1;
	$restart_init .= "<FORM ACTION=/exec/initlapw.pl><INPUT TYPE=HIDDEN NAME=SID VALUE=$SID><INPUT TYPE=HIDDEN NAME=phase VALUE=$tmp><INPUT TYPE=SUBMIT VALUE=\"Restart with phase $oldphase?\"></FORM>";
    }
}

$inactivecolor=$gray;
if ($phase == 21 || $phase == 25  ) {
    $inactivecolor=$gray;
    $endtext="<h2>Initialization done</h2><p class=\"info\"><A HREF=/exec/scf.pl?SID=$SID>Continue with run SCF</A></p>";
    $complex="CHECKED" if ( -e "$DIR/$CASE.in1c" );
    &SaveSession;
}


sub mark {
    $bgcolor="taskactive";
}
sub reset {
    $bgcolor="task";
    $nextinteractive=1;
}

&reset;


sub Execute {
    $next    = "initlapw";
    $ia = "";
    if ($interactive) {
	$ia="<INPUT TYPE=checkbox NAME=nextinteractive CHECKED >interactively"
	} else {
	    $ia="<INPUT TYPE=hidden NAME=nextinteractive VALUE=1>"
	    }
    $nexturl = "/exec/initlapw.pl?SID=$SID&phase=$myphase";
    $OUT .=  <<__STOP__;
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
&PassHiddenParms;
    $OUT .=  <<__STOP__;
<TR><TD class="$bgcolor">
<INPUT TYPE=HIDDEN NAME=precmd VALUE="$precmd">
<INPUT TYPE=hidden NAME=prog VALUE=$prog>
<INPUT TYPE=hidden NAME=c VALUE=$c>
<INPUT TYPE=hidden NAME=spinpol VALUE=$spinpol>
<INPUT TYPE=hidden NAME=spin VALUE=$spin>
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
<INPUT TYPE=SUBMIT VALUE="$title" >
$ia</FORM>$addtext</TD></TR>
__STOP__
     $addtest="";
    &reset;
}

sub Edit {
    $next = "initlapw";
    $dum=$myphase;
    if ($endphase > 0) {
	$dum = $endphase;
    }
    $nexturl = "/exec/initlapw.pl?SID=$SID&phase=$dum";
    $OUT .=  <<__STOP__;
<TR><TD class="$bgcolor">
<FORM ACTION=/util/edit.pl METHOD=POST>
__STOP__
    &PassHiddenParms;
    $OUT .=  <<__STOP__;
<INPUT TYPE=hidden NAME=file VALUE="$file">
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
<INPUT TYPE=SUBMIT VALUE="$title" >
</FORM>$addtext</TD></TR>
__STOP__
    $addtext="";
    &reset;
}

$phase++;
$OUT .=  <<__STOP__;
<h2>Initialize calculation <small>(phase $phase)</small></h2>
$restart_init
<table border=0 cellspacing=0 cellpadding=1 class="$bgcolor"><tr valign=top><td>
<table cellpadding=1 cellspacing=2>
__STOP__

$myphase = 1;
&mark if ($phase == $myphase);
$title = "x nn";
$precmd  = "x";
$prog    = "nn";
$next    = "initlapw";
&Execute();



$myphase = 2;
&mark if ($phase == $myphase);
$title = "view outputnn";
$file    = "$DIR/$CASE.outputnn";
&Edit();


$myphase = 3;
$bb = qx(grep WARNING $DIR/$CASE.outputnn | wc -l );
if ($bb > 0 && $phase == $myphase) {
	&mark if ($phase == $myphase);
	$title = "choose";
	$OUT .= "<tr><td class=\"$bgcolor\">";
	if ($doit) {
		$dum = qx(cp $DIR/$CASE.struct $DIR/$CASE.struct_init);
		$dum = qx(cp $DIR/$CASE.struct_nn $DIR/$CASE.struct);
		$dum = qx(cd $DIR; echo y | instgen);
		$OUT .= <<__STOP__;
<b>$CASE.struct_nn copied to $CASE.struct<br> 
   old struct-file saved as $CASE.struct_init<br>
   $CASE.inst updated</b>
<br>
<FORM ACTION=/util/structgen.pl>
<INPUT TYPE=HIDDEN NAME=SID VALUE=$SID>
<INPUT TYPE=SUBMIT VALUE="Start StructGen ?">
</FORM>
__STOP__
	} else {
		$OUT .= <<__STOP__;
Use new struct-file? 
<ul>
<FORM ACTION=/exec/initlapw.pl>
<INPUT TYPE=HIDDEN NAME=SID VALUE=$SID>
<INPUT TYPE=HIDDEN NAME=phase VALUE=2>
<INPUT TYPE=HIDDEN NAME=doit VALUE=1>
<INPUT TYPE=SUBMIT VALUE="Yes">
</FORM>
<FORM ACTION=/exec/initlapw.pl>
<INPUT TYPE=HIDDEN NAME=SID VALUE=$SID>
<INPUT TYPE=HIDDEN NAME=phase VALUE=$myphase>
<INPUT TYPE=SUBMIT VALUE="No">
</FORM>
</UL>
__STOP__
	}
	$OUT .= "</td></tr>";
	&reset;
}


$myphase = 4;
&mark if ($phase == $myphase);
&mark if ($phase == 3 && $bb == 0);
$title = "x sgroup";
$precmd  = "x";
$prog    = "sgroup";
$next    = "initlapw";
&Execute();

$myphase=5;
&mark if ($phase == $myphase);
$title = "view outputsgroup";
$file    = "$DIR/$CASE.outputsgroup";
$sg1 = qx(grep '\!\!' $DIR/$CASE.outputsgroup  );
$sg = qx(grep 'space group' $DIR/$CASE.outputsgroup | grep -v NOTE );
$sg =~ s/.*://;
$addtext="$sg1<br>sgroup found: $sg" if ($phase == $myphase || $phase == 6);
&Edit();

$myphase = 6;
if ($phase == $myphase) {
$docopy=$myphase-1;
&mark if ($phase == $myphase);
$title = "choose";
if ($doit) {
	$OUT .= "<tr><td class=\"taskoption\">";
	$dum = qx(cp $DIR/$CASE.struct $DIR/$CASE.struct_init);
	$dum = qx(cp $DIR/$CASE.struct_sgroup $DIR/$CASE.struct );
	$dum = qx(cd $DIR; echo y | instgen);

	$OUT .= <<__STOP__;
<b>$CASE.struct_nn copied to $CASE.struct<br> 
   old struct-file saved as $CASE.struct_init<br>
   $CASE.inst updated</b>
<br>
<FORM ACTION=/util/structgen.pl>
<INPUT TYPE=HIDDEN NAME=SID VALUE=$SID>
<INPUT TYPE=SUBMIT VALUE="Start StructGen ?"> optional
</FORM>
__STOP__
	$phase++;
	} else {
		$OUT .= <<__STOP__;
<tr><td class="$bgcolor">
Use struct-file generated by sgroup? (Usually NO, unless WARNINGS appeared above) 
<ul>
<FORM ACTION=/exec/initlapw.pl>
<INPUT TYPE=HIDDEN NAME=SID VALUE=$SID>
<INPUT TYPE=HIDDEN NAME=phase VALUE=$docopy>
<INPUT TYPE=HIDDEN NAME=doit VALUE=1>
<INPUT TYPE=SUBMIT VALUE="Yes">
</FORM>
<FORM ACTION=/exec/initlapw.pl>
<INPUT TYPE=HIDDEN NAME=SID VALUE=$SID>
<INPUT TYPE=HIDDEN NAME=phase VALUE=$myphase>
<INPUT TYPE=SUBMIT VALUE="No">
</FORM>
</UL>
__STOP__
}
$OUT .= "</td></tr>";
&reset;
}

$myphase = 7;
&mark if ($phase == $myphase);
$title = "x symmetry";
$precmd  = "x";
$prog    = "symmetry";
$next    = "initlapw";
&Execute();

$myphase = 8;
&mark if ($phase == $myphase);
$title = "copy struct_st";
$addtext = "and view outputs";
$file    = "$DIR/$CASE.outputs";
&Edit();

$myphase = 9;
$bb = qx(grep SHIFTED $DIR/$CASE.outputs | wc -w );
if ($bb == 9) {
  &mark if ($phase == $myphase);
  $title = "choose";
  $OUT .= <<__STOP__;
<tr><td class="$bgcolor">
<b>Warning!</b><br>
You must move the origin 
<br>
(see file $CASE.outputs)
<FORM ACTION=/util/structgen.pl>
<INPUT TYPE=HIDDEN NAME=SID VALUE=$SID>
<INPUT TYPE=SUBMIT VALUE="Start StructGen ?">
</FORM>
__STOP__
	$OUT .= "</td></tr>";
	&reset;
}

$myphase = 10;
if ($phase == $myphase || ($phase == 9 && $bb!=9) ) {
&mark;
$umps = qx(cp $DIR/$CASE.struct_st  $DIR/$CASE.struct);
}
$title = "x lstart";
$precmd  = "x";
$prog    = "lstart";
$next    = "initlapw";
&Execute();



$myphase=11;
&mark if ($phase == $myphase);
$title = "view outputst";
$file    = "$DIR/$CASE.outputst";
&Edit();

$OUT .= "</table></td><td><table cellpadding=1 cellspacing=0>";

$myphase=12;
$bb = qx(grep IZ $DIR/$CASE.outputst | wc -w);
if ($bb > 2 || $doinst) {
	&mark if ($phase == $myphase);
	$title = "edit inst";
	$file    = "$DIR/$CASE.inst";
	if($doinst) {
		$bb=3; # to cheat next if statemend :-)
	} else {
		$addtext="error: $CASE.inst not consistent with Z<br>edit $CASE.inst and rerun lstart afterwards or change Z in StructGen!";
	}
	$endphase=9 if ($phase == $myphase);
	&Edit();
	$OUT.=<<__STOP__;
<tr><td class="$bgcolor">
<FORM ACTION=/util/structgen.pl>
<INPUT TYPE=HIDDEN NAME=SID VALUE=$SID>
<INPUT TYPE=SUBMIT VALUE="Start StructGen ?">
</FORM>
</td></tr>
__STOP__
}


$myphase = 13;
&mark if ($phase == $myphase);
&mark if ($phase == 12 && $bb <= 2);
$title = "check $CASE.in1_st";
$addtext="check RKmax (usually 5.0-9.0)";
$file    = "$DIR/$CASE.in1_st";
&Edit();

$myphase = 14;
if ($phase ==$myphase ) {
	&mark;
	$umps=qx(cd $DIR; cat $CASE.in2_ls $CASE.in2_sy > $CASE.in2_st);
}
$title = "check $CASE.in2_st";
$file    = "$DIR/$CASE.in2_st";
$addtext="set LM's, GMAX and Fermi-Energy method";
&Edit();


$myphase = 15;
&mark if ($phase == $myphase);
$title = "check $CASE.inm_st";
$file    = "$DIR/$CASE.inm_st";
$addtext="Reduce mixing for localized magnetic systems";
&Edit();


$myphase = 16;
&mark if ($phase == $myphase);
$OUT .= "<tr><td class=\"$bgcolor\">";
$tmp=$myphase+1;
if ($doit && $phase == $tmp) {
	$myphase++;
	$dum = qx(cp $DIR/$CASE.in0_st $DIR/$CASE.in0);
	$dum = qx(cp $DIR/$CASE.in1_st $DIR/$CASE.in1);
	$dum = qx(cp $DIR/$CASE.in2_st $DIR/$CASE.in2);
	$dum = qx(cp $DIR/$CASE.inc_st $DIR/$CASE.inc);
	$dum = qx(cp $DIR/$CASE.inm_st $DIR/$CASE.inm);
	$mout .= "in0, in1, in2, inc and inm files generated";

	$bb = qx(grep INVERSION $DIR/$CASE.outputs | wc -w);
	if ($bb == 7) {
		# no inversion symmetry
		$dum = qx(mv $DIR/$CASE.in1 $DIR/$CASE.in1c);
		$dum = qx(mv $DIR/$CASE.in2 $DIR/$CASE.in2c);
		$mout .= "in0, in1c, in2c, inc and inm files generated";
	}
	$OUT .= $mout;
} else {
	$OUT .= <<__STOP__;
<FORM ACTION=/exec/initlapw.pl>
<INPUT TYPE=HIDDEN NAME=SID VALUE=$SID>
<INPUT TYPE=HIDDEN NAME=phase VALUE=$myphase>
<INPUT TYPE=HIDDEN NAME=doit VALUE=1>
<INPUT TYPE=SUBMIT VALUE="Prepare input files" >
</FORM>
__STOP__
}
&reset;

$myphase = 17;
&mark if ($phase == $myphase);
$title = "x kgen";
$precmd  = "x";
$prog    = "kgen";
$next    = "initlapw";
&Execute();


$myphase=18;
&mark if ($phase == $myphase);
$title = "view klist";
$file    = "$DIR/$CASE.klist";
&Edit();

$myphase = 19;
&mark if ($phase == $myphase);
if ( -e "$DIR/$CASE.in1c" ) {
	$c="on";
	$addtext="<br>complex selected";
}
$spinpol="";
$spin="";
$title = "x dstart";
$precmd  = "x";
$prog    = "dstart";
$interactive = 1;
$next    = "initlapw";
&Execute();


$OUT .= "</table></td><td><table cellpadding=1 cellspacing=0>";

$myphase=20;
&mark if ($phase == $myphase);
$title = "view $CASE.outputd and cp $CASE.in0_std $CASE.in0";
$addtext = "check if gmax>gmin";
$file    = "$DIR/$CASE.outputd";
&Edit();
$dum = qx(cp $DIR/$CASE.in0_std $DIR/$CASE.in0);



$myphase=21;

if ($phase < 22) {
&mark if ($phase == $myphase);
$title = "choose";
$OUT .= <<__STOP__;
<tr><td class="$bgcolor">
Perform spin-polarized calc.?
<ul>
<FORM ACTION=/exec/initlapw.pl>
<INPUT TYPE=HIDDEN NAME=SID VALUE=$SID>
<INPUT TYPE=HIDDEN NAME=phase VALUE=21>
<INPUT TYPE=SUBMIT VALUE="No" >
</FORM>
<FORM ACTION=/exec/initlapw.pl>
<INPUT TYPE=HIDDEN NAME=SID VALUE=$SID>
<INPUT TYPE=HIDDEN NAME=phase VALUE=22>
<INPUT TYPE=SUBMIT VALUE="Yes" >
</FORM>
</UL>
__STOP__
$OUT .= "</td></tr>";
&reset;
}


$myphase=22; 
if ($phase==$myphase) {
	$spinpol="";
	&SaveSession;
#	$prefspace="single";
#	&SavePrefs();
#	$prefspace="scf";
#	&SavePrefs();
}
# dummy for "we are done" if non-sp calc.
#


# ------ spin polarized ---

if ($phase > 22) {
$myphase = 23;
&mark if ($phase == $myphase);
if ( -e "$DIR/$CASE.in1c" ) {
  $c="on";
  $addtext="<br>complex selected";
}

$spinpol="on";
&SaveSession;
# bloody SaveSession corrupts my $spinpol
# easy way out: just set it again to "on" *ggggg
#
#$prefspace="single";
#SavePrefs();
#$prefspace="scf";
#SavePrefs();

$spinpol="on";
$spin="up";
$title = "x dstart -$spin";
$precmd  = "x";
$prog    = "dstart -up";
$interactive = 1;
$next    = "initlapw";
&Execute();
&reset;

$myphase = 24;
&mark if ($phase == $myphase);
if ( -e "$DIR/$CASE.in1c" ) {
  $c="on";
  $addtext="<br>complex selected";
}
$spinpol="on";
$spin="dn";
$title = "x dstart -$spin";
$precmd  = "x";
$prog    = "dstart";
$interactive = 1;
$next    = "initlapw";
&Execute();
&reset;

$myphase=25;
&mark if ($phase == $myphase);
$title = "choose";
$OUT .= <<__STOP__;
<tr><td class="$bgcolor">
Perform antiferromagnetic calc.? 
<ul>
<FORM ACTION=/exec/initlapw.pl>
<INPUT TYPE=HIDDEN NAME=SID VALUE=$SID>
<INPUT TYPE=HIDDEN NAME=phase VALUE=26>
<INPUT TYPE=SUBMIT VALUE="Yes">
</FORM>
<FORM ACTION=/exec/initlapw.pl>
<INPUT TYPE=HIDDEN NAME=SID VALUE=$SID>
<INPUT TYPE=HIDDEN NAME=phase VALUE=25>
<INPUT TYPE=SUBMIT VALUE="No">
</FORM>
</UL>
__STOP__
$OUT .= "</td></tr>";
&reset;


$myphase=26;
$afm="";
&SaveSession;
#$prefspace="scf";
#&SavePrefs();
$afm="";
# dummy: no afm calc

if ($phase > $myphase) {

$myphase=27;
if ($phase == $myphase) {
&mark if ($phase == $myphase);
$afm="on";
&SaveSession;
#$prefspace="scf";
#&SavePrefs();
$afm="on";
$title = "choose";
$OUT .= <<__STOP__;
<tr><td class="$bgcolor">
We recommend that you have a $CASE.struct_supergroup file from the nonmagnetic case.
<p>Is spin flipped on the AF-atom<br>
and zero spin on non-magentic atoms<br>
in $CASE.inst?
<ul>
<FORM ACTION=/exec/initlapw.pl>
<INPUT TYPE=HIDDEN NAME=SID VALUE=$SID>
<INPUT TYPE=HIDDEN NAME=phase VALUE=11>
<INPUT TYPE=HIDDEN NAME=doinst VALUE=1>
<INPUT TYPE=SUBMIT VALUE="edit $CASE.inst">
</FORM>
<FORM ACTION=/exec/initlapw.pl>
<INPUT TYPE=HIDDEN NAME=SID VALUE=$SID>
<INPUT TYPE=HIDDEN NAME=phase VALUE=27>
<INPUT TYPE=SUBMIT VALUE="Yes, continue">
</FORM>
</UL>
__STOP__
$OUT .= "</td></tr>";
&reset;
}


$myphase = 28;
&mark if ($phase == $myphase);
$title = "x afminput";
$precmd  = "x";
$prog    = "afminput";
$next    = "initlapw"; 
&Execute();
&reset;


$myphase = 29;
&mark if ($phase == $myphase);
$title = "view outputafminput";
$file    = "$DIR/$CASE.outputafminput";
&Edit();
&reset;

$myphase = 30;
&mark if ($phase == $myphase);
$title = "view inclmcopy_st";
$file    = "$DIR/$CASE.inclmcopy_st";
&Edit();
&reset;

$myphase=31;
if($phase==$myphase) {
    &mark;
    if ($doit) {
	$dum = qx(cp $DIR/$CASE.inclmcopy_st $DIR/$CASE.inclmcopy);
	$OUT .= "<tr><td class=\"$bgcolor\">inclmcopy created</td></tr>";
	$endtext="<h2>Initialization done</h2><p class=\"info\">Continue with <A HREF=/exec/scf.pl?SID=$SID&spinpol=$spinpol>run SCF</A></p>";
    } else {
	$OUT .= <<__STOP__;
<tr><td class="$bgcolor">
<FORM ACTION=/exec/initlapw.pl>
<INPUT TYPE=HIDDEN NAME=SID VALUE=$SID>
<INPUT TYPE=HIDDEN NAME=phase VALUE=30>
<INPUT TYPE=HIDDEN NAME=doit VALUE=1>
<INPUT TYPE=SUBMIT VALUE="copy inclmcopy_st">
</td></tr>
__STOP__
}
}
}
}

#
#  ------------------------------------------
#  
$OUT .= "</table>";
$OUT .= "</td></tr></table>";
$OUT .= $endtext;

if ($endphase >0) {
    $phase = $endphase;
}
system "echo $phase >$DIR/.initphase";

PrintPage("init_lapw",$OUT);


