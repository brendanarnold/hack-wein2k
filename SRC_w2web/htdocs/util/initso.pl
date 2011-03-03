#!/usr/bin/perl

require "../../libs/w2web.pl";
#$debug=1;

&GetInput;
&GetSession;

$OUT .= "<h2>Initialization of spin-orbit calculations</h2> ";


$next="continue with spin-orbit initialization";
$nexturl="/util/initso.pl?SID=$SID";
$nextinteractive=1;

#&InitTask();
	$filec="";
	if($complex =~ /CHECKED/ ) {
		$mycomplex="-c";
		$filec="c";
	} else {
        $dum = qx(cp $DIR/$CASE.in2 $DIR/$CASE.in2c);
        $OUT .= "<p> $CASE.in2c has been created"; 
    }

&RequiredFile("inso");


$OUT.="<TABLE BGCOLOR=$green>";

$OUT.=<<__STOP__;
</TD></TR>
<TR><TD BGCOLOR=$gray>
<FORM ACTION="/util/edit.pl">
<INPUT NAME=SID VALUE=$SID TYPE=HIDDEN>
<INPUT NAME=file VALUE="$DIR/$CASE.inso" TYPE=HIDDEN>
<INPUT NAME=redir VALUE="/exec/initso.pl?SID=$SID" TYPE=HIDDEN>
<INPUT TYPE=SUBMIT VALUE="edit $CASE.inso">
Select magnetization direction, RLOs, SO on/off
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>

<TR><TD BGCOLOR=$gray>
<FORM ACTION="/util/edit.pl">
<INPUT NAME=SID VALUE=$SID TYPE=HIDDEN>
<INPUT NAME=file VALUE="$DIR/$CASE.in1$filec" TYPE=HIDDEN>
<INPUT NAME=redir VALUE="/exec/initso.pl?SID=$SID" TYPE=HIDDEN>
<INPUT TYPE=SUBMIT VALUE="edit $CASE.in1$filec">
set larger EMAX in energy window
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>

__STOP__

	if($spinpol=~ /CHECKED/ ) {
#if($myspin) {
$OUT.=<<__STOP__;
<TR><TD>
This is a spin-polarized system. SO may reduce symmetry.
</TD></TR>
<TR><TD BGCOLOR=$gray>
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="symmetso">
<INPUT TYPE=SUBMIT VALUE="x symmetso $mycomplex">
<INPUT TYPE=hidden NAME=nextinteractive VALUE=1>
Determines symmetry in spinpolarized case
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>

<TR><TD BGCOLOR=$gray>
<FORM ACTION="/util/edit.pl">
<INPUT NAME=SID VALUE=$SID TYPE=HIDDEN>
<INPUT NAME=file VALUE="$DIR/$CASE.outsymso" TYPE=HIDDEN>
<INPUT NAME=redir VALUE="/exec/initso.pl?SID=$SID" TYPE=HIDDEN>
<INPUT TYPE=SUBMIT VALUE="edit $CASE.outsymso">
view $CASE.outsymso
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>

<TR><TD>
A new setup for SO calculations has been created (_so). If you commit the next step will create new $CASE.struct, in1$filec, in2c, inc, clmsum/up/dn files. PLEASE "save_lapw" any previous calculation.

</TD></TR>
						  
__STOP__

if ($copyit) {  
	$dum = qx(cp $DIR/$CASE.struct_so $DIR/$CASE.struct) if(-s "$DIR/$CASE.struct_so" ) ;
	$dum = qx(cp $DIR/$CASE.inc_so $DIR/$CASE.inc) if(-s "$DIR/$CASE.inc_so" ) ;
	$dum = qx(cp $DIR/$CASE.in2c_so $DIR/$CASE.in2c) if(-s "$DIR/$CASE.in2c_so" ) ;
	$dum = qx(cp $DIR/$CASE.in2_so $DIR/$CASE.in2c) if(! -s "$DIR/$CASE.in2c_so" && -s "$DIR/$CASE.in2_so" ) ;
	$dum = qx(cp $DIR/$CASE.in1c_so $DIR/$CASE.in1c) if(-s "$DIR/$CASE.in1c_so" ) ;
	$dum = qx(cp $DIR/$CASE.in1_so $DIR/$CASE.in1) if(-s "$DIR/$CASE.in1_so" ) ;
	$dum = qx(cp $DIR/$CASE.clmsum_so $DIR/$CASE.clmsum) if(-s "$DIR/$CASE.clmsum_so" ) ;
	$dum = qx(cp $DIR/$CASE.clmup_so $DIR/$CASE.clmup) if(-s "$DIR/$CASE.clmup_so" ) ;
	$dum = qx(cp $DIR/$CASE.clmdn_so $DIR/$CASE.clmdn) if(-s "$DIR/$CASE.clmdn_so" ) ;
	$dum = qx(cp $DIR/$CASE.vspup_so $DIR/$CASE.vspup) if(-s "$DIR/$CASE.vspup_so" ) ;
	$dum = qx(cp $DIR/$CASE.vspdn_so $DIR/$CASE.vspdn) if(-s "$DIR/$CASE.vspdn_so" ) ;
	$dum = qx(cp $DIR/$CASE.vnsup_so $DIR/$CASE.vnsup) if(-s "$DIR/$CASE.vnsup_so" ) ;
	$dum = qx(cp $DIR/$CASE.vnsdn_so $DIR/$CASE.vnsdn) if(-s "$DIR/$CASE.vnsdn_so" ) ;
	$OUT .= <<__STOP__;
<TR><TD>
new files generated (adapt $CASE.indmc or inorb manually)
</TD></TR>
__STOP__

} else {
	$OUT .= <<__STOP__;
<TR><TD BGCOLOR=$gray>
<FORM ACTION=/util/initso.pl>
<INPUT TYPE=HIDDEN NAME=SID VALUE=$SID>
<INPUT TYPE=HIDDEN NAME=copyit VALUE=1>
<INPUT TYPE=SUBMIT VALUE="Prepare new input files">
</FORM>
</TD></TR>
__STOP__
}
#&reset;
	if (! -z "$DIR/$CASE.ksym" && -e "$DIR/$CASE.ksym" ) {
	$OUT .= <<__STOP__;
<TR><TD>
$CASE.ksym detected, which contains proper symmetry operations for KGEN. Note: When changing the k-mesh later, always run "x kgen -so".
</TD></TR>

<TR><TD BGCOLOR=$gray>
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="kgen">
<INPUT TYPE=SUBMIT VALUE="x kgen -so">
<INPUT TYPE=hidden NAME=so VALUE="on">
<INPUT TYPE=hidden NAME=nextinteractive VALUE=1>
Generate k-mesh with proper SO-symmetry
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>

<TR><TD BGCOLOR=$gray>
<FORM ACTION="/util/edit.pl">
<INPUT NAME=SID VALUE=$SID TYPE=HIDDEN>
<INPUT NAME=file VALUE="$DIR/$CASE.klist" TYPE=HIDDEN>
<INPUT NAME=redir VALUE="/exec/initso.pl?SID=$SID" TYPE=HIDDEN>
<INPUT TYPE=SUBMIT VALUE="edit $CASE.klist">
view $CASE.klist
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>

__STOP__

	} 

} else {
	$OUT .= <<__STOP__;
<tr><td> System not $spinpol spinpolarized 
</TD></TR>
__STOP__
}
$OUT.="</TABLE>";

PrintPage("Context",$OUT);

