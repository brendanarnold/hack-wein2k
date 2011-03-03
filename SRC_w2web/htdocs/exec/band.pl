#!/usr/bin/perl

require "../../libs/w2web.pl";
require "../../libs/struct.pl";

#$debug=1;
$spin="up";
$myspinopt="";

&GetInput;
&GetSession;
&StructRead;
$prefspace="scf";
&GetPrefs();


$bcc="";
$fcc="";
$hcp="";
$simple="";

$bcc=" SELECTED" if ($s_lattice =~ /B/);
$fcc=" SELECTED" if ($s_lattice =~ /F/);
$hcp=" SELECTED" if ($s_lattice =~ /H/);
$simple=" SELECTED" if ($s_lattice =~ /P/);

$next="continue with bandstructure";
$nexturl="/exec/band.pl?SID=$SID";
$nextinteractive=1;

$OUT .= "<h2>Band structure</h2>";

&RequiredFile("insp");
&InitTask();

$OUT.="<TABLE>";


if ($plot) {
$OUT .=  <<__STOP__;
<TR><TD class="task">
<A HREF="$nexturl">Show full menu</A>
</TD></TR>
__STOP__

} else {

if ($ENV{'XCRYSDEN_TOPDIR'}) {
$OUT .=  <<__STOP__;
<TR><TD class="task">
<FORM ACTION=/util/bandxcrys.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT TYPE=SUBMIT VALUE="Generate k-mesh using XCrysden"> (save klist as xcrysden.klist)
</FORM>
</TD></TR>
__STOP__
}

$OUT .=  <<__STOP__;

<TR><TD class="task">
<FORM ACTION="/util/addklist.pl">
<SELECT NAME="klist">
<OPTION VALUE="bcc" $bcc>bcc
<OPTION VALUE="fcc" $fcc>fcc
<OPTION VALUE="hcp" $hcp>hcp
<OPTION VALUE="simple_cubic" $simple>simple_cubic
<OPTION VALUE="xcrysden">from xcrysden
</SELECT>
<INPUT NAME=SID VALUE=$SID TYPE=HIDDEN>
<INPUT TYPE=SUBMIT VALUE="create $CASE.klist_band">
&nbsp;&nbsp;<A HREF="http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-table?from=kv" TARGET=spgrp>Brillouinzones from Bilbao Cryst Server</A><br>
</FORM>
</TD></TR>

</FORM>
</TD></TR>




<TR><TD class="task">
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="lapw1 -band">
<INPUT TYPE=SUBMIT VALUE="x lapw1 -band $mycomplex $myspinopt $mypara">
$myspin
Calculate Eigenvalues
<INPUT NAME=orb TYPE=CHECKBOX $PREFS{'orb'}>&nbsp;orb&nbsp;
$ni
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>

__STOP__
if( $PREFS{'so'} ) {
$OUT .=  <<__STOP__;
<TR><TD class="task">
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="lapwso ">
<INPUT TYPE=SUBMIT VALUE="x lapwso $mycomplex $myspinopt $mypara">
$myspin
Calculate Eigenvalues
<INPUT NAME=orb TYPE=CHECKBOX $PREFS{'orb'}>&nbsp;orb&nbsp;
$ni
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>
__STOP__
}
$OUT .=  <<__STOP__;

<TR><TD class="taskoption">
<b>needed only for continuous lines in the plot (not for non-symmorphic
spacegroups)!</b></td></tr>
<TR><TD class="taskoption">
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="irrep ">
<INPUT TYPE=SUBMIT VALUE="x irrep  $myspinopt $mypara">
$myspin
Calculate irreducible representations
$ni
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>

<TR><TD class="taskoption">
<b>for band character plots only!</b></td></tr>

<TR><TD class="taskoption">
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="lapw2 -band -qtl $mycomplex ">
<INPUT TYPE=SUBMIT VALUE="x lapw2 -band -qtl $mycomplex $myspinopt $mypara">
$myspin
Calculate partial charges ("qtl"-file)
$ni
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>
</td></tr>

<TR><TD class="task">
<FORM ACTION="/util/edit.pl">
<INPUT NAME=SID VALUE=$SID TYPE=HIDDEN>
<INPUT NAME=file VALUE="$DIR/$CASE.insp" TYPE=HIDDEN>
<INPUT NAME=redir VALUE="/exec/band.pl?SID=$SID" TYPE=HIDDEN>
<INPUT TYPE=SUBMIT VALUE="edit $CASE.insp">
Insert correct EF
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>

<TR><TD class="task">
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="spaghetti ">
$myspin
<INPUT TYPE=SUBMIT VALUE="x spaghetti $myspinopt $mypara">
Calculate bandstructure
$ni
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>

__STOP__
}

if($plot){
	$OUT .= <<__STOP__;
<b>We are in plot mode</b><br>
<TR><TD class="task">
__STOP__

$myform1 =  <<__STOP__;
<FORM ACTION=/exec/band.pl METHOD=POST>
__STOP__
$myform2 =  <<__STOP__;
<INPUT TYPE=hidden NAME=doit VALUE="1">
<INPUT TYPE=hidden NAME=plot VALUE="1">
$myspin
<INPUT TYPE=SUBMIT VALUE="plot bandstructure">
</FORM>
</TD></TR>
__STOP__

	if($doit) {

		$plotfile="$tempdir/$SID-$$";

	$umps = qx( cp $DIR/$CASE.spaghetti$spin\_ps $plotfile.ps);
	
	$umps = qx(gs -sDEVICE=jpeg -sOutputFile=$plotfile.jpg -dBATCH -dNOPAUSE $DIR/$CASE.spaghetti$spin\_ps);
	$OUT .= "umps=$umps" if $debug;;
	$OUT .= "<br><IMG SRC=/tmp/$SID-$$.jpg><br clear=all><br>";
	$OUT .= "<A HREF=/tmp/$SID-$$.ps>Download hardcopy in PostScript format</A
>";


		$OUT .= $myform1;
		&PassHiddenParms;
		$OUT .= $myform2
	} else {
		$OUT .= $myform1;
		&PassHiddenParms;
		$OUT .= $myform2

	}
} else {
	$OUT .= <<__STOP__;
<TR><TD class="task">
<FORM ACTION=/exec/band.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT TYPE=hidden NAME=plot VALUE="1">
<INPUT TYPE=hidden NAME=doit VALUE="1">
$myspin
<INPUT TYPE=SUBMIT VALUE="plot bandstructure">
Plot bandstructure
</FORM>
</TD></TR>
__STOP__
}
$OUT .= <<__STOP__;
</TABLE>
__STOP__


PrintPage("Context",$OUT);


