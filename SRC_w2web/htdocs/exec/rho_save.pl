#!/usr/bin/perl

require "../../libs/w2web.pl";
#$debug=1;
$min="";
$max="";

&GetInput;
&GetSession;

$OUT .= "<h2>Electron density plots </h2> You must have a valid
$CASE.vector file (from an scf calculation). If you don't have it, you must run \" x lapw1 \"  with an appropriate input.";


$next="continue with electron density";
$nexturl="/exec/rho.pl?SID=$SID";
$nextinteractive=1;

&InitTask();

&RequiredFile("in5$filec");


$OUT.="<TABLE BGCOLOR=$green>";

if($plot) {
$OUT.=<<__STOP__;
<TR><TD BGCOLOR=$gray>
<A HREF="$nexturl">Show full menu</A>
</TD></TR>
__STOP__

} else {

$OUT.=<<__STOP__;
</TD></TR>
<TR><TD BGCOLOR=$gray>
<FORM ACTION="/util/edit.pl">
<INPUT NAME=SID VALUE=$SID TYPE=HIDDEN>
<INPUT NAME=file VALUE="$DIR/$CASE.in2$filec" TYPE=HIDDEN>
<INPUT NAME=redir VALUE="/exec/rho.pl?SID=$SID" TYPE=HIDDEN>
<INPUT TYPE=SUBMIT VALUE="edit $CASE.in2$filec">
change EMIN to truncate semicore
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>


<TR><TD BGCOLOR=$gray>
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="lapw2">
<INPUT TYPE=SUBMIT VALUE="x lapw2 $mycomplex $myspinopt">
<INPUT TYPE=hidden NAME=nextinteractive VALUE=1>
Calculate clmval
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>
__STOP__

if ($ENV{'XCRYSDEN_TOPDIR'}) {
$OUT .=  <<__STOP__;
<TR><TD BGCOLOR=$gray>
<FORM ACTION=/util/rhoxcrys.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT TYPE=HIDDEN NAME="wos_tua_ma" VALUE="calc">
<INPUT TYPE=SUBMIT VALUE="Calculate density with XCrysden">
</FORM>
</TD></TR>
__STOP__
}

$OUT .=  <<__STOP__;
<TR><TD BGCOLOR=$gray>
<FORM ACTION="/util/edit.pl">
<INPUT NAME=SID VALUE=$SID TYPE=HIDDEN>
<INPUT NAME=file VALUE="$DIR/$CASE.in5$filec" TYPE=HIDDEN>
<INPUT NAME=redir VALUE="/exec/rho.pl?SID=$SID" TYPE=HIDDEN>
<INPUT TYPE=SUBMIT VALUE="edit $CASE.in5$filec">
Edit input-file 
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>


<TR><TD BGCOLOR=$gray>
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="lapw5">
<INPUT TYPE=SUBMIT VALUE="x lapw5 $mycomplex $myspinopt">
<INPUT TYPE=hidden NAME=nextinteractive VALUE=1>
Calculate density
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>

__STOP__

if ($ENV{'XCRYSDEN_TOPDIR'}) {
$OUT .=  <<__STOP__;
<TR><TD BGCOLOR=$gray>
<FORM ACTION=/util/rhoxcrys.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT TYPE=HIDDEN NAME="wos_tua_ma" VALUE="preview">
<INPUT TYPE=SUBMIT VALUE="Preview density with XCrysden">
</FORM>
</TD></TR>
__STOP__
}

}
if($plot){
	$OUT .= <<__STOP__;
<TR><TD BGCOLOR=$gray>
<b>We are in rhoplot mode</b><br>
__STOP__

$myform1 =  <<__STOP__;
<FORM ACTION=/exec/rho.pl METHOD=POST>
__STOP__
$myform2 =  <<__STOP__;
<INPUT TYPE=hidden NAME=doit VALUE="1">
<INPUT TYPE=hidden NAME=plot VALUE="1">
Min <INPUT NAME=min value=$min>
Max <INPUT NAME=max value=$max><br>
<br>
<INPUT TYPE=SUBMIT VALUE="plot electron density">
</FORM>
</TD></TR>
__STOP__

	if($doit) {

		$tmp1="$DIR/:rho1";
		$tmp2="$DIR/:rho2";
		$file="$DIR/$CASE.rho";


		$dum = qx( reformat <$file >$tmp1);

		$plotfile="$tempdir/$SID-$$";
		unless(open(FILE,">$tmp2")) {
			&CGIError("Can't write file $fname.\n");
			exit;
		}
		print FILE <<__STOP__;
set title '$title'
set data style lines
set contour
set zrange[$min:$max]
set nokey
set hidden3d
set terminal png
set output '$plotfile.png'
splot '$tmp1'
set terminal postscript
set output '$plotfile.ps'
replot
__STOP__
close(FILE);

	$umps = qx(cd $DIR;gnuplot $tmp2 2>&1);
	$OUT .= "umps=$umps" if $debug;;
	$OUT .= "<br><IMG SRC=/tmp/$SID-$$.png><br clear=all><br>";
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
<TR><TD BGCOLOR=$gray>
<FORM ACTION=/exec/rho.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT TYPE=hidden NAME=plot VALUE="1">
$myspin
<INPUT TYPE=SUBMIT VALUE="rhoplot">
Plot Density
</FORM>
</TD></TR>
__STOP__
}

$OUT.=<<__STOP__;

<TR><TD BGCOLOR=$gray>
<FORM ACTION="/util/edit.pl">
<INPUT NAME=SID VALUE=$SID TYPE=HIDDEN>
<INPUT NAME=file VALUE="$DIR/$CASE.in2$filec" TYPE=HIDDEN>
<INPUT NAME=redir VALUE="/exec/rho.pl?SID=$SID" TYPE=HIDDEN>
<INPUT TYPE=SUBMIT VALUE="edit $CASE.in2$filec">
reset EMIN
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>
</TABLE>
</td></tr>
</TABLE>

__STOP__



PrintPage("Context",$OUT);

