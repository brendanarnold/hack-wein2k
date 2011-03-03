#!/usr/bin/perl

require "../../libs/w2web.pl";
require "../../libs/telnes2.pl";

#$debug=1;

&GetInput;
&GetSession;
&ShowParms if $debug;

$OUT .= "<H2>TELNES.2</H2>";
&RequiredFile("innes");
&InnesRead;

$next="continue with TELNES.2";
$nexturl="/exec/telnes2.pl?SID=$SID";
$nextinteractive=1;

#$OUT .= "<h3 class=red>Expert mode</h3>" if $elnesexpert;
if ( $elnesexpert ) {
    $OUT .= "<h3 class=red>Expert mode</h3>"; } else {
	$OUT .= "<h3 class=red>Simple mode (switch to expert mode in CONFIGURATION)</h3> "; }
if ($t_orient && !$t_xqtl) {
    #check for LXDOS=3
    $checkfile = "$WIENROOT/SRC_lapw2/modules.F";
    @umps = split(/=/,qx(grep LXDOS= $checkfile));
    if (@umps[1]==1) {
	$OUT .= <<__STOP__;
<p class="info red">
<b>WARNING:</b><br>
for orientation dependent spectra using "lapw2 -qtl" you MUST recompile
LAPW2 with LXDOS set to 3 in file $checkfile !!!</b>
</p>
__STOP__
    }
}

&RequiredFile("innes");
&InitTask();

$OUT.="<TABLE>";

if ($plot) {
    $OUT .=  <<__STOP__;
<TR><TD class="task" colspan=2>
<A HREF="$nexturl">Show full menu</A>
</TD></TR>
__STOP__

} else {

    $OUT .=  <<__STOP__;
<TR><TD class="task">
<FORM ACTION="/util/innesgen.pl">
<INPUT NAME=SID VALUE=$SID TYPE=HIDDEN>
<INPUT NAME=file VALUE="$DIR/$CASE.innes" TYPE=HIDDEN>
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
<INPUT TYPE=SUBMIT VALUE="edit $CASE.innes">
Edit input-file for ELNES (InnesGen<sup>TM</sup>)
</FORM>
</TD></TR>

<TR><TD class="taskoption">
<b>Only if you want to include states with higher energy</b>
</TD></TR>
<TR><TD class="taskoption">
<FORM ACTION="/util/edit.pl">
<INPUT NAME=SID VALUE=$SID TYPE=HIDDEN>
<INPUT NAME=file VALUE="$DIR/$CASE.in1$filec" TYPE=HIDDEN>
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
<INPUT TYPE=SUBMIT VALUE="edit $CASE.in1$filec">
Edit in1$filec 
</FORM>
</TD></TR>

<TR><TD class="taskoption">
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="lapw1">
$myspin
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
<INPUT TYPE=SUBMIT VALUE="x lapw1 $mycomplex $myspinopt $mypara">
Calculate eigenvalues 
$ni
</FORM>
</TD></TR>
__STOP__

if ($t_xqtl) {
$OUT .=  <<__STOP__;
<TR><TD class="task">
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="xqtl">
$myspin
<INPUT NAME=qtl  VALUE="on" TYPE=HIDDEN>
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
<INPUT TYPE=SUBMIT VALUE="x xqtl $mycomplex $myspinopt $mypara">
Calculate partial charges
$ni
</FORM>
</TD></TR>
__STOP__
} else {
$OUT .=  <<__STOP__;
<TR><TD class="task">
<FORM ACTION=/util/elneslapw2.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="lapw2">
$myspin
<INPUT NAME=qtl  VALUE="on" TYPE=HIDDEN>
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
<INPUT TYPE=SUBMIT VALUE="x lapw2 -qtl $mycomplex $myspinopt $mypara">
Calculate partial charges
$ni
</FORM>
</TD></TR>
__STOP__
}

$OUT .=  <<__STOP__;
<TR><TD class="task">
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="telnes2">
$myspin
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
<INPUT TYPE=SUBMIT VALUE="x telnes2 $mycomplex $myspinopt">
Calculate ELNES spectra
$ni
</FORM>
</TD></TR>

<TR><TD class="task">
<FORM ACTION="/util/edit.pl">
<INPUT NAME=SID VALUE=$SID TYPE=HIDDEN>
<INPUT NAME=file VALUE="$DIR/$CASE.outputelnes" TYPE=HIDDEN>
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
<INPUT TYPE=SUBMIT VALUE="view $CASE.outputelnes">
display $CASE.outputelnes (optional)
</FORM>
</TD></TR>

<TR><TD class="task">
<FORM ACTION="/util/edit.pl">
<INPUT NAME=SID VALUE=$SID TYPE=HIDDEN>
<INPUT NAME=file VALUE="$DIR/$CASE.inb" TYPE=HIDDEN>
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
<INPUT TYPE=SUBMIT VALUE="edit $CASE.inb">
Edit input-file for BROADENING
</FORM>
</TD></TR>

<TR><TD class="task">
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="broadening">
$myspin
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
<INPUT TYPE=SUBMIT VALUE="x broadening $mycomplex $myspinopt">
Broaden the spectrum
$ni
</FORM>
</TD></TR>
__STOP__
}

if($plot){
  $OUT .= <<__STOP__;
<TR><TD class="task">
__STOP__

$myform1 =  <<__STOP__;
<FORM ACTION=/exec/telnes2.pl METHOD=POST>
__STOP__

$myform2 =  <<__STOP__;
<SELECT NAME="stype">
<option value="broad">ELNES (broadened)
<option value="plain">ELNES (unbroadened)
</SELECT>
<INPUT TYPE=hidden NAME=doit VALUE="1">
<INPUT TYPE=hidden NAME=plot VALUE="1">
$myspin
<INPUT TYPE=SUBMIT VALUE="plot">
Plot ELNES
<br> 
Set ranges (optional):
<br>
xmin=<INPUT NAME=xmin VALUE="$xmin" SIZE=5>
xmax=<INPUT NAME=xmax VALUE="$xmax" SIZE=5>
ymin=<INPUT NAME=ymin VALUE="$ymin" SIZE=5>
ymax=<INPUT NAME=ymax VALUE="$ymax" SIZE=5>
</FORM>
</TD></TR>
__STOP__

if($doit) {
    #$OUT .= "Type of spectrum:<br>";
    if ($stype =~ /plain/) {
	$label="ELNES (plain)";
	$sfile="$DIR/$CASE.elnes$updn";
    } elsif ($stype =~ /broad/) {
	$label="ELNES (broadened)";
	$sfile="$DIR/$CASE.broadspec$updn";
    }
    
    $col=2;
    $tmp="$DIR/:elnes";
    $xlabel="eV";
    $axis="yzeroaxis";
    $ylabel="$atom $label";
    $plotfile="$tempdir/$SID-$$";
    unless(open(FILE,">$tmp")) {
      &CGIError("Can't write file $fname.\n");
      exit;
    }
    print FILE <<__STOP__;
show all
set terminal png
set output '$plotfile.png'
set title '$dosfile'
set data style lines
set xlabel "$xlabel"
set ylabel "$ylabel"
set xrange [$xmin:$xmax]
set yrange [$ymin:$ymax]
set $axis
plot '$sfile' using 1:$col title '$label'
set terminal postscript
set output '$plotfile.ps'
replot
__STOP__
close(FILE);

    $umps = qx(cd $DIR;gnuplot $tmp 2>&1);
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
<TR><TD class="task">
<FORM ACTION=/exec/telnes2.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT TYPE=hidden NAME=plot VALUE="1">
$myspin
<INPUT TYPE=SUBMIT VALUE="plot">
Plot ELNES
</FORM>
</TD></TR>
__STOP__
}

$OUT.="</TABLE>";

PrintPage("Context",$OUT);
