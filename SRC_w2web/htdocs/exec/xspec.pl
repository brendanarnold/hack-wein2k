#!/usr/bin/perl

require "../../libs/w2web.pl";
$spin="up";
&GetInput;
&GetSession;
#&ShowParms;

$next="continue with xspec";
$nexturl="/exec/xspec.pl?SID=$SID";
$nextinteractive=1;

$OUT .= "<H2>XSPEC</H2>";

&RequiredFile(inxs);
&InitTask();


$OUT.="<TABLE>";

if ($plot) {
$OUT .=  <<__STOP__;
<TR><TD class="task">
<A HREF="$nexturl">Show full menu</A>
</TD></TR>
__STOP__
} else {
	$OUT .=  <<__STOP__;

<TR><TD class="taskoption">
<b>If you want to include states with higher energy</b>
</td></tr>
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


<TR><TD class="task">
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="lapw2">
$myspin
<INPUT NAME=qtl  VALUE="on" TYPE=HIDDEN>
<INPUT TYPE=SUBMIT VALUE="x lapw2 -qtl $mycomplex $myspinopt $mypara">
Calculate partial charges
$ni
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>


<TR><TD class="task">
<FORM ACTION="/util/edit.pl">
<INPUT NAME=SID VALUE=$SID TYPE=HIDDEN>
<INPUT NAME=file VALUE="$DIR/$CASE.inxs" TYPE=HIDDEN>
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
<INPUT TYPE=SUBMIT VALUE="edit $CASE.inxs">
Edit input-file for XSPEC
</FORM>
</TD></TR>

<TR><TD class="task">
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="xspec">
$myspin
<INPUT TYPE=SUBMIT VALUE="x xspec $myspinopt">
Calculate X-ray spectra
$ni
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>

__STOP__
}

if($plot){
  $OUT .= <<__STOP__;
<TR><TD class="task">
<b>We are in plot mode</b><br>
__STOP__

$myform1 =  <<__STOP__;
<FORM ACTION=/exec/xspec.pl METHOD=POST>
__STOP__

$atom = qx(head -2 <$file.inxs|tail -1|awk '{print $1}');
$n = qx(head -3 <$file.inxs|tail -1|awk '{print $1}');
$l = qx(head -4 <$file.inxs|tail -1|awk '{print $1}');

if ($l =~ /0/) {
	$label = "K";
} elsif ($i =~ /1/) {
	$label = "L3";
} elsif ($i =~ /2/) {
	$label = "M5";
} else {
	$label ="";
}

$atomname = qx(grep Z: $file.struct|head -$atom|tail -1|cut -f1 -d\" \");
$label = "$atomname $label";

		if ($spectype =~ /^xspec/) {
                        $xsel=" SELECTED";
		        $column="2";
                        $txsel="";
                        $m1sel="";
                        $m2sel="";
                        $coresel="";
		} elsif ($spectype =~ /^txspec/) {
                        $xsel="";
                        $txsel=" SELECTED";
                        $m1sel="";
                        $m2sel="";
                        $coresel="";
		} elsif ($spectype =~ /^m1/) {
		        $column="2";
                        $xsel="";
                        $txsel="";
                        $m1sel=" SELECTED";
                        $m2sel="";
                        $coresel="";
		} elsif ($spectype =~ /^m2/) {
		        $column="2";
                        $xsel="";
                        $txsel="";
                        $m1sel="";
                        $m2sel=" SELECTED";
                        $coresel="";
		} elsif ($spectype =~ /^corewfx/) {
		        $column="2";
                        $xsel="";
                        $txsel="";
                        $m1sel="";
                        $m2sel="";
                        $coresel=" SELECTED";
		}

$myform2 =  <<__STOP__;
<SELECT NAME="spectype">
<option value="xspec" $xsel>broadened spectrum
<option value="txspec" $txsel>unbroadened spectrum
<option value="m1" $m1sel>matrix elements L+1
<option value="m2" $m2sel>matrix elements L-1
<option value="corewfx" $coresel>core wavefunction
</SELECT>
column=<INPUT NAME=column VALUE="$column" SIZE=2>
<INPUT TYPE=hidden NAME=doit VALUE="1">
<INPUT TYPE=hidden NAME=plot VALUE="1">
$myspin
<INPUT TYPE=SUBMIT VALUE="plot">
 (partial DOS must be plotted in DOS-Task)
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
    $xlabel="eV";
    $OUT .= "Spectrum type $spectype selected:<br>";
		if ($spectype =~ /^xspec/) {
		        $column="2";
			$label.=" (broadened)";
		} elsif ($spectype =~ /^txspec/) {
			$label.=" (unbroadened)";
			$label.=" total" if ($column =~ "2");
			$label.=" L+1" if ($column =~ "3");
			$label.=" L-1" if ($column =~ "4");
		} elsif ($spectype =~ /^m1/) {
		        $column="2";
			$label.=" (matrix elemnts L+1)";
		} elsif ($spectype =~ /^m2/) {
		        $column="2";
			$label.=" (matrix elements L-1)";
		} elsif ($spectype =~ /^corewfx/) {
		        $xlabel="Radius (bohr)";
		        $column="2";
			$label.=" (core wavefunction)";
		}
    $tmp1="$DIR/:xspec1";
    $axis="yzeroaxis";
    $ylabel="$atom $label";
    $dosfile="$DIR/$CASE.$spectype$updn";
    $plotfile="$tempdir/$SID-$$";
    unless(open(FILE,">$tmp1")) {
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
plot '$DIR/$CASE.$spectype$updn' using 1:$column title '$label'
set terminal postscript
set output '$plotfile.ps'
replot
__STOP__
close(FILE);

  $umps = qx(cd $DIR;gnuplot $tmp1 2>&1);
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
<FORM ACTION=/exec/xspec.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT TYPE=hidden NAME=plot VALUE="1">
$myspin
<INPUT TYPE=SUBMIT VALUE="plot">
Plot XSPEC
</FORM>
</TD></TR>
__STOP__
}

$OUT.="</TABLE>";



PrintPage("Context",$OUT);


