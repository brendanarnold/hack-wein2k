#!/usr/bin/perl

require "../../libs/w2web.pl";
#$debug=1;
$spin="up";
$myspinopt="";

&GetInput;
&GetSession;

$OUT .= "<h2>Optical properties</h2>";

$next="continue with optics";
$nexturl="/exec/optic.pl?SID=$SID";
$nextinteractive=1;

&InitTask();
&RequiredFile("inop");
&RequiredFile("injoint");
&RequiredFile("inkram");


$OUT.="<TABLE>";

if ($plot) {
	$OUT .=  <<__STOP__;
<TR><TD class="task">
<A HREF="$nexturl">Show full menu</A>
</TD></TR>
__STOP__

} else {

	$OUT .=  <<__STOP__;
<TR><TD class="task">
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="lapw2 -fermi $mycomplex">
<INPUT TYPE=SUBMIT VALUE="x lapw2 -fermi $mycomplex $myspinopt $mypara">
$myspin
Calculate weights
$ni
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>

<TR><TD class="task">
<FORM ACTION="/util/edit.pl">
<INPUT NAME=SID VALUE=$SID TYPE=HIDDEN>
<INPUT NAME=file VALUE="$DIR/$CASE.inop" TYPE=HIDDEN>
<INPUT TYPE=SUBMIT VALUE="edit $CASE.inop">
Edit $CASE.inop
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
<INPUT NAME=prog TYPE=HIDDEN VALUE="optic">
<INPUT TYPE=SUBMIT VALUE="x optic $mycomplex $myspinopt $mypara">
$myspin
Calculate optic
$ni
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>

<TR><TD class="task">
<FORM ACTION="/util/edit.pl">
<INPUT NAME=SID VALUE=$SID TYPE=HIDDEN>
<INPUT NAME=file VALUE="$DIR/$CASE.injoint" TYPE=HIDDEN>
<INPUT TYPE=SUBMIT VALUE="edit $CASE.injoint">
Edit $CASE.injoint
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
<INPUT NAME=prog TYPE=HIDDEN VALUE="joint">
$myspin
<INPUT TYPE=SUBMIT VALUE="x joint  $myspinopt">
Calculate optic
$ni
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>


<TR><TD class="task">
<FORM ACTION="/util/edit.pl">
<INPUT NAME=SID VALUE=$SID TYPE=HIDDEN>
<INPUT NAME=file VALUE="$DIR/$CASE.inkram" TYPE=HIDDEN>
<INPUT TYPE=SUBMIT VALUE="edit $CASE.inkram">
edit $CASE.inkram
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
<INPUT NAME=prog TYPE=HIDDEN VALUE="kram">
$myspin
<INPUT TYPE=SUBMIT VALUE="x kram  $myspinopt">
$ni
Calculate kram
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
<FORM ACTION=/exec/optic.pl METHOD=POST>
__STOP__
$myform2 =  <<__STOP__;
<SELECT NAME="selfile">
<option value="joint">$CASE.joint$updn
<option value="epsilon">$CASE.epsilon$updn
<option value="sigmak">$CASE.sigmak$updn
<option value="sigma_intra">$CASE.sigma_intra$updn
<option value="eloss">$CASE.eloss$updn
<option value="sumrules">$CASE.sumrules$updn
</SELECT>
<SELECT NAME="col">
<option value=2>2
<option value=3>3
<option value=4>4
<option value=5>5
<option value=6>6
<option value=7>7
</SELECT>
<INPUT TYPE=hidden NAME=doit VALUE="1">
<INPUT TYPE=hidden NAME=plot VALUE="1">
<INPUT TYPE=SUBMIT VALUE="plot">
Plot optical properties
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

		$plotfile="$tempdir/$SID-$$";
		$infile="$DIR/$CASE.";
		$titline=1;
		$delline=0;

    if ($selfile =~ /joint/) {
			$infile.="joint$updn";
			$titline=2;
			$delline=3;
		} elsif ($selfile =~ /epsilon/ ) {
			$infile.="epsilon$updn";
			$titline=5;
			$delline=6;
			$plotoptione="set xzeroaxis";
		} elsif ($selfile =~ /sigmak/ ) {
			$infile.="sigmak$updn";
			$titline=7;
			$delline=8;
		} elsif ($selfile =~ /sigma_intra/ ) {
			$infile.="sigma_intra$updn";
			$titline=3;
			$delline=4;
		} elsif ($selfile =~ /eloss/ ) {
			$infile .= "eloss$updn";
			$titline = 4;
			$delline = 6;
		} elsif ($selfile =~ /sumrules/ ) {
			$infile .= "sumrules$updn";
			$titline = 1;
			$delline = 0;
		}


	
		$test = qx(wc -l $infile);
		if ($test == 0) {
			$OUT .= "$infile is empty - no plot produced";
		} else {

			$tmp1 = "$DIR/:opt1";
			$tmp2 = "$DIR/:opt2";

			$units  = "ev";
			$xlabel = "Energy [eV]";
			$axis   = "yzeroaxis";
	
			$umps = qx(sed "1,${delline}d" $infile >$tmp1);
			$tmp  = qx(cat $infile | head -$titline | tail -1|cut -c1-8,15-99);
			@ylabels = split(" ",$tmp);
			$ylabel = $ylabels[$col-1];
			$title = "$infile column $col";

    unless(open(FILE,">$tmp2")) {
      &CGIError("Can't write file $fname.\n");
      exit;
    }
    print FILE <<__STOP__;
show all
set terminal png
set output '$plotfile.png'
set title '$title'
set data style lines
set xrange [$xmin:$xmax]
set yrange [$ymin:$ymax]
set xlabel "$xlabel"
set ylabel "$ylabel"

set $axis
$plotoptions
plot "$tmp1" using 1:$col  title "$fname"
set terminal postscript
set output '$plotfile.ps'
replot
__STOP__
			close(FILE);
			$umps = qx(cd $DIR;gnuplot $tmp2 2>&1);
			$OUT .= "<br><IMG SRC=/tmp/$SID-$$.png><br clear=all><br>";
			$OUT .= "<A HREF=/tmp/$SID-$$.ps>Download hardcopy in PostScript format</A
>";

	
		}
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
<FORM ACTION=/exec/optic.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT TYPE=hidden NAME=plot VALUE="1">
$myspin
<INPUT TYPE=SUBMIT VALUE="opticplot">
Plot optical properties
</FORM>
</TD></TR>
__STOP__
}

$OUT.=<<__STOP__;

</TABLE>
__STOP__


PrintPage("Context",$OUT);


