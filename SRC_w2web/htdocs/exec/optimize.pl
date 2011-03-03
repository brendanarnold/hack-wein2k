#!/usr/bin/perl

require "../../libs/w2web.pl";

$SCFFILES="*";
&GetInput;
&GetSession;
&ExeTypes();

$mytest="$DIR/$CASE.broyd1";
if ( -e $mytest && !$doit) {
    $OUT .= <<__STOP__;
<H2>Broyden files exist from previous scf cycle!</H2>
<p>
You should most likely "save" this calculation before running a volume or c/a 
optimization.
</p>
<p>
Do you want to
</p>
<ul>
<li><A HREF="/util/savelapw.pl?SID=$SID">Save Calculation with save_lapw</A>
</ul>
<ul>
<li><A HREF="/util/delete.pl?SID=$SID&file=$DIR/$CASE.broyd*">Remove the
files $CASE.broyd[1|2] </A>
</ul>
<p>or </p>
<ul>
<li> <A HREF="/exec/optimize.pl?SID=$SID&doit=1">run anyway...</A> (very unlikely, that this is a good choice)
</ul>
__STOP__
} else {

    $OUT .= "<h2>Optimize volume or c/a-ratio</h2>";

    $next="continue with optimizer";
    $nexturl="/exec/optimize.pl?SID=$SID";
    $nextinteractive=1;

$OUT .=  <<__STOP__;
<FORM ACTION=/exec/executor.pl METHOD=post>
__STOP__
    &PassHiddenParms;

    $xxout = "${CASE}.struct";
    if (-e "$DIR/${CASE}_initial.struct") { $xxout = "${CASE}_initial.struct"};
    $OUT .=  <<__STOP__;
<TABLE>

<TR><TD class="task">
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="optimize">
<INPUT TYPE=SUBMIT VALUE="x optimize">
Generate structure files from $xxout  
<INPUT TYPE=hidden NAME=nextinteractive VALUE=1>
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl&doit=2">
</FORM>
</TD></TR>

<TR><TD class="task">
<FORM ACTION="/util/edit.pl">
<INPUT NAME=SID VALUE=$SID TYPE=HIDDEN>
<INPUT NAME=file VALUE="$DIR/optimize.job" TYPE=HIDDEN>
<INPUT TYPE=SUBMIT VALUE="edit optimize.job">
    Uncomment "x dstart" or "cp clmsum"; change options in run_lapw, save_lapw,...
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl&doit=2">
</FORM>
</TD></TR>

<TR><TD class="task">
<FORM ACTION=/exec/executor.pl METHOD=post>
__STOP__
&PassHiddenParms;

$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="">
<INPUT NAME=prog TYPE=HIDDEN VALUE="./optimize.job">
<INPUT TYPE=SUBMIT VALUE="run optimize.job">
$exetypes<br>
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl&doit=2">
</FORM>
</TD></TR>

__STOP__

if($plot){
    $OUT .= <<__STOP__;
<TR><TD class="task">
<FORM ACTION=/exec/optimize.pl METHOD=POST>
__STOP__
&PassHiddenParms();
  $OUT .= <<__STOP__;

<SELECT name=type >
<OPTION VALUE="vol">E vs. volume
<OPTION VALUE="coa">E vs. c/a
</SELECT>
<INPUT TYPE=hidden NAME=doit VALUE="2">
<INPUT TYPE=hidden NAME=plot VALUE="2">
<INPUT TYPE=SUBMIT VALUE="plot">
energy curve using
<INPUT NAME="SCFFILES" VALUE="$SCFFILES">
.scf   files
</FORM>
</TD></TR>
__STOP__

    if($doit == 2) {
	@files = sort (glob ("$DIR/${SCFFILES}.scf"));
	if (-e "$DIR/$CASE.vol") { $umps = qx( rm $DIR/$CASE.vol)};
	if ($type =~ /vol/) {
	    unless(open(FILE,">$DIR/$CASE.vol")) {
		&CGIError("Can't write file $fname.\n");
		exit;
	    }
	    foreach $i (@files) {
		$OUT .= "-> $i<br>";
		$x = qx(grep :VOL $i | tail -1 | cut -f2 -d= );
		chomp($x);
		$y = qx(grep :ENE $i | tail -1 | cut -f2 -d= );
		chomp($y);
		print FILE "$x $y\n";
	    }
	    close(FILE);

	    $umps = qx( cd $DIR && echo " , " | x eosfit );
	    $tmp = qx(grep V0 $DIR/$CASE.outputeos | tail -4 |head -1) ;
	    @murna = split(" ",$tmp);
	    
	    $plotfile="$tempdir/$SID-$$";
	    unless(open(FILE,">$DIR/:eplot")) {
		&CGIError("Can't write file $fname.\n");
		exit;
	    }

	    $t1=$murna[0];
	    shift(@murna);
	    print FILE <<__STOP__;
set terminal png
set output '$plotfile.png'
set format y "%.4f"
set title "$CASE"
set xlabel "Volume [a.u.^3]"
set ylabel "Energy [Ry]"
plot "$CASE.vol" title "Murnaghan: $t1" w p, "$CASE.eosfit" title "@murna" w l 
set terminal postscript
set output '$plotfile.ps'
replot
__STOP__

close(FILE);
	    
	    $umps = qx(cd $DIR;gnuplot ":eplot" 2>&1);
	    $OUT .= "<br><IMG SRC=/tmp/$SID-$$.png><br clear=all><br>";
	    $OUT .= "<A HREF=/tmp/$SID-$$.ps>Download hardcopy in PostScript format</A>";
	    $OUT .= "<br><pre>";
	    $OUT .= qx(cat $DIR/$CASE.outputeos);
	    
	} else {
	    $umps = qx(cd $DIR;grepline_lapw :ene '@files' 1 > $CASE.analysis);
	    $umps = qx(cd $DIR;eplot -t coa -p);
	    $umps = qx(cp $DIR/$CASE.c_over_a.png $tempdir/$SID-$$.png);
	    $umps = qx(cp $DIR/eplot.ps $tempdir/$SID-$$.ps);
	    $OUT .= "<br><pre>";
	    $OUT .= qx(cd $DIR;grepline_lapw :ene '@files' 1);
	    $OUT .= qx(cd $DIR;echo '  ';echo "Fit of:  E = a1 + a2*x + a3*x^2 + a4*x^3 + a5*x^4";tail -5 fit.log);
	    $OUT .= "</pre>";
	    $OUT .= "<IMG SRC=/tmp/$SID-$$.png><br clear=all><br>";
	    $OUT .= "<A HREF=/tmp/$SID-$$.ps>Download hardcopy in PostScript format</A>";
	}
	
    }

} else {
    $OUT .= <<__STOP__;
<TR><TD class="task">
<FORM ACTION=/exec/optimize.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT TYPE=hidden NAME=doit VALUE="1">
<INPUT TYPE=hidden NAME=plot VALUE="1">
<INPUT TYPE=SUBMIT VALUE="plot">
Plot energy curve
</FORM>
</TD></TR>
__STOP__
}

$OUT.="</table>";


}

PrintPage("Context",$OUT);


