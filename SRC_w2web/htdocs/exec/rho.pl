#!/usr/bin/perl

require "../../libs/w2web.pl";
#$debug=1;
$min="";
$max="";
#$lab="";
$CHECKED_la="CHECKED";
$deltaz="0.1";
$CHECKED_3d="CHECKED";
$CHECKED_co="";

&GetInput;
&GetSession;
$prefspace="single";
&GetPrefs;

$OUT .= "<h2>Electron density plots </h2>";
$OUT.= <<__STOP__;
<p class="info">
You must have a valid $CASE.vector file (from an scf calculation). <br>
If you don't have it, you must run "x lapw1"  with an appropriate input.
</p>
__STOP__


$next="continue with electron density";
$nexturl="/exec/rho.pl?SID=$SID";
$nextinteractive=1;

&InitTask();

&RequiredFile("in5$filec");


$OUT.="<TABLE>";

if($plot) {
$OUT.=<<__STOP__;
<TR><TD class="task">
<A HREF="$nexturl">Show full menu</A>
</TD></TR>
__STOP__

} else {

$OUT.=<<__STOP__;
</TD></TR>
<TR><TD class="task">
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


<TR><TD class="task">
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="lapw2">
$myspin
<INPUT TYPE=SUBMIT VALUE="x lapw2 $mycomplex $myspinopt $mypara">
<INPUT TYPE=hidden NAME=nextinteractive VALUE=1>
Calculate clmval
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>

<TR><TD class="taskoption">
<TABLE>
<tr><td class="taskoption" colspan=2 align=center>
<b>For difference densities only !</b></td></tr>
<TR>
<TH class="taskoption">default valence states:</TH>
<TH class="taskoption">non-default valence states:</TH>
</TR><TR>
<TD class="taskoption">&nbsp;</td><TD class="taskoption">
<FORM ACTION="/util/edit.pl">
<INPUT NAME=SID VALUE=$SID TYPE=HIDDEN>
<INPUT NAME=file VALUE="$DIR/$CASE.inst" TYPE=HIDDEN>
<INPUT NAME=redir VALUE="/exec/rho.pl?SID=$SID" TYPE=HIDDEN>
<INPUT TYPE=SUBMIT VALUE="edit $CASE.inst">
put P for all your states 
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</td></tr>
<TR class="taskoption"><TD class="taskoption">
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$PREFS{'sigma'}= 'on';
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="lstart">
<INPUT NAME=sigma TYPE=HIDDEN VALUE=$PREFS{'sigma'}>
<INPUT TYPE=SUBMIT VALUE="x lstart -sigma">
Calculate atomic valence densities
<INPUT TYPE=hidden NAME=nextinteractive VALUE=1> 
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</td><TD class="taskoption">
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="lstart">
<INPUT TYPE=SUBMIT VALUE="x lstart">
Calculate atomic valence densities as defined above
<INPUT TYPE=hidden NAME=nextinteractive VALUE=1> 
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>
</TABLE>

</TD></TR>
__STOP__

if ($ENV{'XCRYSDEN_TOPDIR'}) {
$OUT .=  <<__STOP__;
<TR><TD class="task">
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
<TR><TD class="task">
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


<TR><TD class="task">
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="lapw5">
$myspin
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
<TR><TD class="task">
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
<TR><TD class="task">
<b>We are in rhoplot mode</b><br>
__STOP__

	        if($plttyp =~ /3/ ) {$CHECKED_3d="CHECKED";$CHECKED_co="";}
	        if($plttyp =~ /c/ ) {$CHECKED_co="CHECKED";$CHECKED_3d="";}
	        if($lab eq "on" ) {$CHECKED_la="CHECKED";}
	        if($lab eq "" ) {$CHECKED_la="";}

$myform1 =  <<__STOP__;
<FORM ACTION=/exec/rho.pl METHOD=POST>
__STOP__
$myform2 =  <<__STOP__;
<INPUT TYPE=hidden NAME=doit VALUE="1">
<INPUT TYPE=hidden NAME=plot VALUE="1">
<b>Select plot type:&nbsp; </b> 
3D-plot <INPUT NAME=plttyp TYPE=radio VALUE="3" $CHECKED_3d >&nbsp;&nbsp;&nbsp;
Contur-plot <INPUT NAME=plttyp TYPE=radio VALUE="c" $CHECKED_co >&nbsp;&nbsp;&nbsp;
with labels <INPUT NAME=lab TYPE=CHECKBOX $CHECKED_la ><br>
Min <INPUT NAME=min value="$min" SIZE=6>
Max <INPUT NAME=max value="$max" SIZE=6>
Delta <INPUT NAME=deltaz value="$deltaz" SIZE=6><br>
<INPUT TYPE=SUBMIT VALUE="plot electron density">
</FORM>
</TD></TR>
__STOP__

	if($doit) {

#		$tmp1="$DIR/:rho1";
#		$tmp2="$DIR/:rho2";
		$file="$DIR/$CASE.rho";
		$lab1="n";
		if ($lab eq "on"){$lab1="y"} ;
#		$dum = qx( reformat <$file >$tmp1);

		$plotfile="$tempdir/$SID-$$";
#		unless(open(FILE,">$tmp2")) {
#			&CGIError("Can't write file $fname.\n");
#			exit;
#		}
#		print FILE <<__STOP__;
#set title '$title'
#set data style lines
#set contour
#set zrange[$min:$max]
#set nokey
#set hidden3d
#set terminal png
#set output '$plotfile.png'
#splot '$tmp1'
#set terminal postscript
#set output '$plotfile.ps'
#replot
#__STOP__
#close(FILE);

#	$umps = qx(cd $DIR;gnuplot $tmp2 2>&1);
	
		$umps = qx(cd $DIR;rhoplot $file "$min" "$max" "$plttyp" "$deltaz" $lab1 $plotfile 2>&1);
	$OUT .= "umps=$umps cd $DIR;rhoplot $file $min $max $plttyp $deltaz $plotfile " if $debug;;
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

<TR><TD class="task">
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

