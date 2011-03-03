#!/usr/bin/perl

require "../../libs/w2web.pl";

&GetInput;
&GetSession;
&ExeTypes();

$mytest="$DIR/$CASE.broyd1";
if ( -e $mytest && !$doit) {
  $OUT .= <<__STOP__;
<H2>Broyden files exist from previous scf cycle!</H2>
<p>
This is ok when this scf-run is the first step in the geometry optimization.
However, when any input / struct / kmesh parameters have changed, you should 
"save" this calculation or remove the "broyden" files before running the 
geometry optimization.
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
<li> <A HREF="/exec/min.pl?SID=$SID&doit=1">run anyway...</A> 
</ul>
__STOP__
} else {

$OUT .= "<h2>Geometry minimizer</h2>";

$editit=1;
$file="$DIR/$CASE.inM";
if (-e $file && -s $file ) {} else {
$OUT .= "x pairhess executed and .minrestart / $CASE.inM files generated !<br>\n ";
                $umps=qx(cd $DIR;$WIENROOT/x pairhess);
                $umps=qx(cd $DIR;cp $CASE.inM_st $CASE.inM; cp .minpair .minhess; cp .minpair .minrestart);
#redirectURL("/util/edit.pl?SID=$SID&f=1&file=$DIR/$CASE.inM");
#exit;
    $nexturl = "/exec/min.pl?SID=$SID";
    $next = "min";
&PassHiddenParms;
 $OUT .= <<__STOP__;
<A HREF="/util/edit.pl?SID=$SID&file=$DIR/$CASE.inM&next=$next&nexturl=$nexturl">Optionally edit 
file $CASE.inM </A>
__STOP__
#&RequiredFile("inM");
}


if ($doit == 2) {

	$init= "-NI" if ($initfiles=~/on/);

	$max="";
	$max="-i $maxstruct" if(length($maxstruct)>0);

	$nohes="";
	$nohes="-nohess" if ($nohess=~/on/);

	$save="";
	$save="-s $savestep" if(length($savestep)>0);

	$pp="";
	$pp="-p " if($FORM{'p1'});

	$job="run_lapw -I -i 40 -fc 1.0 $pp";
	$job= "runsp_lapw -I -i 40 -fc 1.0 $pp " if ($spinpol =~ /CHECKED/ );
	$job= $jobfile if (length($jobfile)>0 );

	$cmdline="min $init $nohes $max $save -j \'$job\'";
	$OUT .= <<__STOP__;
The minimizer will be started with the following commandline:<br>
$indent$cmdline
<br>
<br>
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
	&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="">
<INPUT NAME=prog TYPE=HIDDEN VALUE="$cmdline">
$exetypes
<br>
<INPUT TYPE=SUBMIT VALUE="Start it... ">
</FORM>
__STOP__




} else {

$OUT .=  <<__STOP__;
<FORM ACTION=/exec/min.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=initfiles TYPE=checkbox> No initialization of input files (eg. after a crash)<br>
<INPUT NAME=nohess TYPE=checkbox> Remove .minrestart (initial Hessian for PORT option)<br>
max. number of structure changes <INPUT NAME=maxstruct SIZE=4>
<br>
save every <INPUT NAME=savestep SIZE=4> iterations
<br>
<INPUT NAME=p1 TYPE=CHECKBOX $p>parallel<br>
Job-file (optional): <INPUT NAME=jobfile>
<br>
<INPUT TYPE=hidden NAME=doit VALUE=2>
<INPUT TYPE=SUBMIT VALUE="Prepare commandline">
</FORM>

__STOP__


}
}
PrintPage("Context",$OUT);


