#!/usr/bin/perl

require "../../libs/w2web.pl";
require "../../libs/struct.pl";
require "../../libs/telnes2.pl";

#$debug=1;
&GetInput;
&GetSession;
&InitTask();
&ShowParms if $debug;
&StructRead;
&InnesRead;


$mycmd = "x lapw2 -qtl $mycomplex $myspinopt $mypara ";

$OUT.= <<__STOP__;
<h2>Modified lapw2 -qtl for TELNES.2</h2>
__STOP__

if ($t_orient) {
    $OUT.="<p>calculation is orientation sensitive -> we need ISPLIT=99</p>";

    $OUT .= "<ul>";
    $umps = qx( cd $DIR; cp $CASE.struct $CASE.struct_telnes);
    $OUT .= "<li>backing up $CASE.struct</li>";

    # set isplit
    $s_isplit[$t_atom]=99;
    $OUT .= "<li>set ISPLIT=99 for atom $t_atom</li>";

    &StructWrite;
    $umps = qx( cd $DIR; mv $CASE.struct_i $CASE.struct);
    $OUT .= "<li>write modified $CASE.struct</li>";

    if ($nextinteractive) {
	$OUT .= "<li><b>command:</b> $mycmd</p>";
	$umps .= qx(cd $DIR; echo '$mycmd'; $mycmd 2>&1);
	$OUT .= "<p><b>output</b>:<br><pre>$umps</pre></p>\n";
	$OUT .= "</li>";
	system "cd $DIR; mv $CASE.struct_telnes $CASE.struct";
	$OUT .= "<li>restore $CASE.struct</li>";
    } else {
	# background
	system "(cd $DIR;$mycmd;mv $CASE.struct_telnes $CASE.struct) >$DIR/STDOUT 2>&1 &";
	$OUT .= <<__STOP__;
<li>starting $mycmd in background
</ul>

<p class="info"><b>INFO:</b><br>
While TELNES.2 is running there is a modified $CASE.struct file 
in your directory! This will be reset to the original file 
automatically after the calculation has finished.
</p>

<p>	   
	<A HREF="/util/stdout.pl?SID=$SID">View STDOUT</A>
</p>
__STOP__
    }

} else {
    #$OUT.="<p>calculation is not orientation sensitive -> regular run</p>";
    if ($nextinteractive) {
	$OUT.="<p><b>command:</b> $mycmd</p>";
	$umps = qx( cd $DIR; echo "$mycmd"; $mycmd );
	$OUT .= "<p><b>output</b>:<br><pre>$umps</pre></p>\n";
    } else {
	# background
	system "cd $DIR;$mycmd >$DIR/STDOUT 2>&1 &";
	$OUT .= <<__STOP__;
<p>starting $mycmd in background</p>

<p>	   
	<A HREF="/util/stdout.pl?SID=$SID">View STDOUT</A>
</p>
__STOP__
    }
}

    


if($next) {
    $OUT .= <<__STOP__;
<FORM ACTION="/exec/next.pl" METHOD=POST>
<INPUT TYPE=HIDDEN NAME=nexturl VALUE="$nexturl">
<INPUT TYPE=submit VALUE="$FORM{'next'}">
</FORM>
__STOP__
}


PrintPage("Context",$OUT);


