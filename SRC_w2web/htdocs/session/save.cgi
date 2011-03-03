#!/usr/bin/perl

require "../../libs/w2web.pl";
&GetInput;
&GetSession;


$OUT .= "<h2>Saving session parameters: </h2>";

if (-d $dir && -x $dir && -w $dir ) {
    if ($ENV{"HOME"} =~ $dir ) {$OUT .= <<__STOP__;
    Selected directory is your HOME directory ($dir) ! <BR>
    <BR><A HREF="/util/dir.pl?SID=$SID&dir=$dir&cd=1" TARGET="main">Go back and select a different directory!</A> <BR><BR>Save aborted.
__STOP__
    } else {
	$OUT .= "<p>$dir is writeable directory  </p>" if $debug;
	$DIR = $dir;
	&ShowParms if $debug;
	&SaveSession;
 	$OUT .= <<__STOP__;
<p class="info">
<A HREF="/index.pl?SID=$SID" TARGET=_parent>Click to restart session</A>
</p>
__STOP__
 }
} else { 
	$OUT .= "$dir is not a writeable directory!  <BR><BR>Save aborted.";
}
PrintPage("savedir", $OUT);
