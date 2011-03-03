#!/usr/bin/perl

require "../../libs/w2web.pl";
&GetInput;

$OUT .= "<h2>Session parameters saved.</h2>";
&ShowParms if $debug;
&SaveSession;
$OUT .= <<__STOP__;
<p class="info">
<A HREF="/index.pl?SID=$SID" TARGET=_parent>Click to restart session</A>
</p>
__STOP__


PrintPage("savedir", $OUT);
