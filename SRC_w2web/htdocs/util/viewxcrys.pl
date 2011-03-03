#!/usr/bin/perl

require "../../libs/w2web.pl";

&GetInput;
&GetSession;

$OUT = <<__STOP__;
<h3>View structure with XCrysDen</h3>

Requires X-Windows system ...
__STOP__

# $ENV{'DISPLAY'}="$ENV{'REMOTE_HOST'}:0.0";
# $DISPLAY=":0";
$umps = qx(echo $DISPLAY && cd $DIR && xcrysden --wien_struct $DIR  &);

$OUT .= <<__STOP__;
<pre>
$umps
</pre>
__STOP__


PrintPage("xcrysden", $OUT);
