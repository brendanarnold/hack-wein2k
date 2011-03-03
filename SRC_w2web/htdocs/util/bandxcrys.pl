#!/usr/bin/perl

require "../../libs/w2web.pl";

&GetInput;
&GetSession;

$OUT = <<__STOP__;
<h3>K-path with XCrysDen</h3>

Requires X-Windows system ...
__STOP__

#$ENV{'DISPLAY'}="$ENV{'REMOTE_HOST'}:0.0";
$umps = qx(cd $DIR && xcrysden --wien_kpath $DIR  &);


PrintPage("xcrysden", $OUT);
