#!/usr/bin/perl

require "../../libs/w2web.pl";

&GetInput;
&GetSession;

$OUT = <<__STOP__;
<h3>Electron Density with XCrysDen</h3>

Requires X-Windows system ...
__STOP__

#$ENV{'DISPLAY'}="$ENV{'REMOTE_HOST'}:0.0";

if ("$wos_tua_ma" =~ /calc/) {
$umps = qx(cd $DIR && xcrysden --wien_density $DIR  &);
$OUT .="Calc";
} else {
$umps = qx(cd $DIR && xcrysden --wien_renderdensity $DIR  &);
$OUT .="Render";

}


PrintPage("xcrysden", $OUT);
