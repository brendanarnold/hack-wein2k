#!/usr/bin/perl

require "../../libs/w2web.pl";

$in="0.1";
$cal="+0.0+1/3";
$ump = "\$erg = $in + $cal";
eval $ump;


$OUT = <<__STOP__;
in: $in\n
cal: $cal\n
<br>
erg: $erg\n

__STOP__


PrintPage("testpara", $OUT);
