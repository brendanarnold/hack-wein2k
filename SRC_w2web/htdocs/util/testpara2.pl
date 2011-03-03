#!/usr/bin/perl

require "../../libs/w2web.pl";

&GetInput;
&GetSession;

$OUT = <<__STOP__;
<h3>Test parallel execution</h3>

<pre>
__STOP__

$umps = qx( cd $DIR; testpara2 );

$OUT .= <<__STOP__;
$umps
</pre>
__STOP__


PrintPage("testpara", $OUT);
