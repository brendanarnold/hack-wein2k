#!/usr/bin/perl

require "../../libs/w2web.pl";

&GetInput;
&GetSession;

$OUT = <<__STOP__;
<h3>Stop SCF</h3>

<p>
SCF cycle will stop after completion of mixer.
</p>

<p>
No, I don't want that, <A HREF="/util/delete.pl?SID=$SID&file=$DIR/.stop">remove file .stop</A>.
</p>


__STOP__

$umps = qx( cd $DIR; touch .stop );

$OUT .= <<__STOP__;
$umps
</pre>
__STOP__


PrintPage("testpara", $OUT);
