#!/usr/bin/perl

require "../../libs/w2web.pl";

&GetInput;
&GetSession;

$OUT = <<__STOP__;
<h3>Stop MINI</h3>

<p>
MINI will stop after completion of the current SCF cycle.
</p>

<p>
No, I don't want that, <A HREF="/util/delete.pl?SID=$SID&file=$DIR/.minstop">remove file .minstop</A>.
</p>



__STOP__

$umps = qx( cd $DIR; touch .minstop );

$OUT .= <<__STOP__;
$umps
</pre>
__STOP__


PrintPage("testpara", $OUT);
