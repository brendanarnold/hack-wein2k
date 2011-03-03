#!/usr/bin/perl

require "../../libs/w2web.pl";

&GetInput;
&GetSession;

$OUT = <<__STOP__;
<h3>FullDiag</h3>

<p>
Initiate a full diagonalization in the next run of LAPW1.
</p>

<p>
No, I don't want that, <A HREF="/util/delete.pl?SID=$SID&file=$DIR/.fulldiag">remove file .fulldiag</A>.
</p>

__STOP__

$umps = qx( cd $DIR; touch .fulldiag );

$OUT .= <<__STOP__;
$umps
</pre>
__STOP__


PrintPage("testpara", $OUT);
