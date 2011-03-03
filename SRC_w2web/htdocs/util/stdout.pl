#!/usr/bin/perl

require "../../libs/w2web.pl";

&GetInput;
&GetSession;

$OUT = <<__STOP__;
<h3>Show STDOUT</h3>

<table bgcolor=$green>
<tr><td bgcolor=$gray>
<FORM ACTION=/util/stdout.pl METHOD=POST>
<INPUT TYPE=HIDDEN NAME=SID VALUE=$SID>
<INPUT TYPE=SUBMIT VALUE="reload">
</FORM>
</td><td bgcolor=$gray>
<FORM ACTION=/util/stdout.pl METHOD=POST>
<INPUT TYPE=HIDDEN NAME=SID VALUE=$SID>
<INPUT TYPE=HIDDEN NAME=reverse VALUE=1>
<INPUT TYPE=SUBMIT VALUE="reload in reverse order">
</FORM>
</td></tr>
</table>
<pre>
__STOP__

if ($reverse) {
	$umps = qx( tac $DIR/STDOUT );
} else {
	$umps = qx( cat $DIR/STDOUT );
}
#$umps =~ s/\n/<br>/;

$OUT .= <<__STOP__;
$umps
</pre>
__STOP__


PrintPage("Delete file", $OUT);
