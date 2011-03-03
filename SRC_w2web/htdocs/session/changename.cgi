#!/usr/bin/perl

require "../../libs/w2web.pl";
&GetInput;
&GetSession;



$OUT .= "<h2>Change Session information ($SID)</h2>";


$COMMENT =~ s/<br>/\n/g;

$OUT .= <<__STOP__;
<FORM ACTION=/session/save2.cgi METHOD=POST>
<INPUT NAME=SID VALUE=$SID TYPE=HIDDEN>
<INPUT NAME=DIR VALUE=$DIR TYPE=HIDDEN>
<INPUT NAME=ALERT VALUE=$SID TYPE=HIDDEN>
<br>
Session name:
<INPUT NAME=NAME VALUE="$NAME">
<br>
<br>
Comments:
<br>
<TEXTAREA NAME="COMMENT" WRAP=VIRTUAL COLS=40 ROWS=10>
$COMMENT
</TEXTAREA>
<br>
<br>
<INPUT TYPE=checkbox NAME="spinpol" $spinpol>spin polarized calculation<br>
<INPUT TYPE=checkbox NAME="afm" $afm>AFM calculation<br>
<INPUT TYPE=checkbox NAME="complex" $complex>complex calculation (no inversion)<br>
<INPUT TYPE=checkbox NAME="p" $p>k-parallel execution<br>
<br>
<INPUT TYPE=SUBMIT VALUE="save changes">
</FORM>
__STOP__


PrintPage("savedir", $OUT);
