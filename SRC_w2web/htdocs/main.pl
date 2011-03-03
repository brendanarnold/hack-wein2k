#!/usr/bin/perl

require "../libs/w2web.pl";
&GetInput;
&GetSession;

$OUT .=  <<__STOP__;
<H3>w2web, the fully web-enabled interface to WIEN2k</H3>

<table bgcolor=$gray border=1>
<tr valign=top>
<td><b>Session Name:</b></td><td> $NAME</td></tr>
<tr valign=top><td><b>Session ID:</b> </td><td>$SID
<br>
<tr valign=top><td><b>Directory:</b> </td><td>$DIR
<br>
<tr valign=top><td><b>Last changed:</b> </td><td>$TIME
<br>
<tr valign=top><td><b>Comments:</b></td><td>
$COMMENT
</td></tr>

<tr><td colspan=2>
<form action=/session/changename.cgi>
<input type=hidden name=SID value=$SID>
<INPUT TYPE=checkbox NAME="spinpol" $spinpol>spin polarized calculation
<br>
<INPUT TYPE=checkbox NAME="afm" $afm>AFM calculation
<br>
<INPUT TYPE=checkbox NAME="complex" $complex>complex calculation (no inversion)
<br>
<INPUT TYPE=checkbox NAME="p" $p>parallel calculation
</td></tr>
</table>
<br>

<input type=submit value="Change session information">

</FORM>
__STOP__

PrintPage("Context",$OUT);

