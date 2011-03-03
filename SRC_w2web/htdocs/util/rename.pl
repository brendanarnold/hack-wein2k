#!/usr/bin/perl

require "../../libs/w2web.pl";
$doit=0;

&GetInput;
&GetSession;

if (!$doit) {

$OUT = <<__STOP__;
<h3>Rename file</h3>

Rename <b>$file</b> 
<FORM ACTION="/util/rename.pl" METHOD=POST>
<INPUT TYPE=HIDDEN NAME="file" VALUE="$file">
<INPUT TYPE=HIDDEN NAME="SID" VALUE="$SID">
<INPUT TYPE=HIDDEN NAME="doit" VALUE="1">
to <INPUT NAME="newfile"><br>
<br>
<INPUT TYPE=SUBMIT VALUE="rename it">
</FORM>

__STOP__

} else {

	$umps = qx( mv $file $DIR/$newfile );

	$OUT = <<__STOP__;
<h3>Rename file</h3>

<b>$file</b> was renamed to <b>$DIR/$newfile</b>.
<br>
$umps
__STOP__
}


PrintPage("Delete file", $OUT);
