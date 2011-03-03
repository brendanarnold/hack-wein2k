#!/usr/bin/perl

require "../../libs/w2web.pl";
$doit=0;

&GetInput;
&GetSession;

if (!$doit) {

$OUT = <<__STOP__;
<h3>New file</h3>

<FORM ACTION="/util/newfile.pl" METHOD=POST>
<INPUT TYPE=HIDDEN NAME="file" VALUE="$file">
<INPUT TYPE=HIDDEN NAME="SID" VALUE="$SID">
<INPUT TYPE=HIDDEN NAME="doit" VALUE="1">
filename <INPUT NAME="newfile" VALUE="$CASE"><br>
<br>
<INPUT TYPE=SUBMIT VALUE="create">
</FORM>

__STOP__

} else {

	$umps = qx( touch  $DIR/$newfile );
	redirectURL("/util/edit.pl?SID=$SID&file=$DIR/$newfile");
	exit;
}


PrintPage("New file", $OUT);
