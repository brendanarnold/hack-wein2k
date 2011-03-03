#!/usr/bin/perl

require "../../libs/w2web.pl";
$doit=0;

&GetInput;
&GetSession;

if (!$doit) {

$runpid = qx( head -1 $file);
$runpid =~ s/ //;
$runadd = qx( tail -1 $file);

$OUT = <<__STOP__;
<h3>Kill process</h3>

<b>$file</b> exists and should execute $runadd. 
Do you want to kill process ${runpid}?


<FORM ACTION="/util/kill.pl" METHOD=POST>
<INPUT TYPE=HIDDEN NAME="file" VALUE="$file">
<INPUT TYPE=HIDDEN NAME="SID" VALUE="$SID">
<INPUT TYPE=HIDDEN NAME="doit" VALUE="1">
<INPUT TYPE=SUBMIT VALUE="kill this process">
</FORM>

__STOP__

} else {

	#@count=qx(/bin/kill -1 $runpid);
	@count=kill 1, $runpid;
	$umps=qx(rm $file);

	$OUT = <<__STOP__;
<h3>Kill process</h3>

Process <b>$runpid</b> was killed.
<br>
@count processes signalled.
__STOP__
}


PrintPage("Kill process", $OUT);
