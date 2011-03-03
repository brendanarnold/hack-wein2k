#!/usr/bin/perl

require "../../libs/w2web.pl";
$doit=0;

&GetInput;
&GetSession;

if (!$doit) {

$OUT = <<__STOP__;
<h3>Download file</h3>

Download <b>$file</b> 
<FORM ACTION="/util/download.pl" METHOD=POST>
<INPUT TYPE=HIDDEN NAME="file" VALUE="$file">
<INPUT TYPE=HIDDEN NAME="SID" VALUE="$SID">
<INPUT TYPE=HIDDEN NAME="doit" VALUE="1">
<!--
<INPUT NAME="type" TYPE=radio value="ascii" selected>as text
<br>
<INPUT NAME="type" TYPE=radio value="binary">as binary file
-->
<br>
<INPUT TYPE=SUBMIT VALUE="download it">
</FORM>

__STOP__

} else {


	$umps = qx( cat $file );

	$OUT = <<__STOP__;
<h3>downloading file</h3>

<b>$file</b> 
<br>
__STOP__
	print "Content-Disposition: attachment; filename=\"$file\"\n";
  print "Content-type: application/x-octet-stream\n\n";
	open(FILE,$file);
	print <FILE>;
	close;
	exit;


}


PrintPage("Download file", $OUT);
