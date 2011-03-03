#!/usr/bin/perl

require "../../libs/w2web.pl";
$doit=0;

&GetInput;
&GetSession;

if (!$doit) {

$OUT = <<__STOP__;
<h3>clean_lapw</h3>

<FORM ACTION="/util/cleanlapw.pl" METHOD=POST>
<INPUT TYPE=HIDDEN NAME="SID" VALUE="$SID">
<INPUT TYPE=HIDDEN NAME="doit" VALUE="1">
<br>
Deletes all large and all files that are not
necessary for continuation?
<br>
<INPUT TYPE=SUBMIT VALUE="yes, clean directory">
</FORM>

__STOP__

} else {
	$OUT .= "<h3>clean_lapw</h3>";

	$out = qx(cd $DIR; echo "y" | clean_lapw);
	$OUT .= "<pre>$out</pre>" ;

}


PrintPage("save_lapw", $OUT);
