#!/usr/bin/perl

require "../../libs/w2web.pl";
$doit=0;

&GetInput;
&GetSession;

if (!$doit) {

$OUT = <<__STOP__;
<h3>save_lapw</h3>

<FORM ACTION="/util/savelapw.pl" METHOD=POST>
<INPUT TYPE=HIDDEN NAME="SID" VALUE="$SID">
<INPUT TYPE=HIDDEN NAME="doit" VALUE="1">
<INPUT TYPE=CHECKBOX NAME=a> Save all files<br>
<INPUT TYPE=CHECKBOX NAME=f> force save_lapw to overwrite previous saves<br>
<INPUT TYPE=CHECKBOX NAME=d> save calculation in directory specified<br>
<br>
Save name or directory: <br>
$indent <INPUT NAME="savename"><br>
<br>
<INPUT TYPE=SUBMIT VALUE="save">
</FORM>

__STOP__

} else {
	$OUT .= "<h3>save_lapw $savename</h3>";

	if (!$savename) {
		$OUT .= "<b>ERROR</b> no save name given!";
	} else {
		$cmdline = "save_lapw";
		$cmdline .= " -a" if ($a);
		$cmdline .= " -f" if ($f);
		$cmdline .= " -d" if ($d);
		$cmdline .= " $savename";

		$OUT .= "$cmdline";

		$out = qx(cd $DIR; $cmdline);
		$OUT .= "<pre>$out</pre>" ;
	}



}


PrintPage("save_lapw", $OUT);
