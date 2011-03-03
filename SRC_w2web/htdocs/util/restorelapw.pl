#!/usr/bin/perl

require "../../libs/w2web.pl";
$doit=0;

&GetInput;
&GetSession;

if (!$doit) {

$OUT = <<__STOP__;
<h3>restore_lapw</h3>

<FORM ACTION="/util/restorelapw.pl" METHOD=POST>
<INPUT TYPE=HIDDEN NAME="SID" VALUE="$SID">
<INPUT TYPE=HIDDEN NAME="doit" VALUE="1">
<INPUT TYPE=CHECKBOX NAME=f> force restore_lapw to overwrite previous files<br>
<INPUT TYPE=CHECKBOX NAME=d> restore calculation from directory specified<br>
<INPUT TYPE=CHECKBOX NAME=t> only test which files would be restored<br>
<br>
Save name or directory: <br>
$indent <INPUT NAME="savename"><br>
<br>
<INPUT TYPE=SUBMIT VALUE="restore">
</FORM>

__STOP__

} else {
	$OUT .= "<h3>restore_lapw $savename</h3>";

	if (!$savename) {
		$OUT .= "<b>ERROR</b> no save name given!";
	} else {
		$cmdline = "restore_lapw";
		$cmdline .= " -t" if ($t);
		$cmdline .= " -f" if ($f);
		$cmdline .= " -d $savename" if ($d);
		$savename =~ s/^.*\/// ;
		$cmdline .= " $savename";
		$OUT .= "$cmdline";

		$out = qx(cd $DIR; $cmdline);
		$OUT .= "<pre>$out</pre>" ;
	}



}


PrintPage("save_lapw", $OUT);
