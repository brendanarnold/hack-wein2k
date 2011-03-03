#!/usr/bin/perl

require "../../libs/w2web.pl";
$doit=0;

&GetInput;
&GetSession;
&InitTask();

$mytest="$DIR/$CASE.klist_band";
if ( -e $mytest && !$doit) {

$OUT = <<__STOP__;
<h3>Create $CASE.klist_band from $klist template</h3>

<b>$CASE.klist_band</b> already exists !!
Do you want to overwrite it and copy k-list template for <b>$klist</b> to $CASE.klist_band ?
<FORM ACTION="/util/addklist.pl" METHOD=POST>
<INPUT TYPE=HIDDEN NAME="klist" VALUE="$klist">
<INPUT TYPE=HIDDEN NAME="SID" VALUE="$SID">
<INPUT TYPE=HIDDEN NAME="doit" VALUE="1">
<INPUT TYPE=SUBMIT VALUE="Copy it">
</FORM>

__STOP__

} else {

	$from="$WIENROOT/SRC/$klist.klist";
	if ($klist =~  /xcrysden/) {
		$from="$DIR/xcrysden.klist";
	}

#	$umps = qx( cat $from >> $DIR/$CASE.in1$filec );
	$umps = qx( cp $from  $DIR/$CASE.klist_band );
	&redirectURL("/exec/band.pl?SID=$SID");

#	$OUT = <<__STOP__;
#<h3>k-list tempolte</h3>
#
#Template of <b>$klist</b> was appended to <b>$CASE.in1$filec</b>.
#<br>
#$umps
#__STOP__
}


PrintPage("klist-template", $OUT);
