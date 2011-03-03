#!/usr/bin/perl

require "../../libs/w2web.pl";
&GetInput;
&GetSession();



$exe="";
$hosts="";
if (-e "$w2webdir/conf/execution.conf") {
	$exe=qx(cat $w2webdir/conf/execution.conf);
}

if ($supportnodes) {
	if (-e "$w2webdir/conf/hosts.conf") {
		$hosts=qx(cat $w2webdir/conf/hosts.conf);
	}
}

$OUT .=  <<__STOP__;
<h2>Clone settings</h2>

<FORM ACTION="${CLONENODE}/config/save.pl" METHOD=POST>
<INPUT TYPE=HIDDEN NAME=nexturl VALUE="/index.pl">
<INPUT TYPE=HIDDEN NAME="newrefresh" VALUE=$refresh> 
<INPUT TYPE=HIDDEN NAME="newremovetime" VALUE=$removetime>
<INPUT TYPE=HIDDEN NAME="newemail" VALUE=$notification> 
<INPUT TYPE=HIDDEN NAME="execution" value="$exe">
<INPUT TYPE=HIDDEN NAME="hosts" value="$hosts">

<p>
Clone master settings to $CLONENODE
</p>

<INPUT TYPE=SUBMIT VALUE="Save master values"> 
</p>
</FORM>
__STOP__

PrintPage("Context",$OUT);


