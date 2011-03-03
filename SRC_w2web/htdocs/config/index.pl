#!/usr/bin/perl

require "../../libs/w2web.pl";
&GetInput;
&GetSession();

$OUT .=  <<__STOP__;
<h2>Interface configuration</h2>
<p>

<FORM ACTION="/config/save.pl" METHOD=POST>
<INPUT TYPE=HIDDEN NAME=nexturl VALUE=$nexturl>
<table><tr><td>
Refresh status every</td><td> <INPUT SIZE=5 NAME="newrefresh" VALUE=$refresh>
seconds</td></tr>
<tr><td>Remove temporary files every</td>
<td> <INPUT SIZE=5 NAME="newremovetime"
VALUE=$removetime> minutes</td></tr>
<tr><td>E-mail for notification </td>
<td><INPUT SIZE=30 NAME="newemail" VALUE=$notification></td></tr>
<tr><td>Show icons for status</td>
<td><INPUT type=checkbox NAME="newicons" 
__STOP__
$OUT .= "CHECKED" if $icons;
$OUT.= <<__STOP__;
></td></tr>
<tr><td>Show javascript menu</td><td><INPUT type=checkbox NAME="newjsmenu" 
__STOP__
$OUT .= "CHECKED" if $jsmenu;
$OUT.= "></td></tr>";
$OUT.= "<tr><td>Always show commandline</td><td><INPUT type=checkbox NAME=\"newcommand\""; 
$OUT .= "CHECKED" if $cmdl;
$OUT.= "></td></tr>";
$OUT.= "<tr><td>Task default in background</td><td><INPUT type=checkbox NAME=\"newtaskback\""; 
$OUT .= "CHECKED" if $taskback;
$OUT.= "></td></tr>";
$OUT.= "<tr><td>ELNES expert mode</td><td><INPUT type=checkbox NAME=\"newelnesexpert\""; 
$OUT .= "CHECKED" if $elnesexpert;
$OUT.= "></td></tr>";
$OUT.="</table>";

if (-e "$w2webdir/conf/execution.conf") {
	$exe=qx(cat $w2webdir/conf/execution.conf);
} else {
	$exe="";
}

$OUT .=  <<__STOP__;
<h3>Execution types</h3>
<p>
Syntax:<br>
Description=command to queue submission<br>
use %f for script file<br>
<TEXTAREA NAME="execution" COLS=30 ROWS=10>$exe</TEXTAREA>
</p>
<p>
__STOP__

if ($supportnodes) {
if (-e "$w2webdir/conf/hosts.conf") {
	$hosts=qx(cat $w2webdir/conf/hosts.conf);
} else {
	$hosts="";
}

$OUT .=  <<__STOP__;
<h3>Host nodes for w2web</h3>
<p>
URL including http(s) and port to host-nodes<br>
<TEXTAREA NAME="hosts" COLS=30 ROWS=10>$hosts</TEXTAREA>
</p>
<p>
<INPUT TYPE=SUBMIT VALUE="Save these values">
</p>
</FORM>
__STOP__
}

PrintPage("Context",$OUT);


