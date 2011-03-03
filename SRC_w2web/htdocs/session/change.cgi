#!/usr/bin/perl

$selection = "";
require "../../libs/w2web.pl";
&GetInput;


if ($SID) {
	$OUT.=" Tmp. files of $SID removed" if ($debug);
	system("rm $tempdir/$SID-*");	
}

&listdir($sessionpath);


$OUT .= <<__STOP__;
<CENTER>
<table bgcolor=$topcolor>
<tr><td>
<p>
__STOP__

if ($platinum) {
    $OUT.=<<__STOP__;
<img src="/art/platinum.gif">
__STOP__
} else {
$OUT.=<<__STOP__;
<FONT SIZE=+2>Welcome to <i>w2web</i></FONT><br>
<FONT SIZE=+1>
the fully web-enabled interface to WIEN2k
</FONT>
__STOP__
}
$OUT.=<<__STOP__;
</p>

<table><tr valign=top><td bgcolor=$green>
<H2>Select stored session:</H2>
<form action="/session/change.cgi" METHOD=POST>
<input name=selection value="$selection">
<input type=submit value="show only selection">
</form>

<form action="/session/select.cgi" METHOD=POST>
<SELECT NAME="SID" SIZE=15 >
__STOP__

@newfile = qx( /bin/ls -tr $sessionpath );
$OUT .= "$newfile $#newfile $newfile[$#newfile]";

@sessions = ();

for $file (@ascii_files) {
	$SID=$file;
	&GetSession;
	push(@sessions,"$NAME"."³"."$HOSTNODE"."³"."$file");
}
	
my @mysessions = ();
@mysessions = sort @sessions;



foreach $line (@mysessions) {
	my $sid,$mhost,$name;
	($name, $mhost, $sid) = split (/³/, $line);
	if ( ("$name" =~ m/.*$selection.*/i) || ("$mhost" =~ /.*$selection,*/i) ) {
	$OUT .= "<OPTION VALUE=$sid";
	$OUT .= " SELECTED" if ($sid == $newfile[$#newfile]);
	$OUT .= ">$name";
	if ( $mhost ne "") {
		$OUT .= " ($mhost)\n";
	}
	$OUT .= "\n";
	}
}              


if ($supportnodes) {
	$hostfile="$w2webdir/conf/hosts.conf";
	if (-e $hostfile) {
		$selhost="on host-node<br>";
		$selhost.="<SELECT NAME=HOSTNODE SIZE=10>";
		$selhost.="<OPTION VALUE=\"\" SELECTED>master node";
		open(HOST, $hostfile);
		while(<HOST>) {
			chop;
			$selhost.="<OPTION VALUE=\"$_\">$_" if ($_ !~ /^$/);
		}
		$selhost.="</SELECT>";
	}
	$selhost.=<<__STOP__;
</FORM>
<FORM ACTION="/config/index.pl" METHOD=POST>
<INPUT TYPE=HIDDEN NAME="SID" VALUE=$SID>
<INPUT TYPE=HIDDEN NAME="nexturl" VALUE="/session/change.cgi">
<INPUT TYPE=SUBMIT VALUE="edit hosts">
__STOP__
}
$OUT .= <<__STOP__;    
</SELECT>
<INPUT TYPE=SUBMIT VALUE="Select">
</FORM>
</td><td bgcolor=$gray>
<H2>Create new session:</H2>
<FORM action="/session/new.cgi" METHOD=POST>
<INPUT NAME="NEWNAME" VALUE="Session_name">
<INPUT TYPE=SUBMIT VALUE="Create">
<br>
$selhost
</FORM>

</td></tr></table>
</td></tr></table>
</center>

__STOP__


PrintPage("Change SID", $OUT);
