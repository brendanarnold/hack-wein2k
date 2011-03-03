#!/usr/bin/perl

require "../../libs/w2web.pl";
$debug=0;
&GetInput;
&GetSession;


if ("http://$ENV{'SERVER_NAME'}:$ENV{'SERVER_PORT'}" !~ $ENV{'MASTER_URL'}) {
	# we must also delete this session from the master
	
	# first remove locally
	qx|rm $sessionpath/$SID|;

	$splash="deleting also on master";
	&redirectURL("$ENV{'MASTER_URL'}/session/delete2.cgi?SID=$SID");

	exit;
}


$OUT .= <<__STOP__;
<H2>Session $SID was deleted.</H2>

<p>
<A HREF=/index.pl TARGET=_parent>click to continue</A>
__STOP__

qx|rm $sessionpath/$SID|;
PrintPage("Delete SID", $OUT);   
