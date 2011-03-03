#!/usr/bin/perl

require "../../libs/w2web.pl";
&GetInput;

if ($SID != "") {
	&GetSession;
	if ( -e $DIR ) {
	    &PlatinumWrite;
	    &redirectURL("../index.pl?SID=$SID");
	} else {
		# directory does not exist any more
		$splash="$DIR does not exist any more!";
		&redirectURL("/session/delete.cgi?SID=$SID&SDEL=1");
	}
} else {
	&redirectURL('../index.pl');
}
