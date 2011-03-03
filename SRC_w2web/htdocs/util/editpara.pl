#!/usr/bin/perl

require "../../libs/w2web.pl";

&GetInput;
&GetSession;

if ( -e "$DIR/.machines" ) {
	redirectURL("/util/edit.pl?SID=$SID&file=$DIR/.machines");
	# good, machines file found.
} else {
	#$OUT .= "NO .machines file found" if $debug;
	redirectURL("/exec/para.pl?SID=$SID");
	exit 0;
}
