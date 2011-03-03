#!/usr/bin/perl

require "../../libs/w2web.pl";
&GetInput;
&GetSession;

my $infile = "$WIENROOT/SRC_templates/.machines";
$umps = qx( cd $DIR;cp  $infile . );

redirectURL("/util/edit.pl?SID=$SID&f=1&file=$DIR/.machines");


