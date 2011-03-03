#!/usr/bin/perl

require "../../libs/w2web.pl";
&GetInput;
&GetSession;

my $infile = "$WIENROOT/SRC_templates/case.inorb";
$umps = qx( cd $DIR;cp  $infile $CASE.inorb );

redirectURL("/util/edit.pl?SID=$SID&f=1&file=$DIR/$CASE.inorb");


