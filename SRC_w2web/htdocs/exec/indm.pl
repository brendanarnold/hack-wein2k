#!/usr/bin/perl

require "../../libs/w2web.pl";
&GetInput;
&GetSession;
&InitTask();

my $infile = "$WIENROOT/SRC_templates/case.indm$filec";
$umps = qx( cd $DIR;cp  $infile $CASE.indm$filec );

redirectURL("/util/edit.pl?SID=$SID&f=1&file=$DIR/$CASE.indm$filec");


