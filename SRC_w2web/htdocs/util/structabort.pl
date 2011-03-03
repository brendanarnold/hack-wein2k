#!/usr/bin/perl


require "../../libs/w2web.pl";
$debug=0;
&GetInput;
&GetSession;

$umps = qx(rm $DIR/$CASE.struct_i);


redirectURL("/util/structgen.pl?SID=$SID");


