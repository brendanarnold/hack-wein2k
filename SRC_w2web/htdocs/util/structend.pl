#!/usr/bin/perl


require "../../libs/w2web.pl";
require "../../libs/struct.pl";
$debug=0;
&GetInput;
&GetSession;
&StructRead;

$umps = qx(cp $DIR/$CASE.struct_i $DIR/$CASE.struct_ii);
$umps = qx(mv $DIR/$CASE.struct_i $DIR/$CASE.struct);

&InstGen;

redirectURL("/util/structgen.pl?SID=$SID");


