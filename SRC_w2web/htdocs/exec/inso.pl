#!/usr/bin/perl

require "../../libs/w2web.pl";
&GetInput;
&GetSession;

$next="continue with initso";
$nexturl="/exec/inso.pl?SID=$SID";

if ( -e "$CASE.inso" ) {
$OUT = <<__STOP__;
<p> inso file present 
</p>
__STOP__

   } else {
my $infile = "$WIENROOT/SRC_templates/case.inso";
$umps = qx( cd $DIR;cp  $infile $CASE.inso );

redirectURL("/util/edit.pl?SID=$SID&f=1&file=$DIR/$CASE.inso&next=$next&nexturl=$nexturl");
    } 


PrintPage("Context",$OUT);
