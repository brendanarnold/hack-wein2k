#!/usr/bin/perl

require "../../libs/w2web.pl";
$debug=0;
&GetInput;
$randmax=999999;

if($FORM{'NEWNAME'}) {
} else {
	&CGIError("You must specify a name for the duplicate session.\n");
}


$newsid = int(rand($randmax));
while ( -e "$sessionpath/$newsid" ) {
	$newsid = int(rand($randmax));
}
 
$OUT .= "new SID = $newsid<br>" if $debug;

$OUT .= <<__STOP__;
<H2>Session duplicated:</H2>
new SID is $newsid

<p>
<A HREF=/index.pl?SID=$newsid TARGET=_parent>click to restart duplicated session</A>
__STOP__


&GetSession;
$SID=$newsid;
$NAME=$FORM{'NEWNAME'};
&SaveSession;

PrintPage("Duplicate SID", $OUT);   
