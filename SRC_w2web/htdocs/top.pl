#!/usr/bin/perl

require "../libs/w2web.pl";
&GetInput;
&GetSession;

$OUT .=  <<__STOP__;
Content-type: text/html\n
<HTML>
<HEAD>
   <TITLE>W2web</TITLE>
<link href="$css" rel="stylesheet">
</HEAD>
<BODY BGCOLOR=$topcolor>
Session: <a href="$ENV{'MASTER_URL'}/session/change.cgi?SID=$SID"
target=\"_top\"><b>$NAME</b></a><br/>
$DIR
</BODY>
</HTML>                
__STOP__

print $OUT;
