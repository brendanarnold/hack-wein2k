#!/usr/bin/perl

require "../libs/w2web.pl";
&GetInput;

if (!($SID)) {
	if ($ENV{'MASTER_URL'} && $supportnodes) {
		$splash="redirecting to master on <br>$ENV{'MASTER_URL'}";
	}
	&redirectURL("$ENV{'MASTER_URL'}/session/change.cgi");
} else {

	&GetSession();

	# this session runs on a host
	if($HOSTNODE && $supportnodes) {
		$splash="redirecting to<br>$HOSTNODE";
		&redirectURL("$HOSTNODE/index.pl?SID=$SID");
		exit;
	}

	# oops, we're still in the home directory
	if ($ENV{"HOME"} =~ $DIR ) {
		&redirectURL("/util/dir.pl?SID=$SID&dir=$DIR&cd=1");
		exit;
	}

$OUT .=  <<__STOP__;
Content-type: text/html\n
<HTML>
<HEAD>
   <TITLE>$NAME\@$w2webhost</TITLE>
<link href="/w2web.css" rel="stylesheet">
</HEAD>

<FRAMESET COLS=155,* BORDER=0 FRAMEBORDER=0 FRAMESPACING=0 >
<FRAME SRC="navig.pl?SID=$SID" NAME=left MARGINWIDTH=5 MARGINHEIGHT=0
FRAMEBORDER=0 scrolling=no>
__STOP__

if ($cmdl) {
	$OUT.="<FRAMESET ROWS=\"50,*,40\" BORDER=1 FRAMEBORDER=0 framespacing=0>";
} else {
	$OUT.="<FRAMESET ROWS=\"50,*\" BORDER=1 FRAMEBORDER=0 framespacing=0>";
}

$OUT.=<<__STOP__;
		<FRAMESET COLS=*,200 BORDER=0 FRAMEBORDER=0 framespacing=0>
			<FRAME SRC="top.pl?SID=$SID" NAME=oben SCROLLING=no MARGINWIDTH=5 MARGINHEIGHT=0 FRAMEBORDER=0>
			<FRAME SRC="status.pl?SID=$SID" NAME=status SCROLLING=no MARGINWIDTH=5 MARGINHEIGHT=0 FRAMEBORDER=0>
</FRAMESET>
__STOP__

if ($cd == 1) {
	$OUT .= <<__STOP__;
		<FRAME SRC="/util/dir.pl?SID=$SID&dir=$DIR&cd=1" NAME=main MARGINWIDTH=0 MARGINHEIGHT=0 FRAMEBORDER=0>
__STOP__
} else {
	$OUT .= <<__STOP__;
		<FRAME SRC="main.pl?SID=$SID" NAME=main MARGINWIDTH=0 MARGINHEIGHT=0 FRAMEBORDER=0>
__STOP__

	if ($cmdl) {
		$OUT.="
			<FRAME SRC=\"/exec/cmdline2.pl?SID=$SID\" SCROLLING=no MARGINWIDTH=5 MARGINHEIGHT=0 FRAMEBORDER=0>
";
	}
}

$OUT .= <<__STOP__;
</FRAMESET>
</FRAMESET>
__STOP__


$OUT.=<<__STOP__;
<BODY>
<NOFRAMES>
<P>w2web</P>
<P>Sorry, but for w2web you need a frame-enabled web browser.
</NOFRAMES>
</HTML>                           
__STOP__

print $OUT;

}
