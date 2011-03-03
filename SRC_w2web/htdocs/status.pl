#!/usr/bin/perl

$manual=0;

require "../libs/w2web.pl";
&GetInput;
&GetSession;


if ($manual) {
  $refresh=-1;
  $icons=0;
}


$tm = localtime;

use POSIX qw(strftime);
$tm = strftime("%H:%M:%S", localtime);


$flag=1;
$semfile="$tempdir/$SID-$$.find";


if ( -e $semfile ) {
	$flag=0;
}

system "touch $semfile";

@result = glob("$DIR/.running*");
$run = "idle";
$run = "<img src=\"/art/sleep.gif\" hspace=0>" if $icons;
foreach $i (@result) {
    if ( -e $i ) {
	$running = "running";
	$running = "stop SCF" if ( -e "$DIR/.stop" ) ;
	$running = "stop MINI" if ( -e "$DIR/.minstop" ) ;
	$running .= "|fulldiag" if ( -e "$DIR/.fulldiag" ) ;
	$run = "<b><a href=\"/util/dir.pl?SID=$SID&dir=$DIR&ext=run\" target=main>$running</a></b>";
	$runadd = qx( cat $i | tail -1 );
	$runadd =~ s/ *//;
	$run .= " ($runadd)";
	$run .= "<img src=\"/art/run-ani.gif\">" if $icons;
	$topcolor=$red;
    }
}
@result = glob("$DIR/*.error");
foreach $i (@result) {
    if ( -z $i ) {	
	# this is good, no error
    } else {
	$test = qx( exec wc -l $i );
	if ($test > 1) {
	    $run = "<a href=\"/util/dir.pl?SID=$SID&dir=$DIR&ext=error\" target=main>error</a>";
	    $run = "<img src=\"/art/stupid2.gif\"> <a href=\"/util/dir.pl?SID=$SID&dir=$DIR&ext=error\" target=main>error</a>" if $icons;
	    $topcolor=$darkred;
	}
    }
}

# -----------------------------------------
# check for temporary files and clean up
# -----------------------------------------
$umps = qx(find $tempdir -cmin +$removetime -print ) if $flag;
$count=0;
foreach $i (split(/\n/,$umps)) {
    $count++;
    if ( -f $i) {
	system "rm $i";
    }
}
# -----------------------------------------

$run .= " (skipping find)" if !($flag);

# -----------------------------------------
# time-stamping session file
# -----------------------------------------
system "touch $sessionpath/$SID";
# -----------------------------------------

$renew = <<__STOP__;
<META HTTP-EQUIV="refresh" content="$refresh; URL=status.pl?SID=$SID">
<META HTTP-EQUIV="PRAGMA" CONTENT="no-cache">
<META HTTP-EQUIV="EXPIRES" CONTENT="now">
__STOP__

$OUT .=  <<__STOP__;
Content-type: text/html\n
<html>
<head>
   <title>WIEN2k</title>
$renew
</head>                                    
<link href="$css" rel="stylesheet">
</head>
<body bgcolor=$topcolor>
<div align=right>
<font size=-2>
<nobr>$tm $run</nobr><br>
__STOP__

if ($manual) {
    $OUT .=  <<__STOP__;
<a href="status.pl?SID=$SID&manual=1">refresh now</a>
|
<a hreF="status.pl?SID=$SID&manual=0">auto-refresh</a>
__STOP__
} else {
    $OUT .=  <<__STOP__;
<a href="status.pl?SID=$SID">refresh</a>
|
<a href="status.pl?SID=$SID&manual=1">no refresh</a>
__STOP__
}

$OUT .=  <<__STOP__;
</font>
</div>
</body>
</html>                          
__STOP__

print $OUT;

system "rm $semfile" if $flag;
