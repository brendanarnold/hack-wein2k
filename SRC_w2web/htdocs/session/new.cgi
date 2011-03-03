#!/usr/bin/perl

require "../../libs/w2web.pl";
&GetInput;

if ($platinum && $NEWNAME =~ /hideplatinum/) {
    $platinum=-1;
    &PlatinumWrite;
    $splash="hide platinum";
    &redirectURL("/session/change.cgi");
    exit;
} 

if (!$platinum && $NEWNAME =~ /platinumedition/) {
    $splash="converting to the<br><br>platinum edition";
    $file="$w2webdir/conf/platinum.conf";
    unless(open(FILE,">$file")) {
	&CGIError("Can't write file $fname.\n");
	exit;
    }
    $platinum=1;
    $platinumcounter=1;
    &PlatinumWrite;
    &redirectURL("/session/change.cgi");
    exit;
}

if ($platinum && $NEWNAME =~ /platinumedition/) {
    $splash="platinum version already activated!";
    &redirectURL("/session/change.cgi");
    exit;
}


#$test = qx(grep "NAME=" $sessionpath/* | grep $NEWNAME | wc -l);
&listdir($sessionpath);
$test=0;
foreach $i (@ascii_files) {
	$test=1 if ($i eq $NEWNAME);
}

if($test) {
	$splash="Sorry, that session name already exits<br>please choose another name";
	&redirectURL("/session/change.cgi");
	exit;
}
$randmax=999999;

if (!$SID) {
	$newsid = int(rand($randmax));
	while ( -e "$sessionpath/$newsid" ) {
		$newsid = int(rand($randmax));
	}
} else {
	$newsid=$SID;
}



$SID="";
&GetSession;
$SID=$newsid;
$NAME=$FORM{'NEWNAME'};
$DIR=$ENV{'HOME'};

&SaveSession;

if ($HOSTNODE && $supportnodes) {
	$splash="redirecting to node<br>$HOSTNODE";
	&redirectURL("$HOSTNODE/session/new.cgi?NEWNAME=$NAME&SID=$SID");
} else {
	&redirectURL("/util/dir.pl?SID=$SID&dir=$DIR&cd=1");
}
