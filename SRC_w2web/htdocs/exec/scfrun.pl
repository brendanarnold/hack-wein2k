#!/usr/bin/perl

require "../../libs/w2web.pl";
$spinpol=0;

@options = qw/h it0 p so dm renorm in1orig orb/; 
$prefspace="scf";

$debug=0;

&GetInput();
$dum=$spinpol;

&InitTask();

&GetSession();
$spinpol=$dum;

if ($saveonly) {
	&SavePrefs();
	redirectURL("/exec/scf.pl?SID=$SID&doit=1");
	exit;
}


$opts = "";
# check prerequired files and generate them if necessary
#system "echo '$p' $FORM{'p1'} > /tmp/mistfile";
	$p="";
if($FORM{'p1'}) {
	# oh dear, parallel execution!
	$p="on";
	$opts .= " -p ";
        if ( -e "$DIR/.machines" ) {
		# good, machines file found.
		$OUT .= ".machines file found" if $debug;
	} else {
		&SaveSession;
		#$OUT .= "NO .machines file found" if $debug;
		redirectURL("/exec/para.pl?SID=$SID");
		exit 0;
	}
}
if($FORM{'so'}) {
	if ( -e "$DIR/$CASE.inso" ) {
	   if ( -z "$DIR/$CASE.inso" ) {
                # file is 0 byte long
		#$OUT .= "NO .inso file found" if $debug;
		redirectURL("/util/initso.pl?SID=$SID");
		exit 0;
	}
		# good, inso file found.
		$OUT .= ".inso file found" if $debug;
	} else {
		#$OUT .= "NO .inso file found" if $debug;
		redirectURL("/util/initso.pl?SID=$SID");
		exit 0;
	}
}
if($FORM{'dm'}) {
	if ( -e "$DIR/$CASE.indm$filec" && -s "$DIR/$CASE.indm$filec" ) {
		# good, indm$filec file found.
		$OUT .= ".indm$filec file found" if $debug;
	} else {
		#$OUT .= "NO .indm$filec file found" if $debug;
		redirectURL("/exec/indm.pl?SID=$SID");
		exit 0;
	}
}
if($FORM{'orb'}) {
	if ( -e "$DIR/$CASE.inorb" && -s "$DIR/$CASE.inorb"  ) {
		# good, inorb file found.
		$OUT .= ".inorb file found" if $debug;
	} else {
		#$OUT .= "NO .inorb file found" if $debug;
		redirectURL("/exec/inorb.pl?SID=$SID");
		exit 0;
	}
	if ( -e "$DIR/$CASE.indm$filec" && -s "$DIR/$CASE.indm$filec" ) {
		# good, indm$filec file found.
		$OUT .= ".indm$filec file found" if $debug;
	} else {
		#$OUT .= "NO .indm$filec file found" if $debug;
		redirectURL("/exec/indm.pl?SID=$SID");
		exit 0;
	}
}
# create commandline
$cmd = "run_lapw";
if ("$FORM{'spinpol'}" eq "on") {
	$cmd = "runsp_lapw";
}

	$afm="";
if ("$FORM{'afm1'}" eq "on") {
	$cmd = "runafm_lapw";
	$afm="on";
}

if ("$FORM{'fsm'}" eq "on") {
	$cmd = "runfsm_lapw";
}

&SaveSession;

foreach $i (@options) {
	if ("$FORM{$i}" eq "on") {
		$opts .= " -$i ";
	}
}


if ($FORM{'fsm'}) {
	$opts .= " -m $FORM{'fsm_count'}";
}
if ($FORM{'nohns'}) {
	$opts .= " -nohns $FORM{'nohns_count'}";
}
if ($FORM{'it'}) {
	$opts .= " -it $FORM{'it_count'}";
}
if ($FORM{'in1new'}) {
	$opts .= " -in1new $FORM{'in1new_count'}";
}
if ($FORM{'ql'}) {
	$opts .= " -ql $FORM{'ql_count'}";
}
if ($FORM{'itnum'}) {
	$opts .= " -i $FORM{'itnum_count'}";
}
if ($FORM{'conv_ec'}) {
	$opts .= " -ec $convval_ec";
}
if ($FORM{'conv_fc'}) {
	$opts .= " -fc $convval_fc";
}
if ($FORM{'conv_cc'}) {
	$opts .= " -cc $convval_cc";
}
#if ($FORM{'conv'}) {
#    $convval_ = "$FORM{'convval_ec'}";
#    if ( $FORM{'conv'} =~ 'ec' ) { $convval_ = "$FORM{'convval_ec'}"; }
#    if ( $FORM{'conv'} =~ 'fc' ) { $convval_ = "$FORM{'convval_fc'}"; }
#    if ( $FORM{'conv'} =~ 'cc' ) { $convval_ = "$FORM{'convval_cc'}"; }
#	$opts .= " -$FORM{'conv'} $convval_";
#}

$cmd .= $opts;

if ($FORM{'expert'}) {
	$cmd .= " $FORM{'expert'}";
}

$cmd .= " -NI";

&SavePrefs();

$OUT .= <<__STOP__;
<P><B>Commandline: $cmd</B>
<br>
has been sent to system for execution</P>
<HR>
__STOP__

$sendmail="";
if ($notify) {
	$sendmail="echo \"$NAME finished.\" | mail $notification";
}

if ($exetype =~ /interactive/ ) {
} else {
$cmdfile="$DIR/.command.$SID.$$";
unless(open(FILE,">$cmdfile")) {
	&CGIError("Can't write file $fname.\n");
	exit;
}

print FILE <<__STOP__;
#!/bin/sh
date
cd $DIR
$cmd >$DIR/STDOUT 2>&1
$sendmail
rm $cmdfile
__STOP__
close(FILE);
$umps=qx(chmod +x $cmdfile);
}

if ($exetype =~ /background/ ) {
	system "$cmdfile &";
} elsif ($exetype =~ /interactive/ ) {
		$OUT.= "cd $DIR;$cmd \n<pre>";
		$umps = qx( cd $DIR;$cmd 2>&1);
		$OUT .= "$umps\n";
		$OUT .= "</pre>";
} else {
	if ($exetype =~ /%f/) {
		$exetype=~ s/%f/ $cmdfile/;
		$OUT.="we did: cd $DIR;$exetype &";
		system "cd $DIR;$exetype &";
	} else {
		$OUT.= "cd $DIR;$exetype $cmd >$DIR/STDOUT 2>&1 &";
		system "cd $DIR;$exetype $cmd >$DIR/STDOUT 2>&1 &";
		system "cd $DIR;rm $cmdfile";
	}
}

if ($exetype =~ /interactive/ ) {
} else {
$OUT .= <<__STOP__;
<p>
  <A HREF="/util/stdout.pl?SID=$SID">View STDOUT</A>
</p>
__STOP__
}

PrintPage("Context",$OUT);


