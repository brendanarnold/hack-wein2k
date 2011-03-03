#!/usr/bin/perl

require "../../libs/w2web.pl";
&GetInput;
&GetSession();

chomp($execution);
$execution=~ s///g;
chomp($hosts);
$hosts=~ s///g;

$OUT .=  <<__STOP__;
<h2>Interface configuration</h2>
<p>
Refresh status every <b>$newrefresh</b> seconds<br>
Remove temporary files every <b>$newremovetime</b> minutes<br>
E-mail for notification to <b>$newemail</b> <br>
Show icons <b>$newicons</b> <br>
Show jsmenu <b>$newjsmenu</b> <br>
Always show commandline <b>$newcommand</b> <br>
Task default in background <b>$newtaskback</b> <br>
ELNES expert mode <b>$newelnesexpert</b>
</p>
<h3>Execution types</h3>
<pre>$execution</pre>
<h3>remote nodes</h3>
<pre>$hosts</pre>
__STOP__

# write user.conf
$file="$w2webdir/conf/user.conf";
unless(open(FILE,">$file")) {
    &CGIError("Can't write file $fname.\n");
    exit;
}
print FILE <<__STOP__;
refresh=$newrefresh
removetime=$newremovetime
notification=$newemail
icons=$newicons
jsmenu=$newjsmenu
cmdl=$newcommand
taskback=$newtaskback
elnesexpert=$newelnesexpert
__STOP__
close(FILE);

# write execution.conf

$file="$w2webdir/conf/execution.conf";
unless(open(FILE,">$file")) {
    &CGIError("Can't write file $fname.\n");
    exit;
}
print FILE <<__STOP__;
$execution
__STOP__
close(FILE);

# write hosts.conf
if ($supportnodes ) {
    $file="$w2webdir/conf/hosts.conf";
    unless(open(FILE,">$file")) {
	&CGIError("Can't write file $fname.\n");
	exit;
    }
    print FILE <<__STOP__;
    $hosts
__STOP__
    close(FILE);
}

# clone config:
# only available on master

if ($supportnodes && ! $ENV{'MASTER_URL'}) {
    $hostfile="${w2webdir}/conf/hosts.conf";
    if (-e $hostfile) {
	$OUT.= <<__STOP__;
<FORM ACTION="clone.pl" METHOD=POST>
<b class=orange>Clone settings to hostnode</b>
<br>
<INPUT TYPE=HIDDEN NAME=nexturl VALUE="/index.pl">
__STOP__
        $selhost="<SELECT NAME=CLONENODE SIZE=1>";
	open(HOST, $hostfile);
	while(<HOST>) {
	    chop;
	    $selhost.="<OPTION VALUE=\"$_\">$_" if ($_ !~ /^$/);
	}
	$selhost.="</SELECT>";
    }
    $OUT.=<<__STOP__
$selhost
<INPUT TYPE=SUBMIT VALUE="clone settings">
</FORM>
__STOP__
}

$OUT .= "<p class=\"info\"><A HREF=\"/index.pl?SID=$SID\" TARGET=\"_top\">Click for new settings to take effect immediateley</a></p>";

if ($nexturl) {
    $OUT .=" <A HREF=$nexturl>continue ....</A>";
}


PrintPage("Context",$OUT);


