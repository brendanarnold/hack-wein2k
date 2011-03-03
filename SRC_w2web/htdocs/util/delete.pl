#!/usr/bin/perl

require "../../libs/w2web.pl";
$doit=0;
$d=0;
$what = "file";

&GetInput;
&GetSession;
$what = "directory" if $d;

if (!$doit) {
	if ($DIR =~ /$file/){
		$OUT .= <<__STOP__;
<H1>ERROR</h1>
<p>
You can not delete this directory!
</p>
<p>
It is being used in the current session. To remove it first
select a different working directory for this session.
</p>
__STOP__
	} else {
		if ($d) {
		
			opendir(DIR, $file);
			@files = readdir(DIR);
			closedir(DIR);
			if ($#files>1) {
				$OUT .= <<__STOP__;
<p><b>WARNING:<br>
The selected directory contains the following files!</B><br>
These files will also be erased!</B></p>
<pre>
__STOP__
			foreach $i (@files) {$OUT .= "     $i\n";}
			}
			$OUT .= "</pre>";
		}

	$OUT .= <<__STOP__;
<br>
<h3>Delete $what</h3>

<p>Really delete <b>$file</b> ?
<br>
<FORM ACTION="/util/delete.pl" METHOD=POST>
<INPUT TYPE=HIDDEN NAME="file" VALUE="$file">
<INPUT TYPE=HIDDEN NAME="doit" VALUE="1">
<INPUT TYPE=HIDDEN NAME="d" VALUE="$d">
<INPUT TYPE=SUBMIT VALUE="Yes, delete it">
</FORM>
</p>
__STOP__

	}
} else {

	if ($d) {	
		$umps = qx( rm -rf $file );
	} else {
		$umps = qx( rm $file );
	}

	$OUT = <<__STOP__;
<h3>Delete $what</h3>

<b>$file</b> was deleted.
<br>
$umps
__STOP__
}


PrintPage("Delete $what", $OUT);
