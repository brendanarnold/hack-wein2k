#!/usr/bin/perl

require "../../libs/w2web.pl";
&GetInput;
&GetSession;

&listdir($sessionpath);
$OUT .= <<__STOP__;
<H2>Delete stored sessions:</H2>
<br><br>
Current session is:
<b>$NAME</b> (SID=$SID)<br>
<br>
<FORM ACTION="/session/delete2.cgi">
<INPUT NAME=SID VALUE=$SID TYPE=HIDDEN>
<INPUT TYPE=SUBMIT VALUE="Delete current session">
</FORM>
__STOP__

if ($SDEL) {
$OUT .= <<__STOP__;
or
<FORM ACTION="/util/dir.pl">
<INPUT NAME="SID" VALUE="$SID" TYPE="HIDDEN">
<INPUT NAME="dir" VALUE="$ENV{'HOME'}" TYPE="HIDDEN">
<INPUT NAME="cd" VALUE="1" TYPE="HIDDEN">
<INPUT TYPE=SUBMIT VALUE="Select other directory">
</FORM>
__STOP__
} else {
$OUT .= <<__STOP__;
<form action=/session/delete2.cgi METHOD=POST>
<SELECT NAME="SID" SIZE=5>
__STOP__

for $file (@ascii_files) {
	$SID=$file;
	&GetSession;
        $OUT .= "<OPTION VALUE=$file>$NAME\n";
}              
$OUT .= <<__STOP__;    
</SELECT>
<br>
<i>A session is a shortcut to your working directories.
Deleting a session does <b>not</b> delete any
calculation data on disk!</i>
<BR>
<INPUT TYPE=SUBMIT VALUE="Delete selected Session">
</FORM>


<br>
__STOP__

}

PrintPage("Delete SID", $OUT);
