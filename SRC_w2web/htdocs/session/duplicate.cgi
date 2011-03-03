#!/usr/bin/perl

require "../../libs/w2web.pl";
&GetInput;

&listdir($sessionpath);
$OUT .= <<__STOP__;
<H2>Duplicate a stored session:</H2>
<form action=/session/duplicate2.cgi METHOD=POST>
<SELECT NAME="SID" SIZE=5>
__STOP__

@sessions = ();

for $file (@ascii_files) {
  $SID=$file;
  &GetSession;
  push(@sessions,"$NAME"."³"."$file");
}

@newfile = qx( /bin/ls -tr $sessionpath );
my @mysessions = ();
@mysessions = sort @sessions;



foreach $line (@mysessions) {
  my $sid,$name;
  ($name, $sid) = split (/³/, $line);
  if ("$name" =~ m/.*$selection.*/i) {
  $OUT .= "<OPTION VALUE=$sid";
  $OUT .= " SELECTED" if ($sid == $newfile[$#newfile]);
  $OUT .= ">$name\n";
  }
}
$OUT .= <<__STOP__;
</SELECT>
<BR>
<H2>New session name:</H2>
<INPUT NAME="NEWNAME">
<br>
<INPUT TYPE=SUBMIT VALUE="Duplicate">
</FORM>

__STOP__


PrintPage("Change SID", $OUT);
