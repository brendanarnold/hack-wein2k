#!/usr/bin/perl

require "../../libs/w2web.pl";

&GetInput;
&GetSession;

$OUT = <<__STOP__;
<h3>File upload</h3>

<p>
<form method="POST" action="upload2.pl" ENCTYPE="multipart/form-data">
Select file to upload into $DIR<br>
<br>
__STOP__
if ($cif) {
$OUT.=<<__STOP__;
<input type=hidden name="cif" VALUE="$cif">
__STOP__
}
$OUT .= <<__STOP__;
<input type=hidden name="SID" VALUE="$SID">
<input type=file name="FILE1">
<br>
<input type=submit value="Upload">
</form>

__STOP__
PrintPage("upload", $OUT);
