#!/usr/bin/perl

use CGI; 
require "../../libs/w2web.pl";

#&GetInput;
#&GetSession;
$onnum = 1;
my $req = new CGI; 
$SID = $req->param("SID"); 
$cif = $req->param("cif"); 
&GetSession;
my $file = $req->param("FILE$onnum"); 


$OUT .=<<__STOP__;
	<h3>Upload done</h3>
<p>
 File $file was uploaded into $DIR
</p>
__STOP__


if ($file ne "") {
my $fileName = $file; 
$fileName =~ s!^.*(\\|\/)!!; 
$newmain = $fileName;

open (OUTFILE, ">$DIR/$fileName"); 
#$OUT .= "$DIR/$fileName<br>";
while (my $bytesread = read($file, my $buffer, 1024)) { 
	print OUTFILE $buffer; 
} 
close (OUTFILE); 
}

#$OUT = "Content-type: text/html\n";
#$OUT .= "Location:$donepage\n\n";
if ($cif) {
&redirectURL("/util/structask.pl?SID=$SID");
exit;
}
PrintPage("upload", $OUT);
