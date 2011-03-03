#!/usr/bin/perl

# This is the path to the directory that you want to edit.
# Note that all subdirectories in this directory will also be accessible.

require "../../libs/w2web.pl";
#$debug=1;
$mycols=1;
&GetInput;            
&GetSession;


#if($FORM{'redir'}) {
#	$redir = "<INPUT TYPE=hidden NAME=redir VALUE=\"$FORM{'redir'}\">";
#}	

if($FORM{'path'}) {  
	$path = $FORM{'path'};
	$ENV{'SCRIPT_NAME'} = "$ENV{'SCRIPT_NAME'}?$path";
}

if ($FORM{'edit'}) {
	$FORM{'edit'} =~ s/<\;/</g;
	$FORM{'edit'} =~ s/>\;/>/g;
##blaha	$FORM{'edit'} =~ s/"\;/\"/g;
}

if ($FORM{'diradd'} eq '..') {
	$FORM{'dir'} =~ s/\/[^\/]*?$//;
} elsif ($FORM{'diradd'} ne '') {
	$FORM{'dir'} .= "/$FORM{'diradd'}";
}

$path = "$path$FORM{'dir'}/";

if ($FORM{'file'}) {
	unless(open(FILE,"$path$FORM{'file'}")) {
		&system_error("Can't read file $path$FORM{'file'}.\n");
		exit;
	}
	@lines = <FILE>;
	close(FILE);
	chomp(@lines);
	for $line (@lines) {
		$line =~ s/\t/	/g;
##blaha		$line =~ s/\"/"\;/g;
		$line =~ s/</<\;/g;
		$line =~ s/>/>\;/g;

	}

	$OUT .= <<__STOP__;
	<FONT SIZE=+2>File:</font>
  <br> 
  $FORM{'file'}
	<table bgcolor=$gray border=1>
	<tr>
__STOP__

	if($next) {
	$mycols=2;
	$OUT .= <<__STOP__;
	<td bgcolor=$green>
<FORM ACTION="/exec/next.pl" METHOD=POST>
<INPUT TYPE=HIDDEN NAME=nexturl VALUE="$nexturl">
<INPUT TYPE=submit VALUE="$FORM{'next'}">
</FORM>
	</td>
__STOP__
	} else {
	}

	$header="";
	$testfile="$DIR/$CASE";
	if ($file =~ /$testfile.int/) {
		$header = "Header from $CASE.qtl$spin:<br>\n<pre>";
		$header .= qx(grep JATOM $DIR/$CASE.qtl$spin);
		$header =~  s/MULT.* tot/tot/g;
		$header =~  s/JATOM/ATOM/g;
		$header .= "</pre>\n";
	}
	if ($file =~ /$testfile.insp/) {
		$header = "Header from $CASE.qtl$spin and possible FERMI energies:<br>\n<pre>";
		$header .= qx(grep JATOM $DIR/$CASE.qtl$spin);
		$header =~  s/MULT.* tot/tot/g;
		$header =~  s/JATOM/ATOM/g;
		$header .= qx(grepline_lapw :FER "$DIR/*.scf" 1);
		$header =~  s/in.* files://g;
		$header =~  s/:FER  ://g;
		$header =~  s/F E R M I - ENERGY/EF /g;
		$header .= "</pre>\n";
	}


	$OUT .= <<__STOP__;
	<td>
	<form action="/util/edit.pl" METHOD=POST>
	<input type=submit value="Save">
	&nbsp;&nbsp;Download this file: 
	<a href=\"/util/download.pl?SID=$SID&file=$FORM{'file'}" class="img"><img
src=/art/dl.gif BORDER=0 ALT="download"></a>
	</td>
	<tr><td colspan=$mycols>
	$header
	</DIV><DIV CLASS=PRE>
	<textarea rows=20 cols=72 name=edit wrap=off>
__STOP__
#blaha: <pre> war vor textarea	
#	for (@lines) {$OUT .= "$_\n"}
	for (@lines) {$OUT .= "$_\n"}
	$OUT .= <<__STOP__;
</textarea>
	</DIV><DIV CLASS=MAIN>
	<input type=hidden name=filename value="$FORM{'file'}">
	<input type=hidden name=dir value="$FORM{'dir'}">
__STOP__

	if($next) {
	$OUT .= <<__STOP__;
	<INPUT TYPE=HIDDEN NAME=nexturl VALUE="$nexturl">
	<INPUT TYPE=HIDDEN NAME=next VALUE="$next">
__STOP__
	}

	$OUT .= <<__STOP__;
	</form>
	</td></tr>
__STOP__
if (!($next)) {
	$dum=$mycols-1;
	$OUT .= <<__STOP__;
	<tr><td align=center  bgcolor=$red>
	<FORM ACTION="/util/delete.pl" METHOD=POST>
	<INPUT TYPE=HIDDEN NAME=SID VALUE=$SID>
	<INPUT TYPE=HIDDEN NAME=file VALUE="$FORM{'file'}">
	<INPUT TYPE=SUBMIT VALUE="delete this file">
  </td><td colspan=$dum  bgcolor=$gray2>&nbsp;</td>
</tr>
__STOP__
}

	$OUT .= <<__STOP__;
	</FORM>
	</table>
__STOP__

	PrintPage("Edit",$OUT);
	exit;

} elsif ($FORM{'edit'}) {
	unless(open(NEW,">$path$FORM{'filename'}")) {
			&system_error("Couldn't write to $path$FORM{'filename'}.\n");
			exit;
	}
	$FORM{'edit'} =~ s/\r//gs;
#	$FORM{'edit'} =~ s/	/\t/gs;
	print NEW "$FORM{'edit'}";
	close(NEW);
}


if ($FORM{'redir'}) {
	redirectURL($FORM{'redir'});
} else {
	$ShowParms;
	$OUT .= "<p><b>File saved.</b></p>";
	if ($next) {
		$OUT .= <<__STOP__;
<FORM ACTION="/exec/next.pl">
<INPUT TYPE=HIDDEN NAME=nexturl VALUE="$nexturl">
<INPUT TYPE=submit VALUE="$FORM{'next'}">
</FORM>
__STOP__
	}
PrintPage("Context",$OUT);
}

exit;
&listdir($path);

print "Content-type: text/html\n\n";
print "<html><head></head><body>\n";
print "<form action=$ENV{'SCRIPT_NAME'} METHOD=POST>\n";
print "<H2>Path:</H2>\n";
print "$path\n";
print "<hr>\n";
print "<H2>Directories:</H2>\n";
print "<input type=radio name=diradd value=\"..\">.. |\n" if $FORM{'dir'};
for $dir (@subdirs) {
	print "<input type=radio name=diradd value=\"$dir\">$dir |\n";
}
print "<hr>\n";
print "<H2>Files:</H2>\n";
for $file (@ascii_files) {
	print "<input type=radio name=file value=\"$file\">$file |\n";
}
print "<input type=hidden name=dir value=\"$FORM{'dir'}\">\n";
print "<hr><input type=submit></form>\n";
print "</body></html>";

exit;

################################
# Error Subs
sub system_error {
		local($errmsg) = @_;
		&print_header("System Error");
		print $errmsg;
		print "<p><i><a href=\"$ENV{HTTP_REFERER}\">back to last page</a></i></p>\n";
		&print_footer;
}

sub print_header {
		local($title) = @_;
		print "Content-type: text/html\n\n";
		print "<HTML>\n";
		print "<HEAD>\n";
		print "<TITLE>$title</TITLE>\n";
		print "</HEAD>\n";
		print "<BODY>\n";
		print "<H1>$title</H1>\n";
}

sub print_footer {
		print "</BODY>\n";
		print "</HTML>\n";
}
#################################################
# listdir
# Given:
#	The path to the directory
#
# Returns:
#	@subdirs - the list of subdirectories
#	@ascii_files - the list of ascii files
#	@binaries - the list of binary files
#################################################

#.sub listdir {
#	my $dirpath = shift(@_);
#	opendir(DIR, "$dirpath");
#	@raw = sort grep(!/^\./, readdir(DIR));
#	closedir(DIR);
#	@ascii_files = ();
#	@subdirs = ();
#	@binaries = ();
#	for $item(@raw) {
#		if(-d "$dirpath/$item") {
#			push(@subdirs, $item);
#		}elsif(-T "$dirpath/$item") {
#			push(@ascii_files, $item);
#		}else {
#			push(@binaries, $item);
#		}
#	}
#}
