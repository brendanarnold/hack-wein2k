#!/usr/bin/perl

# REQUIRE LIBRARIES #
require "../../libs/w2web.pl";

$debug=0;


&GetInput;
&GetSession;
my $OUT;
#####################
if (!($FORM{'dir'})) {
    &CGIError("No directory or file specified");
}
if (!(-e $FORM{'dir'})) {
    &CGIError("$FORM{'dir'} does not exist");
}


$dir = $FORM{'dir'};
$fileext ="";
if ($FORM{'ext'}) {
    $fileext= $FORM{'ext'};
    $OUT .= "fileext found: $fileext<br>" if $debug;
}

if (-d $FORM{'dir'}) {
    opendir(DIR, $dir);
    @files = sort(readdir(DIR));
    closedir(DIR);
    $ext = "SID=$SID";
    $OUT .= <<__STOP__;
<A HREF="$ENV{'SCRIPT_NAME'}?$ENV{'QUERY_STRING'}" class=\"img\">$reload</A><br>
__STOP__


     if ($FORM{'cd'} == 1) {
	 $OUT .= <<__STOP__;
<h2>Current directory: $FORM{'dir'}</h2>
<FORM ACTION="/session/save.cgi" METHOD=POST>
<p class="info">Please change into or create a new directory and "select" it.
This must not be your home directory!</b><br>
__STOP__
          $ext .= "&cd=1";
	 $OUT .= <<__STOP__;
<INPUT TYPE=HIDDEN NAME=SID VALUE=$SID>
<INPUT TYPE=HIDDEN NAME=dir VALUE=$FORM{'dir'}>
<INPUT TYPE=SUBMIT VALUE="Select current directory">
</FORM>
</p>
__STOP__
}
    if ($FORM{'f'} == 1) {
	$OUT .= "<h2>Files in: $FORM{'dir'}</h2>\n";
	$OUT .= "<a href=\"$ENV{'SCRIPT_NAME'}?$ext&dir=$FORM{'dir'}\">show directories only</a><br><br>\n";
	$ext .= "&f=1";
    } elsif ($FORM{'ext'}) {
	$OUT .= "<h2>Files: $FORM{'dir'}/*.$fileext*</h2>\n";
    } else {
	$OUT .= "<br><a href=\"$ENV{'SCRIPT_NAME'}?$ext&dir=$FORM{'dir'}&f=1\">show files</a> in $FORM{'dir'}<br><br>\n";
    }
    
    if ($FORM{'cd'}) {
	$OUT .= <<__STOP__;
<br>
<hr>
<FORM ACTION="/util/mkdir.pl" METHOD=POST>
<INPUT TYPE=HIDDEN NAME="SID" VALUE="$SID">
<INPUT TYPE=HIDDEN NAME="dir" VALUE="$dir">
<INPUT TYPE=HIDDEN NAME="f" VALUE="$f">
<INPUT TYPE=HIDDEN NAME="cd" VALUE="$cd">
New directory: 
<INPUT NAME="newdir">
<INPUT TYPE=SUBMIT VALUE="create">
</FORM>
<FORM ACTION="/util/dir.pl" METHOD=POST>
<INPUT TYPE=HIDDEN NAME="SID" VALUE="$SID">
<INPUT TYPE=HIDDEN NAME="f" VALUE="$f">
<INPUT TYPE=HIDDEN NAME="cd" VALUE="$cd">
Quick cd: 
<INPUT NAME="dir" VALUE="$dir/" SIZE=40>
<INPUT TYPE=SUBMIT VALUE="cd">
</FORM>
<hr>
__STOP__
}

    foreach (@files) {
	if (!($_ =~ /^\.$/ )) {
	    if ($dir =~ /\/$/) {
		$F = $dir . $_;
	    } else {
		$F = $dir . "/" . $_;
	    }
	    
	    if (!($FORM{'ext'})) {
		if (-d $dir . "/" .  $_) {
		    if ( $_ !~ /\.\./ ) {
			$OUT .= "<a href=\"/util/copy.pl?SID=$SID&file=$F\" class=\"img\">$copy</a>";
			$OUT .= "&nbsp;<a href=\"/util/rename.pl?SID=$SID&file=$F\" class=\"img\">$rename</a>";
			$OUT .= "&nbsp;<a href=\"/util/delete.pl?SID=$SID&file=$F&d=1\" class=\"img\">$delete</a>";
			$OUT .= "&nbsp;<a href=\"/util/kill.pl?SID=$SID&file=$F\"class=\"img\">$kill</a>" if ($F =~ /.running/) ;
			$OUT .= "$indent\n"; 
			$mtime = localtime ((stat($F))[9]);
		    }
		    $OUT .= "$mtime <b><a href=\"$ENV{'SCRIPT_NAME'}?$ext&dir=$F\">$_</a></b><br>\n";
		} elsif (-f $dir . "/" . $_) {
		    if ($FORM{'f'} == 1) {
			$OUT .= "<a href=\"/util/copy.pl?SID=$SID&file=$F\" class=\"img\">$copy</a>";
			$OUT .= "&nbsp;<a href=\"/util/rename.pl?SID=$SID&file=$F\"
class=\"img\">$rename</a>\n";
			$OUT .= "&nbsp;<a href=\"/util/delete.pl?SID=$SID&file=$F\" class=\"img\">$delete</a>\n"; 
			$OUT .= "&nbsp;<a href=\"/util/download.pl?SID=$SID&file=$F\" class=\"img\">$download</a>";
			$OUT .= "&nbsp;<a href=\"/util/kill.pl?SID=$SID&file=$F\" class=\"img\">$kill</a>" if ($F =~ /.running/) ;
			$size = (stat($F))[7];
			$mtime = localtime ((stat($F))[9]);
			$OUT .= "$indent $mtime <a href=\"/util/edit.pl?$ext&file=$F\">$_</a>$indent $size<br>\n"; 
		    }
		} else {
#		$OUT .= " --- $_ <br>\n";
		}
	    } else {
		if ($_ =~ /\.$fileext.*$/) {
		    $OUT .= "<a href=\"/util/copy.pl?SID=$SID&file=$F\" class=\"img\">$copy</a>";
		    $OUT .= "&nbsp;<a href=\"/util/rename.pl?SID=$SID&file=$F\" class=\"img\">$rename</a>\n"; 
		    $OUT .= "&nbsp;<a href=\"/util/delete.pl?SID=$SID&file=$F\" class=\"img\">$delete</a>\n"; 
		    $OUT .= "&nbsp;<a href=\"/util/download.pl?SID=$SID&file=$F\"class=\"img\">$download</a>";
		$OUT .= "&nbsp;<a href=\"/util/kill.pl?SID=$SID&file=$F\" class=\"img\">$kill</a>" if ($F =~ /.running/) ;
		    $size = (stat($F))[7];
		    $mtime = localtime ((stat($F))[9]);
		    $OUT .= "$indent $mtime <a href=\"/util/edit.pl?$ext&file=$F\">$_</a>$indent $size<br>\n"; 
		}
	    }
	}
    }
    if ($FORM{'cd'} == 1) { 
	#nop
    } else {
	$OUT .="<br><br><a href=\"/util/upload.pl?SID=$SID\">$upload</a>";
	$OUT .="<br><a href=\"/util/newfile.pl?SID=$SID\">$new</a><br>";
	&done;
	if (-B $dir) {
	    $OUT .= "<b>$dir</b> is a binary file or empty. Sorry, we cannot open binary files now.";
	}
    }
    if (!(-e $dir)) {
	$OUT .= "<b>$dir</b> does not exist. Sorry";
	&done;
    }

}
&done;
sub done {
	$context = "Directory parser";
	PrintPage("Directories and Files", $OUT);
}
