use Cwd qw(abs_path);
$w2webdir = "$ENV{'W2WEB'}";
$w2webhost = "$ENV{'SERVER_NAME'}";

#directories
$sessionpath = "$w2webdir/sessions";
$outpage = "$ENV{'DOCUMENT_ROOT'}/template.html";
$css = "/w2web.css";
$tempdir = "$w2webdir/tmp";
$WIENROOT = "$ENV{'WIENROOT'}";
$DIR = "";
#
#$debug = 0;
#
$refresh = 30;      # in seconds
#$refresh = 3;      # in seconds
$removetime=5;      # in minutes 
$comment = "";
$indent = "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;";

$supportnodes=1;		# support for remote nodes

# colors
$red="#ffc9c9";
$darkred="#F06262";
$green="#c3ffd2";
$gray="#E0E0E0";
$gray2="#C0C0C0";
$topcolor="#efe883";


# define link names/icons

$icons=0;
$jsmenu=0;
$elnesexpert=0;

$delete = "rm&nbsp;";
$delete = "<img src=/art/trash.gif BORDER=0 ALT=\"delete\">";
$rename = "mv&nbsp;";
$rename = "<img src=/art/move.gif BORDER=0 ALT=\"rename\">";
$copy = "cp&nbsp;";
$copy = "<img src=/art/copy.gif BORDER=0 ALT=\"copy\">";
$kill = "kill&nbsp;";
$kill = "<img src=/art/bomb.gif BORDER=0 ALT=\"kill\">";
$download="dl&nbsp;";
$download = "<img src=/art/dl.gif BORDER=0 ALT=\"download\">";
$upload ="ul&nbsp;";
$upload = "<img src=/art/ul.gif BORDER=0 ALT=\"upload\"> upload a file";
$new = "new&nbsp;";
$new = "<img src=/art/new.gif BORDER=0 ALT=\"new\"> create new file";
$reload = "reload";
$reload = "<img src=/art/reload.gif BORDER=0 ALT=\"reload\"> reload";


@savelist = qw/NAME SID HOSTNODE DIR ALERT TIME spinpol afm complex p COMMENT/;
#@savelist = qw/NAME SID HOSTNODE DIR ALERT TIME COMMENT/;

$inactivecolor=$red;
$activecolor=$green;
$spin="up"; 				# default spin is up
$updn=""; 				# default spin is up
$myspinopt="";
$filec="";


# Program Options:

sub GetSession {
	$session = $SID;
	my $lines;
	my $line;
	my $fname;

	#get user config 
	$fname = "$w2webdir/conf/user.conf";
		if (-e $fname) {
				unless(open(FILE,"$fname")) {
			&CGIError("Can't read file $fname.\n");
			exit;
		}
		@lines = <FILE>;
		close(FILE);
		chomp(@lines);
		for $line (@lines) {
			my($key, $content) = split(/=/, $line, 2);
			$$key = $content;
		}
	}
	#get session info 
	$fname = "$sessionpath/$session";
	unless(open(FILE,"$fname")) {
		&CGIError("Can't read file $fname.\n");
		exit;
	}
	@lines = <FILE>;
	close(FILE);
	chomp(@lines);
	for $line (@lines) {
		my($key, $content) = split(/=/, $line, 2);
		$$key = $content;
	}                              
	$CASE = $DIR;
	$CASE =~ s/.*\///m;
	if ($spinpol=~/on|CHECKED/) {
	} else {
		$spin="";
	}

}


$pname = "$w2webdir/conf/platinum.conf";
sub platinumchange {
    $outpage = "$ENV{'DOCUMENT_ROOT'}/platinum.html";
    $css = "/platinum.css";
}

sub PlatinumRead {
    if (-e $pname) {
	unless(open(FILE,"$pname")) {
	    &CGIError("Can't read file $pname.\n");
	    exit;
	}
	@lines = <FILE>;
	close(FILE);
	chomp(@lines);
	for $line (@lines) {
	    my($key, $content) = split(/=/, $line, 2);
	    $$key = $content;
	}
	&platinumchange if $platinum;
    }
}
sub PlatinumWrite {
    if ($platinum) {
	$platinum=0 if ($platinum ==-1);
	unless(open(FILE,">$pname")) {
	    &CGIError("Can't write file $pname.\n");
	    exit;
	}
	$platinumcount++;
	print FILE "platinum=$platinum\n";
	print FILE "platinumcount=$platinumcount\n";
	close(FILE);
    }
}


&PlatinumRead;


sub RequiredFile {
	my $testfile = shift(@_);
	$file="$DIR/$CASE.$testfile";
	if (-e $file && -s $file ) {
		#good there is an input file
	} else {
		#$OUT .= "$CASE.$testfile generated from template!<br>\n ";
		$umps=qx(cp $WIENROOT/SRC_templates/case.$testfile $file);
		if ($editit) {
			redirectURL("/util/edit.pl?SID=$SID&f=1&file=$DIR/$CASE.$testfile");
			exit;
		}
	}
}

sub SaveSession {
	$session = $SID;
	my $lines;
	my $line;
	my $fname;

	$fname = "$sessionpath/$session";
	unless(open(FILE,">$fname")) {
		&CGIError("Can't write file $fname.\n");
		exit;
	}
	$OUT .= "Save Session<br>" if $debug;

	if ($spinpol =~ /on/) {
		$spinpol = "CHECKED"
	}
	if ($complex =~ /on/) {
		$complex = "CHECKED"
	}
	if ($p =~ /on/) {
		$p = "CHECKED"
	}
	if ($afm =~ /on/) {
		$afm = "CHECKED"
	}

	for $mist (@savelist) {
		$OUT .= "$mist = $$mist<br>" if $debug;
		if ($mist =~ /COMMENT/ ) {
			$$mist =~ s/\n/<br>/g;
		}
		print FILE "$mist=$$mist\n";
	}
	print FILE "TIME=";
	$TIME = localtime;
	print FILE "$TIME\n";
	close(FILE);
		
}

sub PassHiddenParms {
	$exc = shift (@_);            
	for $umps (@savelist) {
		if ($umps eq $exc) {
			$OUT .= "$umps: aetsch<br>" if $debug;
		} else {
			 $OUT .= "<INPUT NAME=\"$umps\" VALUE=\"$$umps\" TYPE=HIDDEN>\n" unless ($except eq $umps);
		}
	}     

}
sub ShowParms {
	for $umps (@savelist) {
	 $OUT .= "$umps -> $$umps<br>\n";
	}     
}

sub GetInput {
	if ($ENV{'REQUEST_METHOD'} ne 'POST') {
		$buffer = $ENV{'QUERY_STRING'};
	} else {
		read(STDIN, $buffer, $ENV{'CONTENT_LENGTH'});
	}
	my @pairs = split(/&/, $buffer);
	my $pair;
	my $myindex;
	my $mykey;
	foreach $pair (@pairs) {
		my($key, $content) = split(/=/, $pair, 2);
		$content =~ tr/+/ /;
		$content =~ s/%([a-fA-F0-9][a-fA-F0-9])/pack("C", hex($1))/eg;
		$key =~ s/%([a-fA-F0-9][a-fA-F0-9])/pack("C", hex($1))/eg;
		$FORM{$key} = $content;
		if ($key =~ m/\[.*\]\[.*\]/) {
			$mykey = $key;
# old			$mykey =~ s/\[.*\]//eg;
                        $mykey =~ s/\[(.*)\]\[(.*)\]//;
# for new perl5.8.0  suse8.1
			$myi1 = $key;
			$myi2 = $key;
			$myi1 =~ s/.*\[(.*)\]\[(.*)\]/$1/;
			$myi2 =~ s/.*\[(.*)\]\[(.*)\]/$2/;
			$$mykey[$myi1][$myi2]=$content;
			$OUT .= "<br>$mykey [ $myi1 ] [ $myi2 ] = $content" if $debug;
		} elsif ($key =~ m/.*\[.*\]/) {
			$myindex = $key;
			$mykey = $key;
			$myindex =~ s/.*\[//;
			$myindex =~ s/\]//;
			$mykey =~ s/\[.*\]//;
			$$mykey[$myindex]=$content;
			$OUT .= "<br>$mykey [ $myindex ] = $content" if $debug;
		} else {
			$OUT .= "<br>$key -> $content <br>" if $debug;
			$$key = $content;
		}
	}
	if ($FORM{'SID'}) {
		$SID = $FORM{'SID'};
	}
	if ($FORM{'dir'}) {
		# preventing troubles coming from automount homes
		# thx to Florent Boucher
		$home = abs_path($ENV{"HOME"});
		$dir = abs_path($FORM{'dir'});
		$dir =~ s|$home|$ENV{"HOME"}|;
		$FORM{'dir'} = $dir;
	}
	foreach $dum (@savelist) {
		if ($FORM{$dum})  {
			$$dum = $FORM{$dum};
		}
	}
}


sub CGIError {
	print "Content-type: text/html\n\n";
	print @_;
	exit;
}

sub PrintPage {
	($title, $body) = @_;
	my(@page);
	open(PAGE, $outpage) || die("CANNOT OPEN $outpage. Reason: $!");
	@page = <PAGE>;
	close(PAGE);
	my($page);
	foreach (@page) {
		$page .= $_;
	}
	$page =~ s/\*\*TITLE\*\*/$title/eg;
	$page =~ s/\*\*COMMENT\*\*/$comment/eg;
	$page =~ s/\*\*BODY\*\*/$body/eg;
	print "Content-type: text/html\n\n";
	print $page;
	exit;
}

sub redirectURL {
	$durl = (shift @_);
	if ($splash) {
	print "Content-type: text/HTML\n
<HTML>
<HEAD>
   <META HTTP-EQUIV=\"refresh\"
               content=\"2; URL=$durl\">
</HEAD>
<BODY>
<br><br><br><br>
<H1 ALIGN=CENTER>$splash</H1>
</BODY>
</HTML>";
	} else {
	print "Content-type: text/HTML\n
<HTML>
<HEAD>
   <META HTTP-EQUIV=\"refresh\"
               content=\"0; URL=$durl\">
</HEAD>
<BODY>
</BODY>
</HTML>";
	}

}
sub listdir {
        my $dirpath = shift(@_);
        opendir(DIR, "$dirpath");
        @raw = sort grep(!/^\./, readdir(DIR));
        closedir(DIR);
        @ascii_files = ();
        @subdirs = ();
        @binaries = ();
        for $item(@raw) {
                if(-d "$dirpath/$item") {
                        push(@subdirs, $item);
                }elsif(-T "$dirpath/$item") {
                        push(@ascii_files, $item);
                }else {
                        push(@binaries, $item);
                }
        }
}                                          


sub InitTask {
    $updn="";
    if($spinpol=~ /CHECKED/ ) {
	$OUT .= "<A HREF=$nexturl&spin=up>Spin UP</A> <A HREF=$nexturl&spin=dn>Spin DOWN</A><br><br>";
	$OUT .= "<b>Spin ";
	$OUT .= "UP " if ($spin =~ /up/) ;
	$OUT .= "DOWN " if ($spin =~ /dn/) ;
	$updn = "up" if ($spin =~ /up/) ;
	$updn = "dn" if ($spin =~ /dn/) ;
		$OUT .= "selected.</b><br>";
	$myspin="<INPUT TYPE=HIDDEN NAME=spinpol VALUE=on><INPUT TYPE=HIDDEN NAME=spin VALUE=$spin>";
	$myspinopt="-$spin";
	$nexturl.="&spin=$spin";
    }
    
    $filec="";
    if($complex =~ /CHECKED/ ) {
	$mycomplex="-c";
	$filec="c";
    }
    $mypara="";
    if($p =~ /on|CHECKED/ ) {
	$mypara="-p";
    }
    $ni="<INPUT TYPE=checkbox NAME=nextinteractive CHECKED> interactively";
    $ni="<INPUT TYPE=checkbox NAME=nextinteractive> interactively" if $taskback;
}


sub ExeTypes {
	my $myfile;
  $exetypes="Type of execution: <SELECT NAME=exetype><OPTION VALUE=\"background\">background<OPTION VALUE=\"interactive\">interactively";
  $exetypes2="Type of execution: <SELECT NAME=exetype><OPTION VALUE=\"interactive\">interactively<OPTION VALUE=\"background\">background";
	$myfile = "$w2webdir/conf/execution.conf";
	if (-e $myfile) {
		open(EXE, $myfile);
		while(<EXE>) {
			chop;
			if (/^#/ || !/\S/) { next; }
			/^([^=]+)=(.*)$/;
			$name = $1; $val = $2;
			$ck="";
			$ck="SELECTED" if ($PREFS{'exetype'} eq $val);
			$exetypes.="<OPTION VALUE=\"$val\" $ck>$name\n";
			$exetypes2.="<OPTION VALUE=\"$val\" $ck>$name\n";
		}
		close(EXE);
	}
	
	$exetypes .= "</SELECT>";
}

sub SavePrefs {
	if ($prefspace eq "single" || $prefspace eq "scf") {
		@prefs = qw(h afm1 it it_count it0 renorm in1orig p1 so dm orb fsm fsm_count nohns nohns_count in1new in1new_count ql ql_count itnum itnum_count expert conv_ec conv_fc conv_cc convval_ec convval_fc convval_cc notify exetype spin qtl band sigma);
	} elsif ($prefspace eq "ana") {
		push @prefs, @anapar1;
		push @prefs, @anapar2;
		push @prefs, @anapar3;
		push @prefs, @anapar4;
	} else {
		@prefs = qw( );

	}
  $fname = "$DIR/.prefs-$prefspace";
  unless(open(PREF,">$fname")) {
    &CGIError("Can't write file $fname.\n");
    exit;
  }
  $OUT .= "Save preferences<br>" if $debug;
	foreach $i (@prefs) {
	$OUT.="$i=$$i<br>\n" if $debug;
		print PREF "$i=$$i\n";
	}
	if ($prefspace eq "scf") {
##		print PREF "conv_=CHECKED" if (!$conv);
#		print PREF "conv_ec=CHECKED" if ($conv eq "ec");
#		print PREF "conv_fc=CHECKED" if ($conv eq "fc");
#		print PREF "conv_cc=CHECKED" if ($conv eq "cc");
	}
	if ($prefspace eq "single") {
		print PREF "prog_$prog=CHECKED";
	}
	if ($prefspace eq "ana") {
		print PREF "c_scffile=on" if ($scffile eq "default");
		print PREF "c_altfile=on" if ($scffile eq "alt");
	}
	close(PREF);
}
sub GetPrefs {
  $fname = "$DIR/.prefs-$prefspace";
	#set default values if not set...
	$PREFS{'itnum_count'}=20;
	$PREFS{'ql_count'}=0.05;
	$PREFS{'in1new_count'}=2;
	$PREFS{'nohns_count'}=6;
	$PREFS{'fsm_count'}=0.;
	$PREFS{'it_count'}=6;
	$PREFS{'altfile'}="$CASE.scf";
	$PREFS{'lines'}=10;
	$PREFS{'altlines'}=1;
	$PREFS{'c_scffile'}="CHECKED";
	$PREFS{'spin'}="up";
	$PREFS{'conv_ec'}="CHECKED";
	$PREFS{'convval_ec'}=0.0001;
	$PREFS{'convval_fc'}=1.0;
	$PREFS{'convval_cc'}=0.001;
	#
	if (-e $fname) {
		unless(open(FILE,"$fname")) {
			&CGIError("Can't read file $fname.\n");
			exit;
		}
		@lines = <FILE>;
		close(FILE);
		chomp(@lines);
		for $line (@lines) {
			my($key, $content) = split(/=/, $line, 2);
			$content="CHECKED" if ($content eq "on");
			$PREFS{$key} = $content;
			$OUT.=" $key -> $content<br>" if $debug;
		}
	}
	$PREFS{'c_scffile'}="" if $PREFS{'c_altfile'};
	$upc="CHECKED" if ($PREFS{'spin'} eq "up");
	$dnc="CHECKED" if ($PREFS{'spin'} eq "dn");
}

1;


