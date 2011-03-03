#!/usr/bin/perl

require "../../libs/w2web.pl";

@needinput = qw(nn lstart kgen optimize afminput supercell lapw3);
@options = qw(c it p spinpol so d qtl h dm orb nohns renorm in1orig band sigma); 

$debug=0;

$interactive=0;
$inputok=1;
$proginput="";
$again="";
$precmd="";
$prog="";

&GetInput();
&GetSession();
&ExeTypes();

if ($saveprefs) {
	$prefspace="single";
	&SavePrefs();
}

if($FORM{'p'}) {
	# oh dear, parallel execution!
	if ( -e "$DIR/.machines" ) {
		# good, machines file found.
		$OUT .= ".machines file found" if $debug;
	} else {
		#$OUT .= "NO .machines file found" if $debug;
		redirectURL("/exec/para.pl?SID=$SID");
		exit 0;
	}
}


if ("$FORM{'prog'}" =~  m/^x.*/) {
	# if we came from "Execute command line, we have to move
	# x to the "precommand" string.

	$FORM{'precmd'} = 'x';
	$FORM{'prog'} =~ s/x //;
}
if ("$FORM{'prog'}" =~  m/^xterm/) {
	$FORM{'precmd'} = '';
#	$ENV{'DISPLAY'}="$ENV{'REMOTE_HOST'}:0.0";
  $FORM{'prog'} = "xterm";
	$interactive=0;
}
if ("$FORM{'prog'}" =~  m/^run.*/) {
    $FORM{'precmd'} = '';
# probably we want an scf? 
    &redirectURL("/exec/scf.pl?SID=$SID");
    exit;
}


#--------------------------------------
sub WriteParms {
	$OUT .= "<FORM ACTION=/exec/executor.pl>";
	$OUT .= <<__STOP__;
<INPUT TYPE=HIDDEN NAME=next VALUE="$next">
<INPUT TYPE=HIDDEN NAME=nexturl VALUE="$nexturl">
<INPUT TYPE=HIDDEN NAME=nextinteractive VALUE=$nextinteractive>           
__STOP__
	foreach $i ("prog", "precmd", @options) {
		$OUT .= qq/<INPUT TYPE=hidden NAME="$i" VALUE=$FORM{$i}>\n/
	}
	&PassHiddenParms();
}
#--------------------------------------
sub nn {
	$comment = "NN needs input";
	&WriteParms();
	# nn needs only one input:
	$OUT .= <<__STOP__;
please specify nn-bondlength factor: (usually=2)<br>
<INPUT NAME="INPUT1" VALUE=2>
__STOP__
}

sub nn_post {
	$umps = qx/cd $DIR;grep ERROR *.outputnn|wc/;
	if ( $umps > 0 ) {
		$OUT .<<__STOP__;
<p><b>An error occured in <i>nn</i>!</B></p>
<p>
You should check RMT radii!
</p>
__STOP__
	}
}
#--------------------------------------
sub lapw3 {
	$comment = "LAPW3 needs input";
	&WriteParms();
	$OUT .= <<__STOP__;
<p>
<b>Specify sin (theta/lambda) until which structure factors should be calculated:  </b>
<INPUT NAME="INPUT1" VALUE="1.0"> <br>
</p>
<p>
<b>Using w2web you can only run lapw3 only with total density!!</b><br>
<INPUT NAME="INPUT2" TYPE=HIDDEN VALUE="TOT">
</p>
__STOP__
}
#--------------------------------------
sub lstart {
	$comment = "LSTART needs input";
	&WriteParms();
	$OUT .= <<__STOP__;
<p>
<b>Select Exchange Correlation Potential:</b><br>
<SELECT NAME="INPUT1">
<OPTION VALUE="13">GGA (Perdew-Burke-Ernzerhof 96)
<OPTION VALUE="5">LSDA
<OPTION VALUE="11">GGA (Wu-Cohen 06)
</SELECT>
</p>
<p>
<b>ENERGY to separate core and valence states:</b><br>
<INPUT NAME="INPUT2" VALUE="-6.0"> (recommended: -6.0 Ry)<br>
(check how much core charge leaks out of MT-sphere)
</p>
__STOP__
}
#--------------------------------------
sub kgen {
	&WriteParms();

	if ( -e "$DIR/$CASE.klist") {
		$mist = qx(cp $DIR/${CASE}.klist $DIR/${CASE}.klist_x);
		$mist = qx(cp $DIR/${CASE}.kgen $DIR/${CASE}.kgen_x);
	}

			if ("$FORM{so}" =~ /on|CHECKED/) {
	$umps = qx/cd $DIR;echo "0\n0\n0\n0\n0"|x kgen -so 2>&1/;
    } else {
	$umps = qx/cd $DIR;echo "0\n0\n0\n0\n0"|x kgen 2>&1/;
    }
	if ( -e "$DIR/$CASE.klist_x") {
		$mist = qx(mv $DIR/$CASE.klist_x $DIR/$CASE.klist);
		$mist = qx(mv $DIR/$CASE.kgen_x $DIR/$CASE.kgen);
	}

	$kpts = qx/cd $DIR;head -1 *.klist |cut -c47-56|sed "s\/ *\/\/"/;
	if ("$umps" =~  m/.*add inversion.*/) {
		$OUT .= "WE DON'T HAVE INVERSION<br>" if $debug;
		$OUT .= <<__STOP__;
Add inversion?
<SELECT NAME="INPUT1">
<OPTION VALUE=0>No
<OPTION VALUE=1>Yes
</SELECT>
__STOP__
	} else {
		$OUT .= "WE HAVE INVERSION" if $debug;
		$OUT .= "<INPUT NAME=INPUT1 VALUE=\" \" TYPE=HIDDEN>"
	}
 
	$OUT .= <<__STOP__;
<p>
<b>Number of k-points:</b>
<INPUT NAME="INPUT2" SIZE=10 VALUE="$kpts">
<INPUT NAME="INPUT3" TYPE=HIDDEN VALUE=0>
<br>
Shift k-mesh (if applicable)
<SELECT NAME="INPUT4">
<OPTION VALUE=1>Yes
<OPTION VALUE=0>No
</SELECT>
<br>
                          
</P>
__STOP__

}
#--------------------------------------
sub supercell {
	&WriteParms();
		$OUT .= <<__STOP__;
<b>Filename of original struct file:</b>
<INPUT NAME="INPUT1" size=20 value="$CASE.struct">
<p> <p><b>Number of cells in x direction:</b>
<INPUT NAME="INPUT2" size=2 value="1">
<br><b>Number of cells in y direction:</b>
<INPUT NAME="INPUT3" size=2 value="1">
<br><b>Number of cells in z direction:</b>
<INPUT NAME="INPUT4" size=2 value="1">
<p> <p><b>Optional shift of all atoms (fractional coordinates) </b> 
<br>in x direction:
<INPUT NAME="INPUT5" size=5 value="0.0">
<br>in y direction:
<INPUT NAME="INPUT6" size=5 value="0.0">
<br>in z direction:
<INPUT NAME="INPUT7" size=5 value="0.0">
__STOP__
#			$OUT .= qq/<INPUT TYPE=SUBMIT VALUE="Execute!">/;

#	$umps = qx/cd $DIR;echo "tic.struct\n1\n1\n1\nP\n0\n0\n0\n0"|x supercell 2>&1/;
#	$umps = qx/cd $DIR;echo "$INPUT1\n$INPUT2\n$INPUT3\n$INPUT4\nP\n0\n0\n0\n0"|x supercell 2>&1/;

#	@output = split(/\n/,$umps);
#	foreach $i (@output) {
#		$OUT .= "$i<br>";
#	if ("$i" =~  m/.*Enter your target lattice type.*/) {
#		$OUT .= "$i<br>";
		$OUT .= <<__STOP__;
<p><b>Enter your target lattice type:</b>
<SELECT NAME="INPUT8">
<OPTION VALUE=P>P
<OPTION VALUE=F>F
<OPTION VALUE=B>B
</SELECT> <br>(Some choices are restricted by symmetry)

__STOP__
#	}
#	    }
		$OUT .= <<__STOP__;
<p><b>For surfaces or isolated molecules:</b> (for P lattice only)
<br>Vacuum in x-direction [bohr]:
<INPUT NAME="INPUT9" size=15 value="0"> 
<br>Repeat atoms at x=0:
<SELECT NAME="INPUT10">
<OPTION VALUE=N>N
<OPTION VALUE=Y>Y
</SELECT>
<p>Vacuum in y-direction [bohr]:
<INPUT NAME="INPUT11" size=15 value="0"> 
<br>Repeat atoms at y=0:
<SELECT NAME="INPUT12">
<OPTION VALUE=N>N
<OPTION VALUE=Y>Y
</SELECT>
<p>Vacuum in z-direction [bohr]:
<INPUT NAME="INPUT13" size=15 value="0"> 
<br>Repeat atoms at z=0:
<SELECT NAME="INPUT14">
<OPTION VALUE=N>N
<OPTION VALUE=Y>Y
</SELECT>

<p>
__STOP__
#	$kpts = qx/cd $DIR;head -1 *.klist |cut -c47-56|sed "s\/ *\/\/"/;
#	if ("$umps" =~  m/.*add inversion.*/) {
#		$OUT .= "WE DON'T HAVE INVERSION<br>" if $debug;
#		$OUT .= <<__STOP__;
#Add inversion?
#<SELECT NAME="INPUT1">
#<OPTION VALUE=0>No
#<OPTION VALUE=1>Yes
#</SELECT>
#__STOP__
#	} else {
#		$OUT .= "WE HAVE INVERSION" if $debug;
#		$OUT .= "<INPUT NAME=INPUT1 VALUE=\" \" TYPE=HIDDEN>"
#	}
# 
#	$OUT .= <<__STOP__;
#<p>
#<b>Number of k-points:</b>
#<INPUT NAME="INPUT2" SIZE=10 VALUE="$kpts">
#<INPUT NAME="INPUT3" TYPE=HIDDEN VALUE=0>
#<br>
#Shift k-mesh (if applicable)
#<SELECT NAME="INPUT4">
#<OPTION VALUE=1>Yes
#<OPTION VALUE=0>No
#</SELECT>
#<br>
#                          
#</P>
#__STOP__
#
}
#--------------------------------------
sub optimize {
	$comment = "optimize needs input";
	&WriteParms();
	$OUT .= <<__STOP__;
<h2>optimizer</h2>
<SELECT NAME=INPUT1>
<OPTION VALUE="1">1) vary VOLUME with constant a:b:c ratio
<OPTION VALUE="2">2) vary C/A-RATIO with constant volume (tetr.,hex.lattic)
<OPTION VALUE="3">3) vary C/A-RATIO with constant volume and b/a (orthor.lattic)
<OPTION VALUE="4">4) vary B/A-RATIO with constant volume and c/a (orthor.lattic)
<OPTION VALUE="5">5) vary A and C (2D-case) (tetragonal or hexagonal lattice)
<OPTION VALUE="6">6) vary A, B and C (3D-case) (orthorhombic lattice)
<OPTION VALUE="7">7) vary A, B, C and Gamma (4D-case) (monoclinic lattice)
</SELECT>
<br>
<b>For option 1-4:</B> specify structure changes in % (each value in separate line)<br>
<TEXTAREA NAME=INPUT3 ROWS=5>
</TEXTAREA>
<br>

<b>For option 5:</B> specify number of structures: 6, 9 (3x3), 16 (4x4), 25 (5x5), 36 
 <br>

<b>For option 6:</B> specify number of structures: 10, 27 (3x3x3), 64 (4x4x4), 125 (5x5x5) 
<INPUT NAME="INPUT2_5" size=4 value="16"> specify the % change:<INPUT NAME="INPUT3_5" size=6 value="1.0"> <br>

<b>For option 7:</B> specify number of structures: 15, 81 (3x3x3x3), 256 (4x4x4x4) 
 <br>


__STOP__
}


#--------------------------------------
sub afminput {
	$comment = "afminput needs input";
	&WriteParms();
	# create dummy input with 200 no answers ...
	unless(open(FILE,">$DIR/.dum")) {
		&CGIError("Can't write file $fname.\n");
		exit;
	}
	for ($i = 1; $i < 200; $i++) {
		print FILE "n\n";
	}
	close(FILE); 
	$out =  qx( cd $DIR ; cat .dum | x afminput );
#	$out =  qx( cd $DIR ; cat .dum | x afminput | grep KLASSENGLEICH);
	@output = split(/\n/,$out);

	if ("$out" =~  m/.*Found a symmetry*/) {
	$OUT .= "<p> case.struct_supergroup present";
	$OUT .= "<p> Symmetry operation found! Click on Execute";
	$count =0;
	 } elsif ("$out" =~  m/.*group present*/ ) {
	$OUT .= "<p> case.struct_supergroup present";
	$OUT .= "<p> The super and subgroups are KLASSENGLEICH";
	$OUT .= "<p> You must specify a TRANSLATION VECTOR which transforms the";
	$OUT .= "<br>spin-up into the spin-dn atom (e.g.: 0.5,0.5,0.5 for AFM bcc Cr)!";
;
	$OUT .= <<__STOP__;
        <p> <INPUT NAME="afm1" SIZE=10 VALUE="0.5">
           <INPUT NAME="afm2" SIZE=10 VALUE="0.5">
           <INPUT NAME="afm3" SIZE=10 VALUE="0.5"> 
__STOP__
	$count =3;
	 } elsif ("$out" =~  m/.*NOT present*/ ) {
	$OUT .= "<p> case.struct_supergroup NOT present";
        $OUT .= "<p> It is strongly recommended that you copy the (nonmagnetic) supergroup"; 
        $OUT .= "<br> struct file to case.struct_supergroup (unless they are KLASSENGLEICH)";
        $OUT .= "<br> Otherwise:";
        $OUT .= "<p> You must specify a symmetry operation (rotation + translation vector)"; 
        $OUT .= "<br> which transforms the spin-up into the spin-dn atom (e.g. for AFM bcc Cr:)!";
        $OUT .= "<p>  1 0 0  0.5";
        $OUT .= "<br>  0 1 0  0.5";
        $OUT .= "<br>  0 0 1  0.5";

	$OUT .= <<__STOP__;
        <p> <INPUT NAME="afm1" SIZE=6 VALUE="1.0">
           <INPUT NAME="afm2" SIZE=6 VALUE="0.0">
           <INPUT NAME="afm3" SIZE=6 VALUE="0.0"> 
           <INPUT NAME="afm4" SIZE=10 VALUE="0.0"> 
        <p> <INPUT NAME="afm5" SIZE=6 VALUE="0.0">
           <INPUT NAME="afm6" SIZE=6 VALUE="1.0">
           <INPUT NAME="afm7" SIZE=6 VALUE="0.0"> 
           <INPUT NAME="afm8" SIZE=10 VALUE="0.0"> 
        <p> <INPUT NAME="afm9" SIZE=6 VALUE="0.0">
           <INPUT NAME="afm10" SIZE=6 VALUE="0.0">
           <INPUT NAME="afm11" SIZE=6 VALUE="1.0"> 
           <INPUT NAME="afm12" SIZE=10 VALUE="0.0"> 
__STOP__
	$count =12;
}
#	$count =0;
#	foreach $i (@output) {
#		$count++;
#		$i =~ s/NO MATCH:/<SELECT NAME=afm$count><OPTION VALUE=n>no<OPTION VALUE=y>yes<\/SELECT>/;
#		$OUT .= "$i<br>";
#	}
	$OUT.= <<__STOP__;
<INPUT TYPE=HIDDEN NAME=INPUT1 VALUE=$count>
__STOP__
	

}



if (!($FORM{'prog'})) {
	$OUT .= <<__STOP__;
<P>
<B>No program selected!</B>
</P>
<A HREF=$ENV{'HTTP_REFERER'}>return to previous page</A>

__STOP__
	PrintPage("Context",$OUT);
	exit;
}


if ($FORM{'INPUT1'} eq "") {
	$need = 0;
	foreach $i (@needinput) {
#		if ( $i =~ /$FORM{'prog'}/){;
#		if ( $FORM{'prog'} =~ /$i/){;
		if ( $FORM{'prog'} eq "$i" ){;
			$need = 1;
			$OUT .= "NEED INPUT!!!" if $debug;
			$inputok=0;
			&$i;
			$OUT .= qq/<INPUT TYPE=SUBMIT VALUE="Execute!">/;
		}
	}
} else {
	if ($FORM{'prog'}=~/optimize/) {
		#in optimizer we need to know how many changes there are ...
		if ( $FORM{'INPUT1'} <= 4 ){
		$INPUT2=qx(echo "$INPUT3"|wc -w);
		$FORM{'INPUT2'}=$INPUT2;
		$FORM{'INPUT3'}=~ s/\n\n/\n/; 
	    } else {
		$FORM{'INPUT2'}=$INPUT2_5;
		$FORM{'INPUT3'}=$INPUT3_5; 
            }
	}
	if ($FORM{'prog'}=~/afminput/) {
		$tmp = "";
		for ($i=1; $i <= $FORM{'INPUT1'}; $i++) {
			$var="afm$i";
			$tmp .= "$FORM{$var}\n";
			
		}
		$FORM{'INPUT1'}=$tmp;
	}
	foreach $i (qw/1 2 3 4 5 6 7 8 9 10 11 12 13 14/) {
		my $var;
		$interactive = 1;
		$var = "INPUT$i";
		if ("$INPUT$i" != " ") {
			$proginput .= "$FORM{$var}\n";
		}
	}
}

if ($inputok) {
	# create commandline
	$cmd = "";
	$opts = "";
	foreach $i (@options) {
		if ($i eq "spinpol" || $i eq "spin") {
			if ($i eq "spinpol") {
				if ("$spinpol" =~ /on|CHECKED/ && $spin) {
				 $opts .= "-$spin ";
				}
			}
		} else {
			if ("$FORM{$i}" =~ /on|CHECKED/) {
				$opts .= "-$i ";
			}
		}
	}

	if ($FORM{'precmd'}) {
		$opts .= " -c" if ($complex =~ /on|CHECKED/);
		$cmd ="$FORM{'precmd'} $FORM{'prog'} $opts";
	} else {
		$cmd = "$FORM{'prog'}";
	}
	$again = <<__STOP__;
<H4>Execute another command line:</H4>
<p>
<table><tr valign=bottom><td>
<FORM ACTION=/exec/executor.pl METHOD=POST>
<INPUT NAME=prog TYPE=INPUT SIZE=50>
<BR>
<INPUT TYPE=HIDDEN NAME=SID VALUE=$SID>
<INPUT TYPE=SUBMIT VALUE="Run!">
$exetypes2
</SELECT>
</FORM>
</td><td>
<FORM ACTION=/exec/executor.pl METHOD=POST>
<INPUT NAME="prog" VALUE="$cmd" TYPE=HIDDEN>
<INPUT TYPE=HIDDEN NAME=SID VALUE=$SID>
<INPUT NAME="interactive" VALUE="$interactive" TYPE=HIDDEN>
<INPUT TYPE=SUBMIT VALUE="Repeat: $cmd">
</FORM>
</td></tr></table>
</p>
__STOP__


# if $next is set, just overwrite the $again....
if ($next) {
	$interactive=$nextinteractive;
	$again = <<__STOP__;
<H4>Continue with</H4>
<table><tr valign=bottom><td>
<FORM ACTION=/exec/next.pl METHOD=POST>
<INPUT TYPE=HIDDEN NAME=next VALUE=$next>
<INPUT TYPE=HIDDEN NAME=nexturl VALUE=$nexturl>
<INPUT TYPE=HIDDEN NAME=interactive VALUE=$nextinteractive>
<INPUT TYPE=SUBMIT VALUE="$next">
</FORM>
</td></tr></table>
__STOP__
}

	$OUT .= <<__STOP__;
<P>
Commandline: <b>$cmd</b><BR>
Program input is: <b>"$proginput"</b>
<br>
<br>
<table border=1><tr><td>
<PRE>
__STOP__

if ($exetype =~ /interactive/ ) {
} elsif ($exetype eq "" ) {
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
echo '$proginput'|$cmd >$DIR/STDOUT 2>&1
#$sendmail
rm $cmdfile
__STOP__
close(FILE);
$umps=qx(chmod +x $cmdfile);
}
	if ($exetype =~ /background/ ) {$exetype =""}
	if ($exetype =~ /interactive/ ) {$interactive=1;}

	system "rm $DIR/STDOUT" if ( -e "$DIR/STDOUT");                     
	if ($interactive) {
		if ($proginput) {
			$OUT .= "proginput \n" if $debug;
			$umps = qx( cd $DIR;echo "$proginput"|$cmd  2>&1);
		}else {
			$OUT .= "no input \n" if $debug;
			$umps = qx( cd $DIR;$cmd 2>&1);
		}
		$OUT .= "$umps\n";
	} else {
	if ($exetype =~ /%f/) {
		$exetype=~ s/%f/ $cmdfile/;
		$OUT.="we did: cd $DIR;$exetype &";
		system "cd $DIR;$exetype &";
	} else {
		$OUT.= "cd $DIR;$exetype echo '$proginput'|$cmd >$DIR/STDOUT 2>&1 &";
		system "cd $DIR;$exetype echo '$proginput'|$cmd >$DIR/STDOUT 2>&1 &";
		system "cd $DIR;rm $cmdfile";
	}
}
#		system "cd $DIR;echo '$proginput'|$cmd >$DIR/STDOUT 2>&1 &";
#		$OUT .= "commandline sent to system for execution";
#
#	}

	$OUT .= <<__STOP__;
</PRE>
</td></tr></table>
</p>
__STOP__

if ($interactive){
} else {
$OUT .= <<__STOP__;
	<A HREF="/util/stdout.pl?SID=$SID">View STDOUT</A>
__STOP__
}
$OUT .= <<__STOP__;
$again
__STOP__

}

PrintPage("Context",$OUT);


