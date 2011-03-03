#!/usr/bin/perl

require "../../libs/w2web.pl";

&GetInput;
&GetSession;

$next="continue with DOS";
$nexturl="/exec/dos.pl?SID=$SID";
$nextinteractive=1;

$OUT .= "<h2>Density of states</h2>";

&RequiredFile("int");
&InitTask();


$OUT.="<TABLE>";

if ($plot) {
} else {
	$OUT .=  <<__STOP__;
<TR><TD class="task">
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="lapw2">
$myspin
<INPUT NAME=qtl  VALUE="on" TYPE=HIDDEN>
<INPUT TYPE=SUBMIT VALUE="x lapw2 -qtl $mycomplex $myspinopt $mypara">
Calculate partial charges
$ni
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>



<TR><TD class="task">
<FORM ACTION="/util/edit.pl">
__STOP__
&PassHiddenParms;

$OUT.=<<__STOP__;
<INPUT NAME=file VALUE="$DIR/$CASE.int" TYPE=HIDDEN>
<INPUT TYPE=HIDDEN NAME="spin" VALUE="$spin">
<INPUT TYPE=SUBMIT VALUE="edit $CASE.int">
Edit input-file for TETRA
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>

<TR><TD class="task">
<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<INPUT NAME=prog TYPE=HIDDEN VALUE="tetra">
$myspin
<INPUT TYPE=SUBMIT VALUE="x tetra $myspinopt">
Calculate partial DOS
$ni
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>

<TR><TD class="task">
<FORM ACTION="/util/edit.pl">
<INPUT NAME=SID VALUE=$SID TYPE=HIDDEN>
<INPUT NAME=file VALUE="$DIR/$CASE.outputt$spin" TYPE=HIDDEN>
<INPUT NAME=redir VALUE="/exec/dos.pl?SID=$SID" TYPE=HIDDEN>
<INPUT TYPE=SUBMIT VALUE="edit $CASE.outputt$spin">
Check output of TETRA
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
</FORM>
</TD></TR>
__STOP__

}

if($plot){
    if($column){}else{$column=1};

	$OUT .= <<__STOP__;
<TR><TD class="task">
<b>We are in Dosplot mode, </b>
__STOP__

$myform1 =  <<__STOP__;
<TR><TD class="task">
<FORM ACTION=/exec/dos.pl METHOD=POST>
__STOP__

$numcol=qx(head -3 $DIR/$CASE.int | tail -1|cut -f1-6 -d' ');

#$myform2 =  "<SELECT NAME=column>";

#for ($i = 1; $i <= $numcol; $i++) {
#  $myform2 .= "<option value=$i>Col $i\n";
#}

#$myform2 .= <<__STOP__;
#</SELECT>

$myform2 = <<__STOP__;
Select column 1-$numcol:<INPUT NAME=column VALUE="$column" SIZE=2>

<INPUT TYPE=hidden NAME=doit VALUE="1">
<INPUT TYPE=hidden NAME=plot VALUE="1">
$myspin
<INPUT TYPE=SUBMIT VALUE="dosplot">
Plot DOS
<SELECT NAME="units">
<OPTION VALUE="ev">eV
<OPTION VALUE="">Ry
</SELECT>
<br>
Set ranges (optional):
<br>
xmin=<INPUT NAME=xmin VALUE="$xmin" SIZE=5>
xmax=<INPUT NAME=xmax VALUE="$xmax" SIZE=5>
ymin=<INPUT NAME=ymin VALUE="$ymin" SIZE=5>
ymax=<INPUT NAME=ymax VALUE="$ymax" SIZE=5>
</FORM>
</TD></TR>
__STOP__

	if($doit) {
		$fileext=1;	
		$fileext = 2 if ($column > 7 ) ;
		$fileext = 3 if ($column > 14) ;
		$fileext = 4 if ($column > 21) ;
		$fileext = 5 if ($column > 28) ;
		$fileext = 6 if ($column > 35) ;
		$fileext = 7 if ($column > 42) ;
		$fileext = 8 if ($column > 49) ;
		$fileext = 9 if ($column > 56) ;


		$grepline = $column+3;
    $line = qx(head -$grepline <$DIR/$CASE.int|tail -1);
		($atom,$qtl)=split(" ", $line);

		$OUT .= "Column $column selected:<br>";

		if($qtl >=0) {
		$grepline = $atom + 4 ;
    $line = qx(head -$grepline <$DIR/$CASE.qtl|tail -1);
		$line =~ s/^.* tot/tot/;
		@line = split (/,/, $line);
		$count = 0;
		foreach $i ( @line) {
			$count++;
			if ($count == $qtl) {
				$label=$i;
				$dum=$label;
				if ($i =~ /0/) {
					$label = "s-DOS";
				} elsif ($i =~ /1/) {
					$label = "p-DOS";
				} elsif ($i == "2") {
					$label = "d-DOS";
				} elsif ($i =~ /3/) {
					$label = "f-DOS";
				} 
			}
		}
		}
		if ($atom == 0) {
			$atom="total DOS";
			$qtl=-1;
                        $label="";
		} else {
			$atom="Atom $atom ";
		}
#		$label="";
		$OUT .= "at=$atom, qtl=$qtl" if ($debug);

		$OUT .= "Plot from atom: $atom, $label (file $CASE.dos$fileext$units$updn)";

	        $col =  ($column-1) % 7;
		$col++;
		$col++;
		$tmp1="$DIR/:dos1";
		$tmp2="$DIR/:dos2";
		$xlabel="Ry";
		$axis="noyzeroaxis";
		if ($units=~/ev/) {
			$xlabel="eV";
			$axis="yzeroaxis";
		}

		$ylabel="$atom $label";
		$dosfile="$DIR/$CASE.dos$fileext$units$updn";
		$plotfile="$tempdir/$SID-$$";
		$umps=qx(sed "1,3d" $dosfile >$tmp1);
		unless(open(FILE,">$tmp2")) {
			&CGIError("Can't write file $fname.\n");
			exit;
		}
		print FILE <<__STOP__;
show all
set terminal png
set output '$plotfile.png'
set title '$dosfile'
set data style lines
set xlabel "$xlabel"
set ylabel "$ylabel"
set xrange [$xmin:$xmax]
set yrange [$ymin:$ymax]
set $axis
plot '$tmp1' using 1:$col title '$label'
set terminal postscript
set output '$plotfile.ps'
replot
__STOP__
close(FILE);

	$umps = qx(cd $DIR;gnuplot $tmp2 2>&1);
	$OUT .= "umps=$umps" if $debug;;
	$OUT .= "<br><IMG SRC=/tmp/$SID-$$.png><br clear=all>";
	$OUT .= "<A HREF=/tmp/$SID-$$.ps>Download hardcopy in PostScript format</A
>";


		$OUT .= $myform1;
		&PassHiddenParms;
		$OUT .= $myform2
	} else {
		$OUT .= $myform1;
		&PassHiddenParms;
		$OUT .= $myform2

	}
	$OUT .=  <<__STOP__;
<TR><TD class="task">
<A HREF="$nexturl">Show full menu</A>
</TD></TR>
__STOP__
} else {
	$OUT .= <<__STOP__;
<TR><TD class="task">
<FORM ACTION=/exec/dos.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT TYPE=hidden NAME=plot VALUE="1">
$myspin
<INPUT TYPE=SUBMIT VALUE="dosplot">
Plot DOS
</FORM>
</TD></TR>
__STOP__
}
$OUT .= <<__STOP__;
</TABLE>
__STOP__


PrintPage("Context",$OUT);


