#!/usr/bin/perl

require "../../libs/w2web.pl";
require "../../libs/struct.pl";
$doit=0;
$debug=0;

@anapar1=qw(ENE FER DIS NEC01 NEC02 MMTOT );
@anapar2=qw(FOR FGL POS QTL EFG ETA CHA DTO CTO NTO);
@anapar3=qw(CUP CDN HFF MMI);
@anapar4=qw(other otherparm lines altlines scffile altfile selectatom selectatom1 selectatom2 selectatom3 selectatom4 selectatom5  SCFMON ATOMX ATOMY ATOMZ);
$prefspace="ana";

&GetInput;
&GetSession;


	&StructRead;
	&GetPrefs();

	$atoms="<SELECT NAME=selectatom>";
	for ($i=0;$i<=$s_ineq;$i++) {
		$mycheck="";
		$mycheck="SELECTED" if ($i == $PREFS{selectatom});
		$atoms.="<OPTION VALUE=$i $mycheck>$i";
	}
	$atoms.="</SELECT>";
	$atoms1="<SELECT NAME=selectatom1>";
	for ($i=0;$i<=$s_ineq;$i++) {
		$mycheck="";
		$mycheck="SELECTED" if ($i == $PREFS{selectatom1});
		$atoms1.="<OPTION VALUE=$i $mycheck>$i";
	}
	$atoms1.="</SELECT>";
	$atoms2="<SELECT NAME=selectatom2>";
	for ($i=0;$i<=$s_ineq;$i++) {
		$mycheck="";
		$mycheck="SELECTED" if ($i == $PREFS{selectatom2});
		$atoms2.="<OPTION VALUE=$i $mycheck>$i";
	}
	$atoms2.="</SELECT>";
	$atoms3="<SELECT NAME=selectatom3>";
	for ($i=0;$i<=$s_ineq;$i++) {
		$mycheck="";
		$mycheck="SELECTED" if ($i == $PREFS{selectatom3});
		$atoms3.="<OPTION VALUE=$i $mycheck>$i";
	}
	$atoms3.="</SELECT>";
	$atoms4="<SELECT NAME=selectatom4>";
	for ($i=0;$i<=$s_ineq;$i++) {
		$mycheck="";
		$mycheck="SELECTED" if ($i == $PREFS{selectatom4});
		$atoms4.="<OPTION VALUE=$i $mycheck>$i";
	}
	$atoms4.="</SELECT>";
	$atoms5="<SELECT NAME=selectatom5>";
	for ($i=0;$i<=$s_ineq;$i++) {
		$mycheck="";
		$mycheck="SELECTED" if ($i == $PREFS{selectatom5});
		$atoms5.="<OPTION VALUE=$i $mycheck>$i";
	}
	$atoms5.="</SELECT>";


	$OUT .= "<h2>Analysis</h2>";
	if (!$doit) {
		$OUT .= <<__STOP__;
<FORM ACTION="/util/analysis.pl" METHOD=POST>
__STOP__

		&PassHiddenParms();
#		&StructRead();

		$OUT .= <<__STOP__;
<table bgcolor=$green>
<TR><TD BGCOLOR=$gray>
<b>atom independent parameters:</b>
<br>
$indent
<INPUT TYPE=CHECKBOX NAME="ENE" $PREFS{'ENE'}>ENE
<INPUT TYPE=CHECKBOX NAME="FER" $PREFS{'FER'}>FER
<INPUT TYPE=CHECKBOX NAME="DIS" $PREFS{'DIS'}>DIS
<INPUT TYPE=CHECKBOX NAME="NEC01" $PREFS{'NEC01'}>NEC-new
<INPUT TYPE=CHECKBOX NAME="NEC02" $PREFS{'NEC02'}>NEC-old
<INPUT TYPE=CHECKBOX NAME="MMTOT" $PREFS{'MMTOT'}>MMTOT
</TD></TR>

<TR><TD BGCOLOR=$gray>
<b>atom dependent parameters:</b>
<br>
$indent
<INPUT TYPE=CHECKBOX NAME="QTL" $PREFS{'QTL'}>QTL
<INPUT TYPE=CHECKBOX NAME="EFG" $PREFS{'EFG'}>EFG
<INPUT TYPE=CHECKBOX NAME="ETA" $PREFS{'ETA'}>ETA
<INPUT TYPE=CHECKBOX NAME="CHA" $PREFS{'CHA'}>CHA
<INPUT TYPE=CHECKBOX NAME="DTO" $PREFS{'DTO'}>DTO
<INPUT TYPE=CHECKBOX NAME="CTO" $PREFS{'CTO'}>CTO
<INPUT TYPE=CHECKBOX NAME="NTO" $PREFS{'NTO'}>NTO
</TD></TR>
<TR><TD BGCOLOR=$gray>
<b>atom dependent vector parameters:</b>
<br>
$indent
<INPUT TYPE=CHECKBOX NAME="FOR" $PREFS{'FOR'}>FOR
<INPUT TYPE=CHECKBOX NAME="FGL" $PREFS{'FGL'}>FGL
<INPUT TYPE=CHECKBOX NAME="POS" $PREFS{'POS'}>POS
&nbsp;&nbsp;(<INPUT TYPE=CHECKBOX NAME="ATOMX" $PREFS{'ATOMX'}>x-
<INPUT TYPE=CHECKBOX NAME="ATOMY" $PREFS{'ATOMY'}>y-
<INPUT TYPE=CHECKBOX NAME="ATOMZ" $PREFS{'ATOMZ'}>z-coordinate for scfmonitor)
</TD></TR>

<TR><TD BGCOLOR=$gray>
<b>for spin polarized systems:</b>
<br>
$indent
<INPUT TYPE=CHECKBOX NAME="CUP" $PREFS{'CUP'}>CUP
<INPUT TYPE=CHECKBOX NAME="CDN" $PREFS{'CDN'}>CDN
<INPUT TYPE=CHECKBOX NAME="HFF" $PREFS{'HFF'}>HFF
<INPUT TYPE=CHECKBOX NAME="MMI" $PREFS{'MMI'}>MMI
</TD></TR>

<TR><TD BGCOLOR=$gray>
<b>other parameter:</b>
<br>
$indent
<INPUT NAME=other TYPE=CHECKBOX $PREFS{'other'}>
<INPUT NAME=otherparm VALUE=$PREFS{'otherparm'}>
</TD></TR>


<TR><TD BGCOLOR=$gray>
<b>Select atom</b> for atom dependent param. (0 means all atoms, up to 6 atoms possible)<br>
$indent $atoms &nbsp;&nbsp;$atoms1 &nbsp;$atoms2 &nbsp;$atoms3 &nbsp;$atoms4 &nbsp;$atoms5 
</TD></TR>


<TR><TD BGCOLOR=$gray>
<b>Analysis of:  
<INPUT TYPE=RADIO NAME=scffile VALUE="default" $PREFS{'c_scffile'}>
$CASE.scf </b> &nbsp;&nbsp;with <INPUT NAME=lines SIZE=3 VALUE=$PREFS{'lines'}> lines.
<br>
<b>or of alternate scf-files:</b> <INPUT TYPE=RADIO NAME=scffile VALUE="alt" $PREFS{'c_altfile'}>
<INPUT NAME=altfile VALUE="$PREFS{'altfile'}">&nbsp;&nbsp;with <INPUT NAME=altlines SIZE=3 VALUE=$PREFS{'altlines'}> lines.
</TD></TR>

<TR><TD BGCOLOR=$gray>
<INPUT NAME=doit VALUE=1 TYPE=hidden>
<INPUT TYPE=SUBMIT VALUE="Analyze scf file">
&nbsp;&nbsp;<INPUT TYPE=CHECKBOX NAME="SCFMON" $PREFS{'SCFMON'}>Graphics using scfmonitor (only for single scf file) 
</TD></TR>
</TABLE>
</FORM>
__STOP__



	} else {
		&SavePrefs();

		$testfile="$CASE.scf";
		if ($scffile=~ /alt/) {
			$testfile="$altfile";
		}

                  if($selectatom < 1) {
                        $testatom="";
		}else {
		    if($selectatom < 10) {
			$testatom="00$selectatom";
		}else {
                  if($selectatom < 100) {
                        $testatom="0$selectatom";
                  }else {
			$testatom="$selectatom";
		    }}}
                  if($selectatom1 < 1) {
                        $testatom1="";
		}else {
		    if($selectatom1 < 10) {
			$testatom1="00$selectatom1";
		}else {
                  if($selectatom1 < 100) {
                        $testatom1="0$selectatom1";
                  }else {
			$testatom1="$selectatom1";
		    }}}
                  if($selectatom2 < 1) {
                        $testatom2="";
		}else {
		    if($selectatom2 < 10) {
			$testatom2="00$selectatom2";
		}else {
                  if($selectatom2 < 100) {
                        $testatom2="0$selectatom2";
                  }else {
			$testatom2="$selectatom2";
		    }}}
                  if($selectatom3 < 1) {
                        $testatom3="";
		}else {
		    if($selectatom3 < 10) {
			$testatom3="00$selectatom3";
		}else {
                  if($selectatom3 < 100) {
                        $testatom3="0$selectatom3";
                  }else {
			$testatom3="$selectatom3";
		    }}}
                  if($selectatom4 < 1) {
                        $testatom4="";
		}else {
		    if($selectatom4 < 10) {
			$testatom4="00$selectatom4";
		}else {
                  if($selectatom4 < 100) {
                        $testatom4="0$selectatom4";
                  }else {
			$testatom4="$selectatom4";
		    }}}
                  if($selectatom5 < 1) {
                        $testatom5="";
		}else {
		    if($selectatom5 < 10) {
			$testatom5="00$selectatom5";
		}else {
                  if($selectatom5 < 100) {
                        $testatom5="0$selectatom5";
                  }else {
			$testatom5="$selectatom5";
		    }}}
		$output ="";
		$showpars="";
		#atom indep.
		foreach $i (@anapar1) {
			if ($$i) {
				$showpars.=":$i ";
				$output .= "--- $i -----------\n";
				if ($scffile=~ /alt/) {
				    $output .= qx(cd $DIR;grepline :$i '$testfile' $altlines);} 
				else {$output .= qx(cd $DIR; grep :$i $testfile | tail -$lines);}
				
			}
		}	
		#atom dep.
		foreach $i (@anapar2 ,  @anapar3) {
			if ($$i) {
				$output .= "--- atom dependend parameter $i -----------\n";
# all atoms
                  if($selectatom < 1) {
		      for ($i1=1;$i1<=$s_ineq;$i1++) {
			  if($i1 < 10) {
			      $i2="00$i1";
			  }else {
			      if($i1 < 100) {
				  $i2="0$i1";
			      }else {
				  $i2="$i1";
			      }}			  
#				$showpars.=":$i$i2 ";
			$lab = "none" ;
			if ( $ATOMX ) {$lab = "x"; $showpars.=":$i$i2$lab ";}
			if ( $ATOMY ) {$lab = "y"; $showpars.=":$i$i2$lab ";}
			if ( $ATOMZ ) {$lab = "z"; $showpars.=":$i$i2$lab ";}
			if ($lab =~ "none" ) {$showpars.=":$i$i2 ";}
			  if ($scffile=~ /alt/) {
			      $output .= qx(cd $DIR; grepline :$i$i2 '$testfile' $altlines);} 
			  else {$output .= qx(cd $DIR; grep :$i$i2 $testfile | tail -$lines);}
		      }
		  } else {
# single atoms
			$lab = "" ;
			if ( $ATOMX ) {$lab = "x";}
			if ( $ATOMY ) {$lab = "y";}
			if ( $ATOMZ ) {$lab = "z";}
			$showpars.=":$i$testatom$lab ";
			if($selectatom1 > 1){$showpars.=":$i$testatom1$lab ";}
			if($selectatom2 > 1){$showpars.=":$i$testatom2$lab ";}
			if($selectatom3 > 1){$showpars.=":$i$testatom3$lab ";}
			if($selectatom4 > 1){$showpars.=":$i$testatom4$lab ";}
			if($selectatom5 > 1){$showpars.=":$i$testatom5$lab ";}
		      if ($scffile=~ /alt/) {
			  $output .= qx(cd $DIR; grepline :$i$testatom '$testfile' $altlines); 
			  if($selectatom1 > 1) {$output .= qx(cd $DIR; grepline :$i$testatom1 '$testfile' $altlines);}
			  if($selectatom2 > 1) {$output .= qx(cd $DIR; grepline :$i$testatom2 '$testfile' $altlines);}
			  if($selectatom3 > 1) {$output .= qx(cd $DIR; grepline :$i$testatom3 '$testfile' $altlines);}
			  if($selectatom4 > 1) {$output .= qx(cd $DIR; grepline :$i$testatom4 '$testfile' $altlines);}
			  if($selectatom5 > 1) {$output .= qx(cd $DIR; grepline :$i$testatom5 '$testfile' $altlines);}
                      } 
		      else {
                          $output .= qx(cd $DIR; grep :$i$testatom $testfile | tail -$lines);
                          if($selectatom1 > 1) {$output .= qx(cd $DIR; grep :$i$testatom1 $testfile | tail -$lines);}
                          if($selectatom2 > 1) {$output .= qx(cd $DIR; grep :$i$testatom2 $testfile | tail -$lines);}
                          if($selectatom3 > 1) {$output .= qx(cd $DIR; grep :$i$testatom3 $testfile | tail -$lines);}
                          if($selectatom4 > 1) {$output .= qx(cd $DIR; grep :$i$testatom4 $testfile | tail -$lines);}
                          if($selectatom5 > 1) {$output .= qx(cd $DIR; grep :$i$testatom5 $testfile | tail -$lines);}
                      }}
			    }
		    }	
		#spinpol .
#		foreach $i (@anapar3) {
#			if ($$i) {
#				$showpars.="$i,";
#				$output .= "--- $i$testatom -----------\n";
#				if ($scffile=~ /alt/) {
#				    $output .= qx(cd $DIR; grepline :$i$testatom '$testfile' $altlines);} 
#				else {$output .= qx(cd $DIR; grep :$i$testatom $testfile | tail -$lines); }
#			}
#		}	
		#other
		if ($other) {
			$showpars.=":$otherparm ";
			$output .= "--- $otherparm -----------\n";
				if ($scffile=~ /alt/) {
				    $output .= qx( cd $DIR;grepline :$otherparm '$testfile' $altlines);} 
			else {$output .= qx(cd $DIR; grep :$otherparm $testfile | tail -$lines);}
		}

		$OUT.= <<__STOP__;
<table>
<tr>
<td >
<FORM ACTION=/util/analysis.pl METHOD=POST>

__STOP__

		&PassHiddenParms;
		foreach $i (@anapar1) {
			$OUT.= "<INPUT TYPE=HIDDEN NAME=$i VALUE=$$i>\n";
		}
		foreach $i (@anapar2) {
			$OUT.= "<INPUT TYPE=HIDDEN NAME=$i VALUE=$$i>\n";
		}
		foreach $i (@anapar3) {
			$OUT.= "<INPUT TYPE=HIDDEN NAME=$i VALUE=$$i>\n";
		}
		foreach $i (@anapar4) {
			$OUT.= "<INPUT TYPE=HIDDEN NAME=$i VALUE=$$i>\n";
		}
		$OUT.= <<__STOP__;
<INPUT TYPE=HIDDEN NAME=parameter VALUE=$parameter>
<INPUT TYPE=HIDDEN NAME=selectatom VALUE=$selectatom>
<INPUT TYPE=HIDDEN NAME=other VALUE=$other>
<INPUT TYPE=HIDDEN NAME=lines VALUE=$lines>
<INPUT TYPE=HIDDEN NAME=altlines VALUE=$altlines>
<INPUT TYPE=HIDDEN NAME=scffile VALUE=$scffile>
<INPUT TYPE=HIDDEN NAME=altfile VALUE=$altfile>
<INPUT TYPE=HIDDEN NAME=doit VALUE=1>
<INPUT TYPE=SUBMIT VALUE="reanalyze scf file">
</FORM>
</td>
<td >
<FORM ACTION=/util/analysis.pl METHOD=POST>
__STOP__
    &PassHiddenParms;
    $OUT.= <<__STOP__;
<INPUT TYPE=HIDDEN NAME=doit VALUE=0>
<INPUT TYPE=SUBMIT VALUE="select other parameter">
</FORM></td>
</tr>
</table>
__STOP__

		
if ($SCFMON) {  $OUT .= "SCFmonitor selected: &nbsp;&nbsp;(only up to 6 parameter)";
# 		$OUT .= "$umps = qx(cd $DIR;scfmonitor_lapw -p -i $lines -f $testfile $showpars);";
            if ($scffile=~ /alt/) {
 		$umps = qx(cd $DIR;scfmonitor_lapw -p -i $altlines -f $testfile $showpars); } 
            else {
		$umps = qx(cd $DIR;scfmonitor_lapw -p -i $lines -f $testfile $showpars); }
		$umps = qx(cp $DIR/scfmonitor.png $tempdir/$SID-$$.png);
#		$OUT .= "$umps und $DIR/scfmonitor.png $tempdir/$SID-$$.png";
		$OUT .= "<IMG SRC=/tmp/$SID-$$.png> <br clear=all>";
             }
		$OUT .= <<__STOP__;
<p>
Analysis of parameter:<br> $showpars<br>
in $testfile
(showing last $lines / $altlines lines)
</p>
<pre>
$output
</pre>

<table>
<tr>
<td >
<FORM ACTION=/util/analysis.pl METHOD=POST>

__STOP__

		&PassHiddenParms;
		foreach $i (@anapar1) {
			$OUT.= "<INPUT TYPE=HIDDEN NAME=$i VALUE=$$i>\n";
		}
		foreach $i (@anapar2) {
			$OUT.= "<INPUT TYPE=HIDDEN NAME=$i VALUE=$$i>\n";
		}
		foreach $i (@anapar3) {
			$OUT.= "<INPUT TYPE=HIDDEN NAME=$i VALUE=$$i>\n";
		}
		foreach $i (@anapar4) {
			$OUT.= "<INPUT TYPE=HIDDEN NAME=$i VALUE=$$i>\n";
		}
		$OUT.= <<__STOP__;
<INPUT TYPE=HIDDEN NAME=parameter VALUE=$parameter>
<INPUT TYPE=HIDDEN NAME=selectatom VALUE=$selectatom>
<INPUT TYPE=HIDDEN NAME=other VALUE=$other>
<INPUT TYPE=HIDDEN NAME=lines VALUE=$lines>
<INPUT TYPE=HIDDEN NAME=altlines VALUE=$altlines>
<INPUT TYPE=HIDDEN NAME=scffile VALUE=$scffile>
<INPUT TYPE=HIDDEN NAME=altfile VALUE=$altfile>
<INPUT TYPE=HIDDEN NAME=doit VALUE=1>
<INPUT TYPE=SUBMIT VALUE="reanalyze scf file">
</FORM>
</td>
<td >
<FORM ACTION=/util/analysis.pl METHOD=POST>
__STOP__
    &PassHiddenParms;
    $OUT.= <<__STOP__;
<INPUT TYPE=HIDDEN NAME=doit VALUE=0>
<INPUT TYPE=SUBMIT VALUE="select other parameter">
</FORM></td>
</tr>
</table>
__STOP__

	}

	PrintPage("Anaylsis", $OUT);

