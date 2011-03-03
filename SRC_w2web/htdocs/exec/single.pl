#!/usr/bin/perl

require "../../libs/w2web.pl";
&GetInput;
&GetSession;
$prefspace="single";
&GetPrefs;
&ExeTypes();

$OUT .=  <<__STOP__;
<H2>Single Programs</H2>

<FORM ACTION=/exec/executor.pl METHOD=POST>
<INPUT TYPE=HIDDEN NAME="saveprefs" VALUE="1">
__STOP__

&PassHiddenParms;

$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="x">
<table border=0 cellspacing=0 cellpadding=1>

<tr valign=top><td colspan=5><H3>Programs for
initialisation</H3></td></tr>
<tr valign=top><td>$indent
</td><td><INPUT NAME=prog TYPE=radio VALUE=nn $PREFS{'prog_nn'}>&nbsp;nn
</td><td><INPUT NAME=prog TYPE=radio VALUE=sgroup $PREFS{'prog_sgroup'}>&nbsp;sgroup
</td><td><INPUT NAME=prog TYPE=radio VALUE=symmetry $PREFS{'prog_symmetry'}>&nbsp;symmetry
</td><td><INPUT NAME=prog TYPE=radio VALUE=lstart $PREFS{'prog_lstart'}>&nbsp;lstart
</td></tr>
<tr valign=top><td>
</td><td><INPUT NAME=prog TYPE=radio VALUE=kgen $PREFS{'prog_kgen'}>&nbsp;kgen
</td><td><INPUT NAME=prog TYPE=radio VALUE=dstart $PREFS{'prog_dstart'}>&nbsp;dstart
</td><td><INPUT NAME=prog TYPE=radio VALUE=afminput $PREFS{'prog_afminput'}>&nbsp;afminput
</td></tr>

<tr valign=top><td colspan=5><H3>SCF Programs</H3></td></tr>
<tr valign=top><td>$indent</TD>
</td><td><INPUT NAME=prog TYPE=radio VALUE=lapw0 $PREFS{'prog_lapw0'}>&nbsp;lapw0
</td><td><INPUT NAME=prog TYPE=radio VALUE=lapw1 $PREFS{'prog_lapw1'}>&nbsp;lapw1
</td><td><INPUT NAME=prog TYPE=radio VALUE=lapwso $PREFS{'prog_lapwso'}>&nbsp;lapwso
</td><td><INPUT NAME=prog TYPE=radio VALUE=lapw2 $PREFS{'prog_lapw2'}>&nbsp;lapw2</td></tr>
<tr valign=top><td>
</td><td><INPUT NAME=prog TYPE=radio VALUE=lapwdm $PREFS{'prog_lapwdm'}>&nbsp;lapwdm
</td><td><INPUT NAME=prog TYPE=radio VALUE=orb $PREFS{'prog_orb'}>&nbsp;orb
</td><td><INPUT NAME=prog TYPE=radio VALUE=lcore $PREFS{'prog_lcore'}>&nbsp;lcore
</td><td><INPUT NAME=prog TYPE=radio VALUE=mixer $PREFS{'prog_mixer'}>&nbsp;mixer
</td></tr>


<tr valign=top><td colspan=5><H3>Other programs</h3></td></tr>
<tr valign=top><td>$indent</TD>
</td><td><INPUT NAME=prog TYPE=radio VALUE=aim $PREFS{'prog_aim'}>&nbsp;aim
</td><td><INPUT NAME=prog TYPE=radio VALUE=clmcopy $PREFS{'prog_clmcopy'}>&nbsp;clmcopy
</td><td><INPUT NAME=prog TYPE=radio VALUE=clminter $PREFS{'prog_clminter'}>&nbsp;clminter
</td><td><INPUT NAME=prog TYPE=radio VALUE=elnes $PREFS{'prog_elnes'}>&nbsp;elnes
</td></tr>
<tr valign=top><td>
</td><td><INPUT NAME=prog TYPE=radio VALUE=filtvec $PREFS{'prog_filtvec'}>&nbsp;filtvec
</td><td><INPUT NAME=prog TYPE=radio VALUE=irrep $PREFS{'prog_irrep'}>&nbsp;irrep
</td><td><INPUT NAME=prog TYPE=radio VALUE=joint $PREFS{'prog_joint'}>&nbsp;joint
</td><td><INPUT NAME=prog TYPE=radio VALUE=lapw3 $PREFS{'prog_lapw3'}>&nbsp;lapw3
</td></tr>
<tr valign=top><td>
</td><td><INPUT NAME=prog TYPE=radio VALUE=lapw5 $PREFS{'prog_lapw5'}>&nbsp;lapw5
</td><td><INPUT NAME=prog TYPE=radio VALUE=lapw7 $PREFS{'prog_lapw7'}>&nbsp;lapw7
</td><td><INPUT NAME=prog TYPE=radio VALUE=optic $PREFS{'prog_optic'}>&nbsp;optic
</td><td><INPUT NAME=prog TYPE=radio VALUE=spaghetti $PREFS{'prog_spaghetti'}>&nbsp;spaghetti
</td></tr>
<tr valign=top><td>
</td><td><INPUT NAME=prog TYPE=radio VALUE=supercell $PREFS{'prog_supercell'}>&nbsp;supercell
</td><td><INPUT NAME=prog TYPE=radio VALUE=telnes $PREFS{'prog_telnes'}>&nbsp;telnes
</td><td><INPUT NAME=prog TYPE=radio VALUE=tetra $PREFS{'prog_tetra'}>&nbsp;tetra
</td><td><INPUT NAME=prog TYPE=radio VALUE=txspec $PREFS{'prog_txspec'}>&nbsp;txspec
</td></tr>
<tr valign=top><td>
</td><td><INPUT NAME=prog TYPE=radio VALUE=xspec $PREFS{'prog_xspec'}>&nbsp;xspec
</td></tr>

<tr valign=top><td colspan=5><HR>
<H3>Options: (<INPUT NAME=h TYPE=CHECKBOX >help)</H3> 
</td></tr>
<tr valign=top><td>
</td><td><INPUT NAME=c TYPE=CHECKBOX $complex>&nbsp;complex
</td><td><INPUT NAME=it TYPE=CHECKBOX $PREFS{'it'}>&nbsp;iterative
</td><td><INPUT NAME=so TYPE=CHECKBOX $PREFS{'so'}>&nbsp;spinorbit
</td></tr>
<tr valign=top><td>
</td><td><INPUT NAME=spinpol TYPE=CHECKBOX $spinpol>&nbsp;spin polarized: 
</td><td><INPUT NAME=spin TYPE=radio VALUE=up $upc>&nbsp;spin up
</td><td><INPUT NAME=spin TYPE=radio VALUE=dn $dnc>&nbsp;spin down
</td><td><INPUT NAME=orb TYPE=CHECKBOX $PREFS{'orb'}>&nbsp;orb
</td></tr>

<tr valign=top><td>
</td><td><INPUT NAME=p TYPE=CHECKBOX $p>&nbsp;parallel
</td><td><INPUT NAME=nohns TYPE=CHECKBOX $PREFS{'nohns'}>&nbsp;nohns
</td><td><INPUT NAME=qtl TYPE=CHECKBOX $PREFS{'qtl'}>&nbsp;qtl
</td><td><INPUT NAME=band TYPE=CHECKBOX $PREFS{'band'}>&nbsp;band
</td></tr>
<tr>
<td>
</td><td><INPUT NAME=d TYPE=CHECKBOX $PREFS{'d'}>&nbsp;only def-file
</td></tr>

<tr valign=top><td colspan=5>
<br>
$exetypes
<INPUT TYPE=SUBMIT VALUE="Execute!">
</td></tr>
</table>
</FORM>
__STOP__

PrintPage("Context",$OUT);


