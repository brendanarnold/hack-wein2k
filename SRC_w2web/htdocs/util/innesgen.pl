#!/usr/bin/perl

require "../../libs/w2web.pl";
require "../../libs/struct.pl";
require "../../libs/telnes2.pl";

#$debug=1;
&GetInput;
&GetSession;
&StructRead;
&InnesRead;


$OUT .=  <<__STOP__;
<H2>InnesGen<font size=-2><sup>TM</sup></font> for TELNES.2 </H2>
<table bgcolor=$gray>
<FORM ACTION="/util/innessave.pl" METHOD=POST>
<input type=hidden name="next" value="$next">
<input type=hidden name="nexturl" value="$nexturl";
<br><br>
__STOP__

&PassHiddenParms;

$OUT .=  <<__STOP__;

<tr><td bgcolor=$gray colspan=2>
<b>Title:&nbsp; </b> <INPUT NAME="t_title" VALUE="$t_title" SIZE=40>
</td></tr>

<tr><td>
<b>Atom:&nbsp;</b>
<select name="t_atom">
__STOP__

for ($i = 1; $i <= $s_nato; $i++) {
    $mysel="";
    $mysel="selected" if ($i == $t_atom);
    $OUT .= "<option value=\"$i\" $mysel>$i: $s_name[$i]$s_nameadd[$i]";
}


$OUT .=  <<__STOP__;
</select>
</td><td>
<b>Edge:</b>&nbsp;
<select name="t_edge">
<option value="">use n and l
<option value="K">K
<option value="L1">L1
<option value="L23">L23
</select>
(n=<input name="t_n" value="$t_n" size=2>
l=<input name="t_l" value="$t_l" size=2>)
</td></tr>
<tr><td>
<b>Edge onset:</b>&nbsp;
<input name="t_ene" value="$t_ene" size=8> eV
</td><td>
<b>Beam energy:</b>&nbsp;
<input name="t_beam" value="$t_beam" size=8> keV
</td></tr>
<tr><td colspan=2>
<b>Energy grid:</b>
<input name="t_egrid1" value="$t_egrid1" size=6>  eV to 
<input name="t_egrid2" value="$t_egrid2" size=6>  eV 
in steps of 
<input name="t_egrid3" value="$t_egrid3" size=6> eV
</td></tr>
<tr><td>
<b>Collection s.a.:</b>&nbsp;
<input name="t_coll" value="$t_coll" size=8> mrad
</td><td>
<b>Convergence s.a.:</b>&nbsp;
<input name="t_conv" value="$t_conv" size=8> mrad
</td></tr>

<tr><td>
<b>Spectrometer broadening</b>&nbsp;
<input name="t_spec" value="$t_spec" size=6> eV
</td><td>
<b>Q-mesh:</b>&nbsp;
NR=<input name="t_nr" value="$t_nr" size=3> 
NT=<input name="t_nt" value="$t_nt" size=3>
</td></tr>
__STOP__

if($elnesexpert) {


    $OUT .= <<__STOP__;
<tr><td bgcolor=$gray2 colspan=2>
<b>Advanced settings:</b>
</td></tr>



<tr><td bgcolor=$gray2 colspan=2>
<table><tr><td>
<b>Branching ratio:</b></td><td>
<input name="t_branch1" value="$t_branch1" size=6> (statistical if empty)
</td></tr>
<tr><td>
<b>Split:</b>
</td><td>
<input name="t_split1" value="$t_split1" size=6> (calculated if empty)
</td></tr></table>

</td></tr>

<tr><td bgcolor=$gray2 colspan=2>
__STOP__

$co=0;
$co = " checked" if $t_orient;

    $OUT .= <<__STOP__;
<input type=checkbox name="t_orient" $co>
<b>Orientation sensitive: </b>
&alpha;=<input name="t_orient1" value="$t_orient1" size=6>°,
&beta;=<input name="t_orient2" value="$t_orient2" size=6>°,
&gamma;=<input name="t_orient3" value="$t_orient3" size=6>°
</td></tr>

<tr><td colspan=2 bgcolor=$gray2>

<b>Integrate over eq. atoms:</b>
<input name="t_atoms1" value="$t_atoms1" size=3>  to
<input name="t_atoms2" value="$t_atoms2" size=3> (all eq. atoms if empty)
</td></tr>

<tr><td colspan=2 bgcolor=$gray2>
<b>Detector position:</b>
theta_x <input name="t_det1" value="$t_det1" size=6> mrad, 
theta_y <input name="t_det2" value="$t_det2" size=6> mrad
</td></tr>

<tr><td colspan=2 bgcolor=$gray2>

__STOP__

$esel="";
$esel=" selected" if ($t_modus1 =~ /energy/);
$asel="";
$asel=" selected" if ($t_modus1 =~ /angle/);

$OUT.=<<__STOP__;

<b>Modus:</b>
<select name="t_modus">
<option value="">----
<option value="energy" $esel>energy
<option value="angle" $asel>scattering angle
</select>

</td></tr>
<tr><td bgcolor=$gray2 colspan=2>

__STOP__
$mysel="";
$mysel=" selected" if ($t_xqtl == 1);

$OUT .= <<__STOP__;
<b>Use charges from: </b>
<select name="t_xqtl">
<option value="">x lapw2 -qtl (default)
<option value="xqtl" $mysel>x xqtl
</select>
</td></tr>
<tr><td bgcolor=$gray2 colspan=2>
__STOP__

my $c1 = "";
my $c2 = "";
my $c3 = "";
my $c4 = "";
$c1 = " checked" if ($t_init1 =~ /y/);
$c2 = " checked" if ($t_init2 =~ /y/);
$c3 = " checked" if ($t_init3 =~ /y/);
$c4 = " checked" if ($t_init4 =~ /y/);

$OUT .= <<__STOP__;
<table><tr valign=top><td rowspan=2>
<b>Initialization:</b>
</td><td>
<input type=checkbox name="t_init1" $c1> Calculate DOS
</td><td>
<input type=checkbox name="t_init2" $c2> write DOS
</td></tr>
<tr><td>
<input type=checkbox name="t_init3" $c3> Calculate rotation matrices
</td><td>
<input type=checkbox name="t_init4" $c4> write rotation matrices
</td></tr></table>

</td></tr>

<tr><td bgcolor=$gray2>

__STOP__



$mysel="";
$mysel=" selected" if ($t_noheaders == 1);
    $c0=$c1=$c2=0;
    $c0=" selected" if ($t_output1==0);
    $c1=" selected" if ($t_output1==1);
    $c2=" selected" if ($t_output1==2);

$OUT .= <<__STOP__;
<b>Verbosity:</b>
<select name="t_output">
<option value="">----
<option value="0" $c0>low
<option value="1" $c1>medium
<option value="2" $c2>crazy
</select>
</td><td bgcolor=$gray2>
<b>File headers: </b>
<select name="t_noheader">
<option value="">Write headers (default)
<option value="NOHEADERS" $mysel>no headers
</select>
</td></tr>

<tr><td bgcolor=$gray2 colspan=2>
__STOP__

$mysel="";
$mysel=" selected" if ($t_nrel);

$OUT .= <<__STOP__;
<b>Relatvistic correction:</b>
<select name="t_nrel">
<option value="">Use relativistic correction (recommended)
<option valut="NREL" $mysel>do NOT use relativistic correction
</select>
</td></tr>

<tr><td bgcolor=$gray2 colspan=2>
<b>Q-grid:</b>
$t_qgrid1
<select name="t_qgrid">
<option value="">----
__STOP__

$mysel="";
$mysel=" selected" if ($t_qgrid1 =~ /^U/);
    $OUT.="<option value=\"U\" $mysel>uniform";
$mysel="";
$mysel=" selected" if ($t_qgrid1 =~ /^L/);
    $OUT.="<option value=\"L\" $mysel>logarithmic";
$mysel="";
$mysel=" selected" if ($t_qgrid1 =~ /^1/);
    $OUT.="<option value=\"1\" $mysel>1-dimensional";

$OUT .= <<__STOP__;
</select>
Th0=<input name="t_qgrid2" value="$t_qgrid2" size=6>
(not used for uniform grid)
<br>
</td></tr>


<tr><td bgcolor=$gray2>
<b>Selection rule:</b>
__STOP__



my $cm=" selected" if ($t_srule1 =~ /m/);
my $cd=" selected" if ($t_srule1 =~ /p/);
my $cq=" selected" if ($t_srule1 =~ /q/);
my $co=" selected" if ($t_srule1 =~ /o/);
my $cn=" selected" if ($t_srule1 =~ /n/);
my $c1=" selected" if ($t_srule1 =~ /1/);
my $c2=" selected" if ($t_srule1 =~ /2/);



$OUT.=<<__STOP__;
<select name="t_srule">
<option value="">dipole (default)
<option value="m" $cm>monopole
<option value="d" $cd>dipole
<option value="q" $cq>quadrupole
<option value="o" $co>octopole
<option value="n" $cn>all transitions
<option value="1" $c1>mono+dipole
<option value="2" $c2>mono+di+quadrupole

</select>
</td>
<td bgcolor=$gray2>
<b>L(final)-selection rule:</b>
__STOP__

$cm=$cd=$cq=$co=$cn=$c1=$c2="";

$cm=" selected" if ($t_lrule1 =~ /m/);
$cd=" selected" if ($t_lrule1 =~ /p/);
$cq=" selected" if ($t_lrule1 =~ /q/);
$co=" selected" if ($t_lrule1 =~ /o/);
$cn=" selected" if ($t_lrule1 =~ /n/);
$c1=" selected" if ($t_lrule1 =~ /1/);
$c2=" selected" if ($t_lrule1 =~ /2/);

$OUT.=<<__STOP__;
<select name="t_lrule">
<option value="">----
<option value="m" $cm>L=l
<option value="d" $cd>L=l +/- 1
<option value="q" $cq>L=l +/- 2
<option value="o" $co>L=l +/- 3
<option value="n" $cn>all final states
<option value="1" $c1>|L-l| &le; 1
<option value="2" $c2>|L-l| &le; 2
</select>

</td></tr>


__STOP__

}

$OUT .= "<tr><td bgcolor=$gray><input type=submit value=save></FORM>";
$OUT .= "</td></tr></table>";

&PrintPage("Context",$OUT);


