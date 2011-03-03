#!/usr/bin/perl

require "../../libs/w2web.pl";
require "../../libs/struct.pl";
require "../../libs/telnes2.pl";


#$debug=1;
&GetInput;
&GetSession;

#correct a few things:
if($t_edge) {
    if ($t_edge =~ /K/) {$t_n=1;$t_l=0;}
    if ($t_edge =~ /L1/) {$t_n=2;$t_l=0;}
    if ($t_edge =~ /L23/) {$t_n=2;$t_l=1;}
}

$t_output1=$t_output;
if ($t_atoms1 || $t_atom2) {$t_atoms=1;}
$t_modus1=$t_modus;

if ($t_egrid1 || $t_egrid2 || $t_egrid3) {$t_egrid=1;}

if ($t_det1 || $t_det2) {$t_det=1;}
$t_init=0;
@ilist=qw(t_init1 t_init2 t_init3 t_init4);
    foreach $i (@ilist) {
	if ($$i =~ /on/){
	    $$i = "y";
	    $t_init++;
	} else {
	    $$i = "n";
	}
    }
$t_branch=1 if $t_branch1;
$t_split=1 if $t_split1;


$t_qgrid1=$t_qgrid;
$t_srule1=$t_srule;
$t_lrule1=$t_lrule;

$OUT .=  <<__STOP__;
<H2>InnesGen<font size=-2><sup>TM</sup></font> for TELNES.2</H2>

<p class=red>$CASE.innes saved</p>

__STOP__

&InnesWrite;

if($next) {
    $OUT .= <<__STOP__;
<FORM ACTION="/exec/next.pl" METHOD=POST>
<INPUT TYPE=HIDDEN NAME=nexturl VALUE="$nexturl">
<INPUT TYPE=submit VALUE="$FORM{'next'}">
</FORM>
__STOP__
}

PrintPage("Context",$OUT);


