#!/usr/bin/perl
 
require "../../libs/w2web.pl";
#$spinpol=0;
&GetInput();
&GetSession();
$prefspace="scf";
&GetPrefs();
#$spinpol=1 if ($PREFS{'spinpol'});
&ExeTypes();

$mytest="$DIR/$CASE.broyd1";
if ( -e $mytest && !$doit) {
  $OUT .= <<__STOP__;
<H2>Broyden files exist from previous scf-cycle!</H2>
<p>
When you have made changes to  input / struct / kmesh you should most likely "save" the previous calculation or at least remove the "broyden"-files.
</p>
<p> Do you want to
</p>
<ul>
<li><A HREF="/util/savelapw.pl?SID=$SID">Save Calculation with save_lapw</A>
</ul>
<ul>
<li><A HREF="/util/delete.pl?SID=$SID&file=$DIR/$CASE.broyd*">Remove the
files $CASE.broyd[1|2] </A>
</ul>
<p>or </p>
<ul>
<li> <A HREF="/exec/scf.pl?SID=$SID&doit=1">run anyway...</A> (only if you want to continue with the current scf-cycle)
</ul>
__STOP__
} else {

$OUT .=  <<__STOP__;
<H2>SCF Cycle</H2>

<FORM ACTION=/exec/scfrun.pl METHOD=POST>
__STOP__

&PassHiddenParms;


$OUT .=  <<__STOP__;
<TABLE><TR valign=TOP><TD bgcolor=$gray>
<fieldset>
<H3>Options: (<INPUT NAME=h TYPE=CHECKBOX >help)</H3> 
<table><tr valign=top><td>
<INPUT NAME=spinpol TYPE=CHECKBOX $spinpol>spin&nbsp;polarized<br>
<INPUT NAME=afm1 TYPE=CHECKBOX $afm>AFM calc.<br>
<INPUT NAME=it0 TYPE=CHECKBOX $PREFS{'it0'}>iterative(0)<br>
<INPUT NAME=in1orig TYPE=CHECKBOX $PREFS{'in1orig'}>in1orig<br>
<INPUT NAME=renorm TYPE=CHECKBOX $PREFS{'renorm'}>renorm rho<br>
<INPUT NAME=p1 TYPE=CHECKBOX $p>parallel<br>
<INPUT NAME=so TYPE=CHECKBOX $PREFS{'so'}>spinorbit<br>
<INPUT NAME=dm TYPE=CHECKBOX $PREFS{'dm'}>dm<br>
<INPUT NAME=orb TYPE=CHECKBOX $PREFS{'orb'}>orbital pot<br>
</td><td>
<INPUT NAME=fsm TYPE=CHECKBOX $PREFS{'fsm'}>FSM <INPUT NAME=fsm_count
SIZE=3 VALUE=$PREFS{'fsm_count'}><br>
<INPUT NAME=nohns TYPE=CHECKBOX $PREFS{'nohns'}>no HNS <INPUT NAME=nohns_count
SIZE=2 VALUE=$PREFS{'nohns_count'}><br>
<INPUT NAME=it TYPE=CHECKBOX $PREFS{'it'}>iterative diag  <INPUT NAME=it_count SIZE=2 VALUE=$PREFS{'it_count'}><br>
<INPUT NAME=in1new TYPE=CHECKBOX $PREFS{'in1new'}>in1new <INPUT NAME=in1new_count SIZE=2 VALUE=$PREFS{'in1new_count'}><br>
<INPUT NAME=ql TYPE=CHECKBOX $PREFS{'ql'}>q-limit <INPUT NAME=ql_count
SIZE=4 VALUE=$PREFS{'ql_count'}><br>
<INPUT NAME=itnum TYPE=CHECKBOX $PREFS{'itnum'}>It-number <INPUT NAME=itnum_count
SIZE=1 VALUE=$PREFS{'itnum_count'}><br>
</td></tr></table>
</fieldset>

</td>
<TD bgcolor=$gray>

<fieldset>
<H3>Expert options:</H3>
<INPUT NAME=expert SIZE=22>
</fieldset>

<br>
<fieldset>
<H3>Convergence criteria:</H3>
<table><tr valign=top><td>
<input name="conv_ec" type=checkbox  $PREFS{'conv_ec'}> Energy:
</td><td>
<input name="convval_ec" size=8 value="$PREFS{'convval_ec'}">Ry
</td></tr>
<tr valign=top><td>
<input name="conv_fc" type=checkbox  $PREFS{'conv_fc'}> Force:
</td><td>  
<input name="convval_fc" size=8 value="$PREFS{'convval_fc'}">mRy/au
</td></tr>
<tr valign=top><td>
<input name="conv_cc" type=checkbox  $PREFS{'conv_cc'}> Charge:
</td><td>
<input name="convval_cc" size=8 value="$PREFS{'convval_cc'}">e
</td></tr>
</td></tr></table>
</fieldset>
</td></tr>

<tr><td colspan=2 bgcolor=$gray>
$exetypes
</td></tr>
<tr><td colspan=2 bgcolor=$gray>
E-mail notification <INPUT TYPE=CHECKBOX NAME="notify" $PREFS{'notify'}>
to <INPUT NAME="email" SIZE=30 VALUE="$notification">
</td></tr>
<tr><td colspan=1 bgcolor=$gray>
<INPUT TYPE=SUBMIT VALUE="start SCF cycle">
<INPUT TYPE=RESET VALUE="Clear entries">
</td>
<td bgcolor=$red>
<form action=saveparam.pl>
<input type=checkbox name=saveonly>
only save parameters
</td></tr>
</form>
</table>

__STOP__
}

PrintPage("Context",$OUT);


#<tr><td colspan=2 bgcolor=$gray>
#<b>Expert options</b>
#<INPUT NAME=expert SIZE=40>
#</td></tr>
