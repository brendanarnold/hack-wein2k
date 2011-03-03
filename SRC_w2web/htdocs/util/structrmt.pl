#!/usr/bin/perl


require "../../libs/w2web.pl";
$debug=0;
&GetInput;
&GetSession;

if (!$doit) {

$OUT = <<__STOP__;
<h3>Automatic determination of RMTs</h3>

Please specify the desired RMT reduction compared to almost touching spheres.
<br>Typically use:
<br><br>for a single calculation: 0 %
<br>for force minimization: 1-5 % 
<br>for volume effects you may need even larger reductions.
<br><br>


<FORM ACTION="/util/structrmt.pl" METHOD=POST>
<INPUT TYPE=HIDDEN NAME="doit" VALUE="1">
<INPUT TYPE=HIDDEN NAME="SID" VALUE="$SID">
<b>Reduce RMTs by </b> <INPUT NAME="reduc" VALUE="0" size=3 > % <br>
<br>
<INPUT TYPE=SUBMIT VALUE="do it">
</FORM>

__STOP__

} else {

$umps = qx(cp $DIR/$CASE.struct_i $DIR/setrmt.struct);
$umps = qx(cd $DIR;$WIENROOT/setrmt_lapw setrmt -r $reduc);
$umps = qx(cp $DIR/setrmt.struct_setrmt $DIR/$CASE.struct_i);


redirectURL("/util/structgen.pl?SID=$SID");

}


PrintPage("Set RMT", $OUT);
