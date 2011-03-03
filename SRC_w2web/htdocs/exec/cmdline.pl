#!/usr/bin/perl

require "../../libs/w2web.pl";
&GetInput;
&GetSession;
&ExeTypes();

$OUT .=  <<__STOP__;
<H2>Execute a command line:</H2>

<FORM ACTION=/exec/executor.pl METHOD=POST>
__STOP__

&PassHiddenParms;

$OUT .=  <<__STOP__;
<INPUT NAME=prog TYPE=INPUT SIZE=50>
<BR>
<INPUT TYPE=SUBMIT VALUE=" Run ">
$exetypes2
</select>
</FORM>
__STOP__

PrintPage("Context",$OUT);


