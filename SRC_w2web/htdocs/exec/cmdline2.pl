#!/usr/bin/perl

require "../../libs/w2web.pl";
&GetInput;
&GetSession;
&ExeTypes();

$OUT .=  <<__STOP__;
<p>
<FORM ACTION=/exec/executor.pl METHOD=POST TARGET="main">
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT NAME=prog TYPE=INPUT SIZE=50 value="command">
<INPUT TYPE=SUBMIT VALUE=" Run ">
$exetypes2
</SELECT>
</FORM>
</p>
__STOP__

PrintPage("Context",$OUT);


