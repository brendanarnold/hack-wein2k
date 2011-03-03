#!/usr/bin/perl

require "../../libs/w2web.pl";
require "../../libs/struct.pl";
$debug=0;
&GetInput;
&GetSession;
	$generateeq = 1 if (length $s_lattice > 3);

&StructWrite;

$OUT .=  <<__STOP__;
<H2>StructGen<font size=-2><sup>TM</sup></font></H2>

<H3>
<ul>

<li><A HREF="structrmt.pl?SID=$SID">set automatically RMT and continue editing</A> (do it at least once!)
<br>
<br>
<li><A HREF="structend.pl?SID=$SID">save file and clean up</A> (when you are done)
<br>
<br>
<li><A HREF="structgen.pl?SID=$SID">continue editing</A>
<br>
<br>
<li><A HREF="structabort.pl?SID=$SID">abort editing and restore 
original file</A>
<br>
</H3>

__STOP__




PrintPage("Context",$OUT);


