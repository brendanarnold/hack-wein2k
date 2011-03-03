#!/usr/bin/perl

require "../libs/w2web.pl";

my $indent = q/<br>&nbsp;/;

&GetInput;
&GetSession;

$OUT .=  <<__STOP__;
Content-type: text/html\n
<HTML>
<HEAD>
   <TITLE>WIEN</TITLE>
<link href="$css" rel="stylesheet">
__STOP__

if ($jsmenu) {
	$OUT.= <<__STOP__;
<style type="text/css">
div.clSlideMenu{ /* All slideMenu2 DIV elements */
  position:absolute;
  font-family:helvetica,arial,sans-serif;
  font-size:10px;
  overflow:hidden;
  width:150;
  height:22;
}
a.clA0{ /* All top level links */
  color:white;
  font-size:12px;
  text-decoration:none;
  font-weight:bold;
}
a.clA1{ /* All sub level links */
  color:black;
  font-size:11px;
  font-weight:bold;
  text-decoration:none;
}
a.clA0:before {content: "";}
a.clA0:after {content: "";}
a.clA1:before {content: "";}
a.clA1:after {content: "";}

div.slideMenuBG{position:absolute; left:0; top:0; z-index:1}
div.slideMenuText{position:absolute; left:2; top:2; z-index:200;}
/* NEEDED ----- This class should be named like this:
  #NAME_OF_YOUR_MENUcont
*/
#slideMenucont{position:absolute; height:600; width:200; visibility:hidden;}
</style>
<script language="JavaScript" src="slidemenu.js" type="text/javascript">
</script>

__STOP__
}

if ($platinum) {
$jstop=80;
$OUT.= <<__STOP__;
</HEAD>
<BODY>
<div class="main">
<img src="/art/platinum.gif" width=132 height=68 alt="w2web platinum edition" border=0>
</div>
<br>
__STOP__
} else {

$jstop=150;
$OUT.= <<__STOP__;
</HEAD>
<BODY BGCOLOR=$topcolor>
<A HREF="http://www.wien2k.at" TARGET="main" class="img"><IMG SRC="/art/w2k.gif" WIDTH=124 HEIGHT=139 ALT="WIEN2k Homepage" BORDER=0></A>
<FONT SIZE=-1>
__STOP__
}

if($jsmenu) {
$xcrys="";
if ($ENV{'XCRYSDEN_TOPDIR'}) {
  $xcrys="slideMenu.makeMenu('sub','view structure','/util/viewxcrys.pl?SID=$SID','main')";
}
	$OUT .= <<__STOP__;
<script>
slideMenu = new createSlideMenu("slideMenu")

//Variables to set:
slideMenu.menuy=$jstop //The top placement of the menu.
slideMenu.menux=0 //The left placement of the menu
slideMenu.useImages = 1 //Are you using images or not?
slideMenu.pxspeed=40 //The pixel speed of the animation
slideMenu.timspeed=25 //The timer speed of the animation
slideMenu.inset = 10 //How much the selected items should pop to the left
slideMenu.arrow = 0 //Set this to className that has font-family:webdings
                    //if you want to use the arrow feature. Note:
                    //This ONLY works on DOM capable browsers, and with
                    //useImages set to 0 - It's basically just a test I
                    //did.
                    //I hope to improve it later on.

//Needed dummy classes - leave in the stylesheet!
slideMenu.bgClass = "slideMenuBG"
slideMenu.txtClass = "slideMenuText"

/*******************************************************************************
Level properties - ALL properties have to be specified in level 0
This works the same way as the CM4 script (if you have used it)

The level[0] values will be the default value. Any value not specified
in higher levels will be inherited from level[0]. If anything
is spesified in level[1], but not in level[2], level[2] (sub2 menus)
will get the property value from level[1] and so on.

The availble values to control the menu by level are:

left           - The left placement of all items in the current level ( px
value )
width          - The width of all items in the current level  ( px value )
height         - The height of all items in the current level  ( px value )
between        - The number of pixels between each item in  the current
level ( px value)
className      - A name of a class specified in the stylesheet to control
the
                 look of all items in this level. 
                 NOTE: The class MUST be in a stylesheet, and it most have
position:absolute.
classNameA     - A name of a class specified in the stylesheet that will
control the 
                 Look of the TEXT for all items in this level. (you can
also specify 
                 a hover class for this className to get a mouseover effect
on the */

slideMenu.level[0] = new slideMenu_makeLevel(
  left = 0,
  width = 138,
  height = 21,
  between = 5,
  className = "clSlideMenu",
  classNameA = "clA0",
  regImage = "art/level0_regular.gif",
  roundImg = "art/level0_round.gif",
  roundImg2 = "",
  subImg = "",
  subRound= "")
  
slideMenu.level[1] = new
slideMenu_makeLevel(10,127,20,2,"clSlideMenu","clA1","art/level1_regular.gif","art/level1_round2.gif","art/level1_round.gif","art/level1_sub.gif",
"level1_sub_round.gif")

//Image preload --- leave this
for(var i=0;i<slideMenu.level;i++){
  var l = slideMenu.level[i]
  new
preLoadBackgrounds(l.regImage,l.roundImg,l.roundImg2,l.subImg,l.subRound)
}

//Menu 1 -----------------------
slideMenu.makeMenu('top','Excecution')
  slideMenu.makeMenu('sub','StructGen','/util/structgen.pl?SID=$SID','main')
	$xcrys
  slideMenu.makeMenu('sub','initialize calc.','/exec/initlapw.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','run SCF','/exec/scf.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','single prog.','/exec/single.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','optimize (V, c/a)','/exec/optimize.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','mini. positions','/exec/min.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','command line','/exec/cmdline.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','frozen phonons','/exec/fphonons.pl?SID=$SID','main')
slideMenu.makeMenu('top','Utilities')
  slideMenu.makeMenu('sub','show dayfile','/util/dayfile.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','show STDOUT','/util/stdout.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','analysis','/util/analysis.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','save_lapw','/util/savelapw.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','restore_lapw','/util/restorelapw.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','initso_lapw','/util/initso.pl?SID=$SID','main')
	$xcrys
  slideMenu.makeMenu('sub','stop SCF','/util/stopscf.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','stop mini','/util/stopmini.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','full diag.','/util/fulldiag.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','edit .machines','/util/editpara.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','testpara','/util/testpara.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','testpara1','/util/testpara1.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','testpara2','/util/testpara2.pl?SID=$SID','main')
slideMenu.makeMenu('top','Tasks')
  slideMenu.makeMenu('sub','El. Dens.','/exec/rho.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','DOS','/exec/dos.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','XSPEC','/exec/xspec.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','TELNES.2','/exec/telnes2.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','OPTIC','/exec/optic.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','Bandstructure','/exec/band.pl?SID=$SID','main')
slideMenu.makeMenu('top','Files')
  slideMenu.makeMenu('sub','struct file(s)','/util/dir.pl?SID=$SID&dir=$DIR&ext=struct','main')
  slideMenu.makeMenu('sub','input file(s)','/util/dir.pl?SID=$SID&dir=$DIR&ext=in','main')
  slideMenu.makeMenu('sub','output file(s)','/util/dir.pl?SID=$SID&dir=$DIR&ext=out','main')
  slideMenu.makeMenu('sub','SCF file(s)','/util/dir.pl?SID=$SID&dir=$DIR&ext=scf','main')
  slideMenu.makeMenu('sub','CLM file(s)','/util/dir.pl?SID=$SID&dir=$DIR&ext=clm','main')
  slideMenu.makeMenu('sub','def file(s)','/util/dir.pl?SID=$SID&dir=$DIR&ext=def','main')
  slideMenu.makeMenu('sub','error file(s)','/util/dir.pl?SID=$SID&dir=$DIR&ext=error','main')
  slideMenu.makeMenu('sub','running file(s)','/util/dir.pl?SID=$SID&dir=$DIR&ext=run','main')
  slideMenu.makeMenu('sub','script file(s)','/util/dir.pl?SID=$SID&dir=$DIR&ext=sh','main')
  slideMenu.makeMenu('sub','upload file','/util/upload.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','new file','/util/newfile.pl?SID=$SID','main')
  slideMenu.makeMenu('sub','show all files','/util/dir.pl?SID=$SID&dir=$DIR&f=1','main')
slideMenu.makeMenu('top','Session Mgmt.')
  slideMenu.makeMenu('sub','change session','$ENV{'MASTER_URL'}/session/change.cgi?SID=$SID','_parent')
  slideMenu.makeMenu('sub','change dir.','/util/dir.pl?SID=$SID&dir=$DIR&cd=1','main')
  slideMenu.makeMenu('sub','change info','/session/changename.cgi?SID=$SID','main')
  slideMenu.makeMenu('sub','duplicate session','/session/duplicate.cgi?SID=$SID','main')
  slideMenu.makeMenu('sub','delete session','/session/delete.cgi?SID=$SID','main')
slideMenu.makeMenu('top','Configuration')
  slideMenu.makeMenu('sub','User config','/config/index.pl?SID=$SID','main')
slideMenu.makeMenu('top','Usersguide')
  slideMenu.makeMenu('sub','HTML Version','/usersguide/usersguide.html','_new')
  slideMenu.makeMenu('sub','PDF Version','/usersguide/usersguide.pdf','_new')
  slideMenu.makeMenu('sub','WIEN2k website','http://www.wien2k.at','_new')
__STOP__

$OUT.="slideMenu.makeMenu('sub','online expert','/util/expert.pl','main')" if $elnesexpert;

$OUT.=<<__STOP__;
//Initiating the menu !! 
slideMenu.init()    
</script>
__STOP__

} else {
if($exec) {
	$OUT .= <<__STOP__;
<A NAME=a>
<br>
<A HREF="/navig.pl?SID=$SID&exec=0"><B>&lt;&lt; Execution </B></A>
__STOP__
} else {
	$OUT .= <<__STOP__;
<A HREF="/navig.pl?SID=$SID&exec=1"><B>Execution &gt;&gt; </B></A>
__STOP__
}
	$OUT .= <<__STOP__;
$indent<A HREF=/util/structgen.pl?SID=$SID TARGET="main">StructGen<font size=-2><sup>TM</sup></font></A> 
$indent<A HREF=/exec/initlapw.pl?SID=$SID TARGET="main">initialize calc.</A>
$indent<A HREF=/exec/scf.pl?SID=$SID TARGET="main">run SCF</A>
$indent<A HREF=/exec/single.pl?SID=$SID TARGET="main">single prog.</A>
$indent<A HREF=/exec/optimize.pl?SID=$SID TARGET="main">optimize(V,c/a)</A>
$indent<A HREF=/exec/min.pl?SID=$SID TARGET="main">mini. positions</A>
__STOP__
if($exec) {
  $OUT .= <<__STOP__;
$indent<A HREF=/exec/cmdline.pl?SID=$SID TARGET="main">command line</A>
$indent<A HREF=/exec/fphonons.pl?SID=$SID TARGET="main">frozen phonons</A>
__STOP__
}
$OUT .= <<__STOP__;
<BR>
<BR>
__STOP__

$xcrys="";
if ($ENV{'XCRYSDEN_TOPDIR'}) {
	$xcrys="$indent<A HREF=/util/viewxcrys.pl?SID=$SID TARGET=\"main\">view structure</A>";
}
if ($utils) {
$OUT .= <<__STOP__;
<A HREF="/navig.pl?SID=$SID&utils=0"><B>&lt;&lt; Utils. </B></A>
$indent<A HREF=/util/dayfile.pl?SID=$SID TARGET="main">show dayfile</A>
$indent<A HREF=/util/stdout.pl?SID=$SID TARGET="main">show STDOUT</A>
$indent<A HREF=/util/analysis.pl?SID=$SID TARGET="main">analysis</A>
$indent<A HREF=/util/savelapw.pl?SID=$SID TARGET="main">save_lapw</A>
$indent<A HREF=/util/restorelapw.pl?SID=$SID TARGET="main">restore_lapw</A>
$indent<A HREF=/util/cleanlapw.pl?SID=$SID TARGET="main">clean_lapw</A>
$indent<A HREF=/util/initso.pl?SID=$SID TARGET="main">initso_lapw</A>
$xcrys
$indent<A HREF=/util/stopscf.pl?SID=$SID TARGET="main">stop SCF</A>
$indent<A HREF=/util/stopmini.pl?SID=$SID TARGET="main">stop mini</A>
$indent<A HREF=/util/fulldiag.pl?SID=$SID TARGET="main">full diag.</A>
$indent<A HREF=/util/editpara.pl?SID=$SID TARGET="main">edit .machines</A>
$indent<A HREF=/util/testpara.pl?SID=$SID TARGET="main">testpara</A>
$indent<A HREF=/util/testpara1.pl?SID=$SID TARGET="main">testpara1</A>
$indent<A HREF=/util/testpara2.pl?SID=$SID TARGET="main">testpara2</A>

__STOP__
} else {
	$OUT .= <<__STOP__;
<A HREF="/navig.pl?SID=$SID&utils=1"><B>Utils. &gt;&gt;</B></A>
__STOP__
}


$OUT .= <<__STOP__;
<BR>
<BR>
__STOP__
if ($tasks) {
$OUT .= <<__STOP__;
<A HREF="/navig.pl?SID=$SID&tasks=0"><B>&lt;&lt; Tasks</B></A>
$indent<A HREF=/exec/rho.pl?SID=$SID TARGET="main">El. Dens.</A>
$indent<A HREF=/exec/dos.pl?SID=$SID TARGET="main">DOS</A>
$indent<A HREF=/exec/xspec.pl?SID=$SID TARGET="main">XSPEC</A>
$indent<A HREF=/exec/telnes2.pl?SID=$SID TARGET="main">TELNES.2</A>
$indent<A HREF=/exec/optic.pl?SID=$SID TARGET="main">OPTIC</A>
$indent<A HREF=/exec/band.pl?SID=$SID TARGET="main">Bandstructure</A>
__STOP__
} else {
$OUT .= <<__STOP__;
<A HREF="/navig.pl?SID=$SID&tasks=1"><B>Tasks &gt;&gt;</B></A>
__STOP__
}
$OUT .= <<__STOP__;
<br>
<br>
__STOP__
if ($files) {
$OUT .= <<__STOP__;
<A HREF="/navig.pl?SID=$SID&files=0"><b>&lt;&lt; Files</b></A>
__STOP__
} else {
$OUT .= <<__STOP__;
<A HREF="/navig.pl?SID=$SID&files=1"><b>Files &gt;&gt;</b></A>
__STOP__
}
$OUT .= <<__STOP__;
$indent<A HREF=/util/dir.pl?SID=$SID&dir=$DIR&ext=struct TARGET="main">struct file(s)</A>     
$indent<A HREF=/util/dir.pl?SID=$SID&dir=$DIR&ext=in TARGET="main">input files</A>     
$indent<A HREF=/util/dir.pl?SID=$SID&dir=$DIR&ext=out TARGET="main">output files</A> 
$indent<A HREF=/util/dir.pl?SID=$SID&dir=$DIR&ext=scf TARGET="main">SCF files</A> 
__STOP__
if ($files) {
$OUT .= <<__STOP__;
$indent<A HREF=/util/dir.pl?SID=$SID&dir=$DIR&ext=clm TARGET="main">CLM files</A> 
$indent<A HREF=/util/dir.pl?SID=$SID&dir=$DIR&ext=def TARGET="main">def files</A> 
$indent<A HREF=/util/dir.pl?SID=$SID&dir=$DIR&ext=error TARGET="main">error files</A> 
$indent<A HREF=/util/dir.pl?SID=$SID&dir=$DIR&ext=run TARGET="main">running files</A> 
$indent<A HREF=/util/dir.pl?SID=$SID&dir=$DIR&ext=sh TARGET="main">script files</A> 
$indent<A HREF=/util/upload.pl?SID=$SID TARGET="main">upload file</A> 
$indent<A HREF=/util/newfile.pl?SID=$SID TARGET="main">new file</A> 
<br>
$indent<A HREF=/util/dir.pl?SID=$SID&dir=$DIR&f=1 TARGET="main">show all files</A>     
<br>
__STOP__
}


$OUT .= <<__STOP__;
<br>
<br>
__STOP__
if ($FORM{'session'}) {
$OUT .= <<__STOP__;
<A HREF="/navig.pl?SID=$SID&session=0"><B>&lt;&lt; Session Mgmt.</B></A>
__STOP__
} else {
$OUT .= <<__STOP__;
<A HREF="/navig.pl?SID=$SID&session=1"><B>Session Mgmt. &gt;&gt;</B></A>
__STOP__
}

$OUT .= <<__STOP__;
$indent<A HREF=$ENV{'MASTER_URL'}/session/change.cgi?SID=$SID TARGET=_parent>change session</A>
$indent<A HREF=/util/dir.pl?SID=$SID&dir=$DIR&cd=1 TARGET="main">change dir</A>
$indent<A HREF=/session/changename.cgi?SID=$SID TARGET="main">change info</A>
__STOP__

if ($FORM{'session'}) {
$OUT .= <<__STOP__;
<!--
$indent<A HREF=/session/duplicate.cgi?SID=$SID TARGET="main">duplicate session</A>
-->
$indent<A HREF=/session/delete.cgi?SID=$SID TARGET="main">delete session</A>
<br>
__STOP__
}

$OUT .= <<__STOP__;
<br>
<br>
<A HREF="/config/index.pl?SID=$SID" TARGET="main"><B>Configuration </B></A>
__STOP__


$OUT .= <<__STOP__;
<br>
<br>
<B>Usersguide</B></A>
$indent<A HREF="/usersguide/usersguide.html" TARGET="_new">html-Version</A>
$indent<A HREF="/usersguide/usersguide.pdf" TARGET="_new">pdf-Version</A>
</FONT>
<br><br>
__STOP__
}

if ($platinum) {
  $w2k="Compute engine<br>powered by <a href=\"http://www.wien2k.at\" target=\"_blank\">WIEN2k</A>.<br><br><br>";
} else {
$w2k="";
}

$OUT .=<<__STOP__;

<div class="copy">
$w2k
Idea and realization by<br>
<A HREF="http://www.luitz.at/w2web" TARGET="_blank"><i>luitz.at</i></A>
&copy; 2001-2006
</div>

</BODY>
</HTML>
__STOP__

print $OUT;

