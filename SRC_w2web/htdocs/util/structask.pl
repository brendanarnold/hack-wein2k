#!/usr/bin/perl

require "../../libs/w2web.pl";
require "../../libs/struct.pl";
$debug=0;
&GetInput;
&GetSession;

$upload ="ul&nbsp;";
$upload = "<img src=/art/ul.gif BORDER=0 ALT=\"upload\">";

if ($doit) {
# check for cif file
	$file="$DIR/$ciffil";
#        $umps = qx(echo $file > $DIR/mistfile);
	if ( -e $file && -s $file && ($ciffil ne "") ) {  
		$umps = qx( cd $DIR; cp $ciffil $CASE.cif; $WIENROOT/cif2struct $CASE.cif > $CASE.outputcif2struct);
		$numatoms=qx(cd $DIR; grep ATOMS $CASE.struct | cut -c28-30);
        } else {
	# generate struct file...
	$umps = qx( cd $DIR;cp $WIENROOT/SRC_templates/case.struct $CASE.struct);
        }

# immediately copy it to _i so we directly start in edit mode
	$structfile = "$DIR/$CASE.struct";
	$testfile = "$DIR/$CASE.struct_i";
	$umps = qx( cd $DIR;cp $structfile $testfile );
	&StructRead;
#set s_lattice, otherwise spacegroup will not be called
        $s_spacegr=~ s/ /_/;
	$s_nato=$numatoms;
        foreach $v (@lattype) {
          if ( $s_spacegr =~ /$v/ ) {
		$s_lattice = $s_spacegr;
	  } 
        }
        for ($i = 1; $i <= $s_nato; $i++) {
		$s_npt[$i]=781;
		$s_r0[$i]=0.0001;
		$s_rmt[$i]=2.0;
		$s_zz[$i]="";
		$s_mult[$i] = 1;
		$s_isplit[$i] = 8;
	}
	&StructWrite;
	redirectURL("/util/structgen.pl?SID=$SID");
} else {
	$OUT .=  <<__STOP__;
<H2>StructGen<font size=-2><sup>TM</sup></font></H2>

<FORM ACTION="/util/structask.pl" METHOD=POST>
<INPUT TYPE=HIDDEN NAME="doit" VALUE=1>
__STOP__
	&PassHiddenParms;
        $ciffil=none;
	$OUT .=  <<__STOP__;
<h3>You do not have a $CASE.struct file yet. </h3>
<p class="info">
<b>You can create it using STRUCTGEN. Please specify the number of independent atoms of your  initial structure!</b>
<br><br>
Number of atoms: 
<INPUT NAME="numatoms" VALUE=2 SIZE=4>
<INPUT TYPE=SUBMIT VALUE="Generate template">
</p>

<p>
<b>Alternatively you can convert a "cif" file: </b> 
(e.g. from the Inorganic crystal structure database). 
</p>
__STOP__

    # only show these lines if we have cif-files
    @ciffiles = qx(cd $DIR; /bin/ls -c *.cif);
    if (@ciffiles) {
	$OUT.="<p>Select one of the following \"cif\" files:<br>";
	foreach $i (@ciffiles) {
        $ii = $i;
	chop($ii);
	$OUT .=  <<__STOP__;
$i <INPUT TYPE=RADIO NAME="ciffil" VALUE="$ii" ><br>
__STOP__
}
$OUT.="</p><p><INPUT TYPE=SUBMIT VALUE=\"Use selected CIF file\">";
}
	$OUT .=  <<__STOP__;
<p>Here you can  <a href=\"/util/upload.pl?SID=$SID&cif=1\">$upload upload </a> a "cif" file from your local computer. 
</p>
</form>
__STOP__
	&PrintPage("Context",$OUT);
}
