#!/usr/bin/perl

require "../../libs/w2web.pl";
require "../../libs/struct.pl";

$SCFFILES="*";
&GetInput;
&GetSession;
&ExeTypes();

#The program will always save the current CASE.struct by CASE_initial.struct
#at the end it will restore the original CASE.struct by mv from CASE_initial.struct!

#values for doit:
#0 default: check for broyden files and ask for their removal
#1 display menu and run fphonons anyway
#2 display the input options for phonon modes and the distortion parameters
#3 LAYOUT plot mode (just input, nothing will be plotted right now!)
#4 generate new struct files from the input
#5 LAYOUT save current CASE.mode as CASE_***.mode (just specify the name, nothing will be saved!)
#6 execute of 5 (save files)
#7 execute of 3 (plot the graph)

#first we check the existing .mode files

$MODEFILES="";
unless(opendir(VERZ,"$DIR"))
  {
    &CGIError ("Could not read from directory!");
    exit;
  }
while($line=readdir(VERZ))  
  {
  if($line =~/(.*?)(\.mode)/i)
    {
      $x=$1;
      chomp($x);
      $selected="";
      if($INPUT1 eq "$DIR/$x")
      {
        $selected="selected";
      }
      $MODEFILES.="<OPTION value=\"$DIR/$x$2\" $selected>$x$2</OPTION>";
    }
  }
closedir(VERZ);

if(!$doit)
{
  #initial state
  $doit=0;
}

if($doit==0)
{
  #checking for broyden files.
  $mytest="$DIR/$CASE.broyd1";
  if ( -e $mytest)
    {
    $OUT .= <<__STOP__;
<H2>Broyden files exist from previous scf cycle!</H2>PUT NAME="afm" VALUE="" TYPE=HIDDEN>
<INPUT NAME="COMPLEX" VALUE="" TYPE=HIDDEN>
<INPUT NAME="p" VALUE="" TYPE=HIDDEN>
<INPUT NAME="COMMENT" VALUE="" TYPE=HIDDEN>
Title: <INPUT NAME="s_title"  VALUE="TiC" SIZE=40>
<br>
Treatment of core: 
<SELECT NAME="s_rels"><OPTION VALUE=RELA SELECTED>relativistic
<OPTION VALUE=NREL >non relativistic
</SELECT>

<br>
<br>
<p>
You should most likely "save" this calculation before running frozen phonon calculation.
</p>
<p>
Do you want to
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
<li> <A HREF="/exec/fphonons.pl?SID=$SID&doit=1">run anyway...</A> (very unlikely, that this is a good choice)
</ul>
__STOP__
    }
else
  {
  $doit=1;
  #if no broyden Files were detected we go to the next step
  }
}

if($doit==1)
{
  $OUT .= "<h2>Introduce Phonon Modes</h2>";
  $next="continue with frozen phonon calculation";
  $nexturl="/exec/fphonons.pl?SID=$SID";
  $nextinteractive=1;
  $xxout = "${CASE}.struct";
  $disabled1="disabled";
  $disabled2="disabled";
  $disabled3="disabled";
  $NEWID=$SID;#is overwritten by the ID stored in __NEWID__
  $changebutton="<INPUT TYPE=BUTTON VALUE=\"change and init session\" disabled>";
  if (-e "$DIR/fphonons.job")
  {
    $disabled1="";
  }
  if (-e "$DIR/__NEWID__")
  {
    unless(open(FILE,"$DIR/__NEWID__"))
    {
      &CGIError("Can't read file $fname.\n");
      exit;
    }
    @lines = <FILE>;
    close(FILE);
    chomp(@lines);
    $NEWID=$lines[0];#the first entry!
    if($NEWID<1||$NEWID>999999)
      {$NEWID=$SID;}
  $changebutton="<a href=\"../session/select.cgi?SID=$NEWID\" TARGET=_parent><INPUT TYPE=BUTTON VALUE=\"change and init session\" onclick=\"parent.location.href='../session/select.cgi?SID=".$NEWID."'\"></a>";
  }
  if (-e "$DIR/${CASE}.mode")
  {
    $disabled2="";
  }
  if ($MODEFILES ne "")
  {
    $disabled3="";
  }
  $OUT .=  <<__STOP__;
<TABLE BGCOLOR=$green>

<TR><TD BGCOLOR=$gray>
<FORM ACTION=/exec/fphonons.pl METHOD=POST>
__STOP__
&PassHiddenParms();
 $OUT .=  <<__STOP__;
<INPUT TYPE=SUBMIT VALUE="init fphonons">
Generate structure files from $xxout
<INPUT type=checkbox NAME=break_symm VALUE=1 checked> break symmetry
<INPUT TYPE=hidden NAME=nextinteractive VALUE=1>
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl&doit=2">
<INPUT TYPE=hidden NAME=doit VALUE="2">
</FORM>
</TD></TR>
</FORM>

<TR><TD BGCOLOR=$gray>
<FORM ACTION=/exec/fphonons.pl METHOD=POST>
__STOP__
&PassHiddenParms();
 $OUT .=  <<__STOP__;
<INPUT TYPE=SUBMIT VALUE="save data" $disabled2>
Save the current fphonon files to a new directory
<INPUT TYPE=hidden NAME=nextinteractive VALUE=1>
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl&doit=5">
<INPUT TYPE=hidden NAME=doit VALUE="5">
</FORM>
</TD></TR>

<TR><TD BGCOLOR=$gray>
<FORM ACTION="/util/edit.pl">
<INPUT NAME=SID VALUE=$SID TYPE=HIDDEN>
<INPUT NAME=file VALUE="$DIR/fphonons.job" TYPE=HIDDEN>
<INPUT TYPE=SUBMIT VALUE="edit fphonons.job" $disabled1>
    Uncomment "x dstart" or "cp clmsum"; change options in run_lapw, save_lapw,...
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl&doit=1">
<INPUT TYPE=hidden NAME=doit VALUE="1">
</FORM>
</TD></TR>

<TR><TD BGCOLOR=$gray>
<FORM NAME="quargl">
$changebutton
change session and init the calculation
</FORM>
</TD></TR>

<TR><TD BGCOLOR=$gray>
<FORM ACTION=/exec/executor.pl METHOD=post>
__STOP__
&PassHiddenParms;
#this part must be redirected to the executor!
$OUT .=  <<__STOP__;
<INPUT NAME=precmd TYPE=HIDDEN VALUE="">
<INPUT NAME=prog TYPE=HIDDEN VALUE="fphonons.job">
<INPUT TYPE=SUBMIT VALUE="run fphonons.job" $disabled1>
$exetypes<br>
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl&doit=2">
</FORM>
</TD></TR>

<TR><TD BGCOLOR=$gray>
<FORM ACTION=/exec/fphonons.pl METHOD=POST>
__STOP__
&PassHiddenParms;
$OUT .=  <<__STOP__;
<INPUT TYPE=hidden NAME=doit VALUE="3">
<INPUT TYPE=SUBMIT VALUE="plot" disabled>
Plot energy curve (not implemented yet)
</FORM>
</TD></TR>

<TR><TD BGCOLOR=$gray>
<FORM ACTION=/util/edit.pl METHOD=POST>
<INPUT TYPE=SUBMIT VALUE="edit mode" $disabled3>
select a .mode file
<SELECT name=file>$MODEFILES</SELECT>
<INPUT TYPE=hidden NAME=nextinteractive VALUE=1>
<INPUT TYPE=hidden NAME=next VALUE="$next">
<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl&doit=1">
<INPUT TYPE=hidden NAME=doit VALUE="1">
</FORM>
</TD></TR>

</table>
__STOP__
}

if ($doit==2)
{
    qx(rm $DIR/__NEWID__);#delete old SESSIONID
    $mode_file="$DIR/$CASE.struct";
    my $filename="$CASE.struct";
    if($INPUT1 ne '0' && $INPUT1 ne '')
    {
      $mode_file=$INPUT1;
      unless(open(FILE,"$mode_file"))
      {
        &CGIError("could not open $mode_file\n");
        exit;
      }
      @lines = <FILE>;
      close(FILE);
      $fphon="";
      &StructRead;#read the structure
      my $i,$j,$temp;
      $i=1;#init of atom counter
      $j=1;#init of multiplicity counter
      foreach(@lines)
      {
        $temp=$_;
        if ( $temp =~ s/^([a-z]{1,2} ?[^#]*?)#.*\n/$1/ig || $temp =~ s/^([a-z]{1,2} ?[^#]*?)\n/$1/ig )
	  {
            $comment="";
	    if($i>$s_nato || $j>$s_mult[$i])
	    {$fphon.=$temp."    #!!! warning, redundant line!";}
	    else
	    {
              $s_name[$i] =~ s/\W//g;
	      $s_nameadd[$i] =~ s/\W//g;
	      if($s_nameadd[$i] eq '')
	        {$regexp="^$s_name[$i] ";}
	      else
	        {$regexp="^$s_name[$i] $s_nameadd[$i] ";}
	      if($break_symm == 1)
	        {$regexp.="\\[".$j."\\]";}
	      if($1 !~ m/$regexp/i)
	      {
		$temp2=$temp;
		$brackets="  ";
		if($break_symm == 1)
		  {$brackets=" [$j]";}
		#$temp =~ s/([a-z]{1,2} *[0-9]* |[a-z]{1,2}[0-9]* )( *[\]\[0-9]* )? +/$s_name[$i] $s_nameadd[$i]$brackets /i;
		#$temp2 =~ /([a-z]{1,2} *[0-9]* |[a-z]{1,2}[0-9]* )( *[\]\[0-9]* )?/ig;
		$temp =~ s/([a-z]{1,2}\W*)([0-9]+[^.0-9]*\W*)?(\[[0-9]*\])? +/$s_name[$i] $s_nameadd[$i]$brackets /i;
		$temp2 =~ /([a-z]{1,2}\W*)([0-9]+[^.0-9]*\W*)?(\[[0-9]*\])? +/ig;
		$temp2=$1.$2;
		$comment=" !!! warning, originally $temp2";
	      }
	      $temp =~ s/\] +/] /ig;#to smoothen all lines
              if($break_symm != 1)
                {$temp =~ s/\[[0-9]*\]/ /ig;}
	      $fphon.=$temp."    #".sprintf("%1.6f",$s_x[$i][$j])."  ".sprintf("%1.6f",$s_y[$i][$j]);
	      $fphon.="  ".sprintf("%1.6f",$s_z[$i][$j])."$comment";
	    }
	    $fphon.="\n";
	    if($break_symm == 1)
	    {
		$j++;
		if($j>$s_mult[$i])
		{
		    $j=1;
		    $i++;
		}
	    }
	    else
	    {$i++;}
          }
        else
	  {
	    $fphon.=$temp;
	  }
      }
    }
    else
    {
    #before we modify anything, we save the original .struct file!
    #at the end of doit=2 the original will be restored. The modified one will be located in ${CASE}_brokensymm.struct
    qx(cp $DIR/$CASE.struct $DIR/${CASE}_initial.struct);
    if($break_symm == 1)
    {
      #run supercell and nn to get all atoms in the unit cell. Only centering symmetry is retained
      qx(cp $DIR/$CASE.struct $DIR/${CASE}_initial.struct); #to retain the original struct file
      unless(open(FILE,">$DIR/temp.tmp"))
        {
        &CGIError("Context","could not create input file\n");
        exit;
        }
      print FILE "$filename\n1\n1\n1\nP\n0\n0\n0";
      close(FILE);
      qx(cp $DIR/temp.tmp $DIR/supercell.tmp);
      qx(cd $DIR; supercell < $DIR/temp.tmp);
      qx(cp $DIR/${CASE}_super.struct $DIR/super.struct);
      qx(cp $DIR/${CASE}_super.struct $DIR/${CASE}.struct);
      &StructRead;
      $s_lattice="P";
      $s_spacegr="1_P1";
      #change all angles to 89°
      $s_aa=89;
      $s_bb=89;
      $s_cc=89;
      &StructWrite;
      qx(mv $DIR/$CASE.struct_i $DIR/$CASE.struct);
      qx(cp $DIR/${CASE}.struct $DIR/angles89.struct);
      qx(cd $DIR;x nn -d);
      qx(cd $DIR; /bin/echo 2 | nn nn.def);
      unless(open(FILE,">$DIR/temp.tmp"))
        {
        &CGIError("Context","could not create input file\n");
        exit;
        }
      print FILE "#!/bin/csh -f\nif ( `grep -ce \"[0-9]\" X3frozen.struct_nn` > 10 ) then\nmv $CASE.struct_nn $CASE.struct\nendif";
      close(FILE);
      qx(cd $DIR; chmod +x $DIR/temp.tmp; temp.tmp;);
      &StructRead;
      $s_lattice="P";
      $s_spacegr="1_P1";
      #change all angles back to 90°
      $s_aa=90;
      $s_bb=90;
      $s_cc=90;
      &StructWrite;
      qx(mv $DIR/$CASE.struct_i $DIR/$CASE.struct);
      qx(cp $DIR/${CASE}.struct $DIR/angles90.struct);
      qx(cp $DIR/$CASE.struct $DIR/${CASE}_brokensymm.struct);
      #now we should have a file with individual positions!
    }
    else
    {
      #in the other case we already have the correct atom-positions. Here we have to retrieve them initially
      &StructRead;
    }
    #now we open the file and extract all positions
    unless(open(STRUCTFILE,"$mode_file"))
      {
      &CGIError("Context","could not open $mode_file\n");
      exit;
      }
    @lines = <STRUCTFILE>;
    close(STRUCTFILE);
    #now we restore again the old file. the brokensymm contains the other one!
    qx(cp $DIR/${CASE}_initial.struct $DIR/${CASE}.struct);
    $fphon="#    Leave at least one blank between the numbers!\n";
    $fphon.="#    This content will be saved to ${CASE}.mode\n";
    $fphon.="#    all Text after the # is meant as comment\n";
    $fphon.="#    redundant lines will be ignored\n";
    $fphon.="#    if you load a mode file, its lines will be modified\n";
    $fphon.="#    to fit the current .struct file!\n\n";
    $fphon.="#        dx        dy        dz      # current positions x y z\n\n";
    $atom_mult=1;#the atom-multiplicity of the preceding lines. Always 1 as default!
    my $atomcount;
    $atomcount=1;
    my $i,$j,$temp;
    for($i=1;$i<=$s_nato;$i++)  
      {
        $s_name[$i] =~ s/\W//g;
        $s_nameadd[$i] =~ s/\W//g;
        if($break_symm == 1)
        {
	  for($j=1;$j<=$s_mult[$i];$j++)
	  {
	    $temp="$s_name[$i] $s_nameadd[$i] [$j]";
	    $fphon.="$temp  0.000000  0.000000  0.000000   #".sprintf("%1.6f",$s_x[$i][$j]);
	    $fphon.="  ".sprintf("%1.6f",$s_y[$i][$j])."  ".sprintf("%1.6f",$s_z[$i][$j])."\n";
          }
	}
	else
	{
	  $temp="$s_name[$i] $s_nameadd[$i]";
	  $fphon.="$temp  0.000000  0.000000  0.000000   #".sprintf("%1.6f",$s_x[$i][1]);
	  $fphon.="  ".sprintf("%1.6f",$s_y[$i][1])."  ".sprintf("%1.6f",$s_z[$i][1])."\n";
        }
      }
    }
    $INPUT1 =~ s/$DIR\///;
    $OUT .= <<__STOP__;
    <FORM ACTION=/exec/fphonons.pl METHOD=POST>
__STOP__
    &PassHiddenParms();
    $OUT .=  <<__STOP__;
    <INPUT TYPE=hidden NAME=nextinteractive VALUE=1>
    <INPUT TYPE=hidden NAME=next VALUE="$next">
    <INPUT TYPE=hidden NAME=break_symm VALUE="$break_symm">
    <INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl&doit=2">
    <INPUT TYPE=hidden NAME=doit VALUE="2">
    <h2>Frozen Phonon calculation</h2>
	choose a saved mode eigenvector or enter a new one.<br>The current input in the textarea will be saved to ${CASE}.mode !<br><br>
    <SELECT NAME=INPUT1><br>
    <OPTION value=0>new mode eigenvector</OPTION>$MODEFILES
    </SELECT>
    <INPUT TYPE=submit value=LOAD><br>
    </FORM>
__STOP__
    $OUT .= <<__STOP__;
    <FORM ACTION=/exec/fphonons.pl METHOD=POST>
__STOP__
    &PassHiddenParms();
    $fphon =~ s/ {3,}/   /ig;
    $disabled="disabled";
    $chkboxcolor="#AAAAAA";
    if($break_symm == 1)
    {
      $disabled="";
      $chkboxcolor="#000000";
    }
    if($INPUT1 ne '' && $INPUT1 !=0)
    {
      $cleanvar="$INPUT1";
      $cleanvar =~ s/\./_/;
    }
    else
    {
      $cleanvar="${CASE}_mode";
    }
    $OUT .=  <<__STOP__;
    <INPUT TYPE=hidden NAME=nextinteractive VALUE=1>
    <INPUT TYPE=hidden NAME=next VALUE="$next">
    <INPUT TYPE=hidden NAME=break_symm VALUE="$break_symm">
    <INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl&doit=4">
    <INPUT TYPE=hidden NAME=doit VALUE="4">
    <INPUT TYPE=hidden NAME=INPUT1 VALUE="$INPUT1">
    <TEXTAREA NAME=INPUT2 ROWS=12 cols=80 wrap=off>$fphon</TEXTAREA>
    <br>
    mode amplitudes (values on separate lines)<br>
    <TEXTAREA NAME=INPUT3 ROWS=5 cols=80></TEXTAREA><br>
    Type of displacement parameter <INPUT TYPE=radio value=2 name="displacement_format" checked> fractional <INPUT TYPE=radio value=1 name="displacement_format">bohr<br>
    <INPUT TYPE=checkbox value=1 name="show_xyz" checked> display complete output<br>
    <INPUT TYPE=checkbox value=1 name="run_sgroup" $disabled><font color="$chkboxcolor"> run sgroup on every new struct-file</font><br>
    <INPUT TYPE=checkbox value=1 name="del_files"> delete old ${CASE}_mode*.struct files<br>
    <INPUT TYPE=submit value="create the .struct files">
__STOP__
}

if($doit==3)
{
  $OUT .= <<__STOP__;
  <b>We are in plot mode</b><br><br>
  <FORM ACTION=/exec/fphonons.pl METHOD=POST>
__STOP__
  &PassHiddenParms();
  $OUT .= <<__STOP__;
  <INPUT TYPE=hidden NAME=doit VALUE="7">
  <INPUT TYPE=SUBMIT VALUE="plot">
  energy curve using
  <INPUT NAME="SCFFILES" VALUE="$SCFFILES">.scf files<br>
  </FORM>
__STOP__
}

if($doit==4)
{
  if($INPUT2 eq "" || $INPUT3 eq "")
  {
    $OUT.= <<__STOP__;
    "No valid input. Please specify a mode eigenvector and at least one amplitude!\n\n";
    <FORM ACTION=/exec/fphonons.pl METHOD=POST>
__STOP__
    &PassHiddenParms();
    $OUT .=  <<__STOP__;
    <INPUT TYPE=hidden NAME=nextinteractive VALUE=1>
    <INPUT TYPE=hidden NAME=next VALUE="$next">
    <INPUT TYPE=hidden NAME=doit VALUE="2">
    <INPUT TYPE=hidden NAME=break_symm VALUE="$break_symm">
    <INPUT TYPE=hidden NAME=INPUT4 VALUE="$INPUT4">
    <INPUT TYPE=hidden NAME="del_files" VALUE="$del_files">
    <INPUT TYPE=hidden NAME="show_xyz" VALUE="$show_xyz">
    <INPUT TYPE=submit value="back">
    </FORM>	
__STOP__
    #exit;
  }
  else
  {
    $structname="${CASE}_mode";
    #in any case we save the current mode 
    #the case.mode file is written to the disk for later use.
      unless(open(NEW,">$DIR/$CASE.mode"))
        {
        &CGIError("Couldn't write to $DIR/$CASE.mode\n");
        exit;
        }
      $FORM{'INPUT2'} =~ s/\r//gs;
      print NEW "$FORM{'INPUT2'}";
      close(NEW);
      $OUT.="$DIR/$CASE.mode created (You may now save this mode by another name)... <br><br>";
    if($break_symm == 1)
      {$OUT.="<b>The mode eigenvector will break symmetry</b><br><br>";}
    if($del_files eq '1' || $del_files == 1)
    { 
      if($structname ne "${CASE}")
      {
        qx(rm $DIR/$structname*.struct);
        $OUT.="Old $structname*.struct files deleted!<br><br>";
      }
      else
      {
	$OUT.="Error: Attempt to delete ${CASE}.struct!<br><br>";
      }
    }
    #if we broke the symmetry, there must exist a CASE_brokensymm.struct file. So we must overwrite the current
    #struct file and restore it afterwards!
    if($break_symm == 1)
    {
      if (-e "$DIR/${CASE}_brokensymm.struct")
        {qx(cp $DIR/${CASE}_brokensymm.struct $DIR/${CASE}.struct);}
      else
        {
        &CGIError("Couldn't find $DIR/${CASE}_brokensymm.struct\n");
        exit;
        }
    }
    #we need to know the structure first, so we can display the mode-changes!
    &StructRead;
    $s_spacegr=~ s/ /_/;
    foreach $v (@lattype)
    {
      if ( $s_spacegr =~ / */ && $s_lattice =~ /$v/ || $s_spacegr =~ /$v/)
	  {$s_lattice=$v;}
    }
    #like in optimizer we need the number of changes
    $INPUT3=qx(echo "$INPUT3"|wc -w);
    $FORM{'INPUT3'}=~ s/\n\n/\n/;
    $INPUT4=$FORM{'INPUT3'};
    @modelines=split('\n',$INPUT2);
    @amplitudes=split(' ',$INPUT4);
    $mode_count=1;
    $OUT.="<table>";
    $OUT.="<tr><td colspan=4><b>Eigenvector</b></td><td width=30>&nbsp;</td>";
    if($displacement_format==2)
    {
	$OUT.="<td colspan=3 nowrap><b>struct modifier (bohr)</b></td></tr>";
    }
    else
    {
	$OUT.="<td colspan=3 nowrap><b>struct modifier (bohr)<br><i>devided by cell-parameters $s_a/$s_b/$s_c</i></b></td></tr>";
    }
    my $i,$j,$temp2;
    $i=1;#init of atom counter
    $j=1;#init of multiplicity counter
    # $modes[0][0][0] always holds the total number of modes
    # $modes[$i][0][0] holds the number of multiplicities
    foreach $mline (@modelines)
    {
      if($mline =~/([a-z]{1,2} +[0-9]*) *\[?([0-9]*)\]? +(-?[0-9.]{1,}) +(-?[0-9.]{1,}) +(-?[0-9.]{1,})/i && $i<=$s_nato && $j<= $s_mult[$i])
      {
	$temp2="";
	if($break_symm == 1)
	{$temp2="[$j]";}
	#$s_name[$i] =~ s/\W//g;
	#$s_nameadd[$i] =~ s/\W//g;
	$OUT.="<tr><td>&nbsp;&nbsp;$s_name[$i] $s_nameadd[$i] $temp2</td>";
	$OUT.="<td>&nbsp;&nbsp;".sprintf("%1.8f",$3)."&nbsp;&nbsp;</td>";
	$OUT.="<td>&nbsp;&nbsp;".sprintf("%1.8f",$4)."&nbsp;&nbsp;</td>";
	$OUT.="<td>&nbsp;&nbsp;".sprintf("%1.8f",$5)."&nbsp;&nbsp;</td><td width=30>&nbsp;</td>";
	$modes[0][0][0]=$i; #the number of atoms
	$modes[$i][0][0]=$j; #the number of multiplicities for overflow-check!
        $modes[$i][$j][0]=$1; #the atom-name from the .mode file
        $modes[$i][$j][1]=$3; #x-part of the vector for atom $i $j
        $modes[$i][$j][2]=$4; #y -#-
        $modes[$i][$j][3]=$5; #z -#-
	if($displacement_format==2)
	{
	    $OUT.="<td>&nbsp;&nbsp;".sprintf("%1.8f",$3)."&nbsp;&nbsp;</td><td>&nbsp;&nbsp;".sprintf("%1.8f",$4);
	    $OUT.="&nbsp;&nbsp;</td><td>&nbsp;&nbsp;".sprintf("%1.8f",$5)."</td></tr>";
	}
	else
	{
	    $OUT.="<td>&nbsp;&nbsp;".sprintf("%1.8f",$3/$s_a)."&nbsp;&nbsp;</td><td>&nbsp;&nbsp;".sprintf("%1.8f",$4/$s_b);
	    $OUT.="&nbsp;&nbsp;</td><td>&nbsp;&nbsp;".sprintf("%1.8f",$5/$s_c)."</td></tr>";
	}
	if($break_symm == 1)
	{
	  $j++;
	  if($j>$s_mult[$i])
	  {
	    $j=1;
	    $i++;
	  }
	}
	else
	  {$i++;}
      }
    }
    $OUT.="</table><br>";
    unless(open(MODE,">$DIR/fphonons.job"))
    {
      &CGIError("Couldn't write to $DIR/fphonons.job\n");
      exit;
    }
    $FPHONONS.= <<__STOP__;
#!/bin/csh -f
    #   Modify this script according to your needs:
    #      Uncomment one of the lines ...
    #      Change run_lapw to runsp_lapw or use different convergence criterium
    #      Change save_lapw -d XXX

    #cd DIR #do not uncomment or delete this line!

    cp ${CASE}.struct ${CASE}_initial.struct #backup

__STOP__
    $filelist="";
    foreach $amp (@amplitudes)
    {
      &StructRead;
      if($break_symm == 1)
      {
        $s_lattice="1_P1";
        $s_spacegr="1_P1";
      }
      else
      {
        $s_spacegr=~ s/ /_/;
        foreach $v (@lattype)
        {
          if ( $s_spacegr =~ / */ && $s_lattice =~ /$v/ || $s_spacegr =~ /$v/)
	  {$s_lattice=$v;}
        }
      }
      if($show_xyz eq '1' || $show_xyz == 1)
      {
	$OUT.="<table>";
        $OUT.="<tr><td colspan=4><b>original struct</b></td><td width=30>&nbsp;</td>";
	$OUT.="<td colspan=3><b>modified by Amplitude:$amp</b></td></tr>";
      }
      my $i,$j,$my_mult,$temp,$tempname;
      for($i=1;$i<=$s_nato;$i++)
      {
	$my_mult=1;
	for($j=1;$j<=$my_mult;$j++)
	{
	  if($break_symm == 1)
	  {
	    $my_mult=$s_mult[$i];
	    $temp=" [$j]";
          }
	  if($modes[0][0][0]<$i || $modes[$i][0][0]<$j)
	  {
	    #prevent an error by array overflow!
            $s_name[$i] =~ s/\W//g;
            $s_nameadd[$i] =~ s/\W//g;
	    $tempname="$s_name[$i] $s_nameadd[$i]";
	    if($break_symm == 1)
              {$tempname.=" [$j]";}
	    $modes[$i][$j][0]=$tempname;
	    $modes[$i][$j][1]=0;
	    $modes[$i][$j][2]=0;
	    $modes[$i][$j][3]=0;
	    $temp="";#delete the $temp content, or you will see some strange results... this affects only the displayed atomname
	  }
	  if($show_xyz == 1)
	  {
            $OUT.="<tr><td>$modes[$i][$j][0]$temp</td><td>&nbsp;&nbsp;$s_x[$i][$j]&nbsp;&nbsp;</td>";
	    $OUT.="<td>&nbsp;&nbsp;$s_y[$i][$j]&nbsp;&nbsp;</td><td>&nbsp;&nbsp;$s_z[$i][$j]</td><td width=30>&nbsp;</td>";
	  }
	  if($displacement_format==2)
	  {
	      $s_x[$i][$j]=$s_x[$i][$j]+$modes[$i][$j][1]*$amp;
	      $s_y[$i][$j]=$s_y[$i][$j]+$modes[$i][$j][2]*$amp;
	      $s_z[$i][$j]=$s_z[$i][$j]+$modes[$i][$j][3]*$amp;
	  }
	  else
	  {
	      $s_x[$i][$j]=$s_x[$i][$j]+$modes[$i][$j][1]/$s_a*$amp;
	      $s_y[$i][$j]=$s_y[$i][$j]+$modes[$i][$j][2]/$s_b*$amp;
	      $s_z[$i][$j]=$s_z[$i][$j]+$modes[$i][$j][3]/$s_c*$amp;
	  }
	  if($show_xyz == 1)
	  {
            $OUT.="<td>".sprintf("%1.8f",$s_x[$i][$j])."&nbsp;&nbsp;</td>";
	    $OUT.="<td>&nbsp;&nbsp;".sprintf("%1.8f",$s_y[$i][$j])."&nbsp;&nbsp;</td>";
	    $OUT.="<td>&nbsp;&nbsp;".sprintf("%1.8f",$s_z[$i][$j])."</td></tr>";
	  }
        }
      }
      if($show_xyz == 1)
        {$OUT.="</table>";}
        my $x;
        $filler="";
        for($x=0;$x<6-length($amp);$x++)
        {
          $filler.="_";
        }
        $filelist.="$structname$filler$amp \\\n";
        &StructWrite;
        if($run_sgroup == 1)
	{
          qx(cp $DIR/${CASE}.struct $DIR/${CASE}_struct.bak);
          qx(cd $DIR; sgroup -wi -wo ${CASE}.struct_i ${CASE}.struct);
          qx(mv $DIR/${CASE}.struct_i $DIR/${CASE}.backupstruct_$amp;);
	  qx(cd $DIR; x symmetry -d; symmetry symmetry.def);
	  qx(mv $DIR/$CASE.struct_st $DIR/$structname$filler$amp.struct);
          qx(cp $DIR/${CASE}_struct.bak $DIR/${CASE}.struct);
          $OUT.="<br>executing sgroup, symmetry....";
	}
        else
	{qx(mv $DIR/$CASE.struct_i $DIR/$structname$filler$amp.struct);}
        $OUT.="<b>$structname$filler$amp.struct created</b><BR><br>";
      }
      #qx(mv $DIR/$CASE.struct_i $DIR/$structname$filler$amp.struct_ii);#the same as save file and clean up in structwrite.pl!
      $OUT.="<br>fphonons.job created<BR><br>";
      $FPHONONS.= <<__STOP__;
      foreach i (\\\n$filelist)
      cp  \$i.struct ${CASE}.struct
   #  cp  \$i.clmsum ${CASE}.clmsum
   #  x dstart
   #  run_lapw -ec 0.0001 -in1new 3 -in1orig -renorm
      run_lapw -ec 0.0001
      set stat = \$status
      if (\$stat) then
         echo ERROR status in \$i
         exit 1
      endif
      save_lapw  \$i
   #  save_lapw  -f -d XXX \$i
   end
      cp ${CASE}_initial.struct ${CASE}.struct #restore
__STOP__
    print  MODE $FPHONONS;
    close(MODE);
    qx(cp $DIR/${CASE}_initial.struct $DIR/${CASE}.struct);
    $OUT.="<br>original ${CASE}.struct restored<BR>";
    qx(chmod +x "$DIR/fphonons.job");
    $OUT.= <<__STOP__;
    <FORM ACTION=/exec/fphonons.pl METHOD=POST>
__STOP__
    &PassHiddenParms();
    $OUT .=  <<__STOP__;
    <INPUT TYPE=hidden NAME=doit VALUE="1">
    <INPUT TYPE=submit value="$next">
    </FORM>
__STOP__
  }
}

if($doit==5)
{
  $OUT.= <<__STOP__;
  <FORM ACTION=/exec/fphonons.pl METHOD=POST>
__STOP__
  &PassHiddenParms();
  $OUT .=  <<__STOP__;
  <b>Save current data</b><br><br>
  <i>please notice: This function will save the current ${CASE}.mode file<br> by the new name. All struct files may be renamed accordingly (checkbox)<br>as well as the entries in fphonons.job will be adapted<br>optionally you may create a new case (recommended)</i><br><br>
  <INPUT TYPE=hidden NAME=doit VALUE="6">
  <INPUT TYPE=hidden NAME=next VALUE="$next">
  <INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl&doit=2">
  <INPUT TYPE=text NAME=mode_name><br>
  <INPUT type=checkbox name=chfilename value=1 checked> rename <b>all</b> current <i>${CASE}_mode.struct</i> files by this name and update fphonons.job<br>
  <INPUT type=checkbox name=create_session value=1 checked>move all current ${CASE}_mode.struct files into a new session<br>
  <INPUT TYPE=submit value="Save">
  <input type=hidden value=0 name=movetosession>
  </FORM>
__STOP__
}

if($doit==6)
{
  my $error;
  $error=0;
  if(chomp($FORM{'mode_name'}) eq '' || $mode_name eq '')
  {
    $nextdoit=5;
    $OUT .= "No valid name specified!";
  }
  else
  {
    $nextdoit=1;
    &listdir($DIR);#we need that for listing the directories as well as for the single files!
    if($create_session == 1 && $movetosession == 0)
    {
      $test=0;
      foreach $i (@subdirs)
        {$test=1 if ($i eq "${CASE}_$FORM{'mode_name'}");}
      if($test)
      {
        $error=1;
        $OUT.= <<__STOP__;
	The name you specified does already exist as a session<br>
	Do you want to move the current files into it?<br>
	<FORM ACTION=/exec/fphonons.pl METHOD=POST>
__STOP__
        &PassHiddenParms();
	$OUT .=  <<__STOP__;
	<input type=hidden name=create_session value=1>
	<input type=hidden name=mode_name value=$mode_name>
	<input type=hidden name=doit value=6>
	<INPUT TYPE=hidden NAME=next VALUE="$next">
	<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
	<INPUT TYPE=hidden NAME=movetosession value=1>
	<INPUT TYPE=hidden NAME=create_session value=1>
	<INPUT TYPE=hidden NAME=chfilename value=$chfilename>
        <input type=submit value="YES">
	</FORM>
	<FORM ACTION=/exec/fphonons.pl METHOD=POST>
__STOP__
        &PassHiddenParms();
	$OUT .=  <<__STOP__;
	<input type=hidden name=doit value=5>
	<INPUT TYPE=hidden NAME=next VALUE="$next">
	<INPUT TYPE=hidden NAME=nexturl VALUE="$nexturl">
        <input type=submit value="NO">
	</FORM>
__STOP__
      }
      else
      {
        qx(mkdir  "$DIR/${CASE}_$FORM{'mode_name'}"); # create a new directory
        $randmax=999999;
	$newsid = int(rand($randmax));
	while ( -e "$sessionpath/$newsid" )
          {$newsid = int(rand($randmax));}
        #save the old Values!
        $NAME_temp=$NAME;
        $SID_temp=$SID;
        $DIR_temp=$DIR;
        $COMMENT_temp=$COMMENT;
	unless(open(FILE,">$DIR/__NEWID__"))
	{
	&CGIError("Can't write file $fname.\n");
	exit;
        }
	print FILE $newsid;
	close(FILE);
        #overwrite!
        $NAME="${CASE}_$FORM{'mode_name'}";
        $SID=$newsid;
        $DIR="$DIR/${CASE}_$FORM{'mode_name'}";
        $COMMENT="This session is a fphonon session of ${CASE} created by fphonons.pl!";
        #restore (if anything went wrong, we face a problem, because we work with the wrong SID!)
	&SaveSession;
        $NAME=$NAME_temp;
        $SID=$SID_temp;
        $DIR=$DIR_temp;
        $COMMENT=$COMMENT_temp;
      }
    }
    if($chfilename == 1 && !$error)
    {
      foreach $line(@ascii_files)
      {
        if($line =~/${CASE}_mode([^ ]+)\.struct/i)
        { 
          qx(mv "$DIR/${CASE}_mode$1.struct" "$DIR/$FORM{'mode_name'}_mode$1.struct");
	  $OUT.= "<i>${CASE}_mode$1.struct</i> renamed to <i>$FORM{'mode_name'}_mode$1.struct</i><br>";
	  $found="$FORM{'mode_name'}_mode$1.struct";#at least one file found and stored for later use!
	  if($create_session == 1)
          {
	    qx(mv "$DIR/$FORM{'mode_name'}_mode$1.struct" "$DIR/${CASE}_$FORM{'mode_name'}");
	    $OUT.="$FORM{'mode_name'}_mode$1.struct moved to /${CASE}_$FORM{'mode_name'}<br>";
	  }
        }  
      }  
      if($found)
      {
        qx(sed "s/${CASE}_mode/$FORM{'mode_name'}_mode/g" $DIR/fphonons.job > $DIR/_tmp_ );
	if($create_session == 1)
        {
	  qx(cd "$DIR/${CASE}_$FORM{'mode_name'}"; cp $found "${CASE}_$FORM{'mode_name'}.struct");#we create an initial struct file
	  qx(sed "s/${CASE}/${CASE}_$FORM{'mode_name'}/g" "$DIR/_tmp_"  > "$DIR/__tmp__" );
	  qx(sed "s/#cd DIR/cd ${CASE}_$FORM{'mode_name'}    /" "$DIR/__tmp__" > $DIR/fphonons.job);
	  qx(rm "$DIR/*_tmp_*");
	  $OUT.="fphonons.job updated with session-directory!<br>";
	}
	else
          {qx(mv "$DIR/_tmp_" $DIR/fphonons.job);}
        qx(chmod +x "$DIR/fphonons.job");
	$OUT.= "<br>entries in fphonons.job changed to $FORM{'mode_name'}_mode*.struct<br><br>";
	qx(cp "$DIR/${CASE}.mode" "$DIR/$FORM{'mode_name'}.mode");
	$OUT.= "<i>${CASE}.mode</i> copied to <i>$FORM{'mode_name'}.mode</i>";
      }
      else
      {
	$OUT.="<br>no files found!</br>, no changes made<br>";
      }
    }
  }
  if(!$error)
  {
    $OUT.= <<__STOP__;
    <FORM ACTION=/exec/fphonons.pl METHOD=POST>
__STOP__
    &PassHiddenParms();
    $OUT .=  <<__STOP__;
    <INPUT TYPE=hidden NAME=doit VALUE="$nextdoit">
    <INPUT TYPE=hidden NAME=next VALUE="$next">
    <INPUT TYPE=submit value="$next">
    </FORM>
__STOP__
  }
}

if($doit == 7)
{
  @files = sort (glob ("$DIR/${SCFFILES}.scf"));
  if (-e "$DIR/$CASE.vol")
  {
    $umps = qx( rm $DIR/$CASE.vol);
  }
  $OUT.="$type";
  $type="test.vol";#just to override some strange construction
  if ($type =~ /vol/)
  {
    unless(open(FILE,">$DIR/$CASE.vol"))
    {
      &CGIError("Can't write file $fname.\n");
      exit;
    }
    my $x;
    foreach $i (@files)
    {
      $line=$i;#save $i as a copy
      $line =~/${CASE}.*_mode_*([^ ]+)\.scf/i;
      #$i= qx(grep :VOL $i | tail -1 | cut -f2 -d= );
      #chomp($x);
    #$x=$1;#to set the x-axis as the Amplitude!
      $y = qx(grep :ENE $i | tail -1 | cut -f2 -d= );
      chomp($y);
      $OUT.=$y;
      print FILE "$1 $y\n";
    }
    close(FILE);
    $umps = qx( cd $DIR && echo "" | x eosfit );
    $tmp = qx(grep V0 $DIR/$CASE.outputeos | tail -1) ;
    $murna = split(" ",$tmp);
    $plotfile="$tempdir/$SID-$$";
    unless(open(FILE,">$DIR/:eplot"))
    {
      &CGIError("Can't write file $fname.\n");
      exit;
    }
    $t1=$murna[0];
    shift(@murna);

    print FILE <<__STOP__;
set terminal png
set output '$plotfile.png'
set format y "%.4f"
set title "$CASE"
set xlabel "Amplitude"
set ylabel "Energy [Ry]"
plot "$CASE.vol" title "Murnaghan: $t1" w p, "$CASE.eosfit" title "@murna" w l 
set terminal postscript
set output '$plotfile.ps'
replot
__STOP__

    close(FILE);
    $umps = qx(cd $DIR;gnuplot ":eplot" 2>&1);
    $OUT .= "<br><IMG SRC=/tmp/$SID-$$.png><br clear=all><br>";
    $OUT .= "<A HREF=/tmp/$SID-$$.ps>Download hardcopy in PostScript format</A>";
    $OUT .= "<br><pre>";
    $OUT .= qx(cat $DIR/$CASE.outputeos);
  }
}

PrintPage("Context",$OUT);
