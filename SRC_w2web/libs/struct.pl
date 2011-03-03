sub StructRead {
	my $i;
	my $j;
	$structedit=0;

	$structfile = "$DIR/$CASE.struct";
	$testfile = "$DIR/$CASE.struct_i";

	if ( -e $structfile ) {
		# good, we have a struct-file
	} else {
		# we need a template
		redirectURL("/util/structask.pl?SID=$SID");
		exit;
	}

	if ( -e $testfile ) {
		$structedit=1;
		$structfile .= "_i"
	}

	open (STRUCT, $structfile);
	# read title
	$_ = <STRUCT>;
	chomp;
	$s_title = $_;
#	$s_title =~ s/^ *([^ ]*) */\1/;
	$s_title =~ s/^ *//;
	$s_title =~ s/ *$//;

	# read lattice, # ineq. atoms and spacegroup
	$_ = <STRUCT>;
	chomp;
	$s_lattice = substr $_, 0,4;
	$s_ineq = substr $_,27,3;
	$s_ineq =~ s/^ *//;
	$s_spacegr = substr $_,30,29;
	$s_spacegr =~ s/^ *//;

	$generateeq = 0;
	$generateeq = 1 if (length $s_spacegr > 0);
	$s_nato = $s_ineq;

	#read relativistic
	$_ = <STRUCT>;
	chomp;
	$s_rels = substr $_,13,4;
  $unit = substr $_,23,4;
	if  ( $unit =~ /ang/ ) {} else { $unit ="bohr"; }

	#lattice parameters
	  #init angles
	$s_aa = $s_bb = $s_cc = 90;
	if ("$s_lattice" =~ m/H  /) {
		$s_cc = 120;
	}

	$_ = <STRUCT>;
	chomp;
	$s_a  = substr $_,0,10;
	$s_a  =~ s/^ *//;
	$s_b  = substr $_,10,10;
	$s_b  =~ s/^ *//;
	$s_c  = substr $_,20,10;
	$s_c  =~ s/^ *//;
	$s_aa = substr $_,30,10;
	$s_aa =~ s/^ *//;
	$s_bb = substr $_,40,10;
	$s_bb =~ s/^ *//;
	$s_cc = substr $_,50,10;
	$s_cc =~ s/^ *//;

	if ($unit =~ /ang/) {
		$s_a *= 0.529177;
		$s_b *= 0.529177;
		$s_c *= 0.529177;
	}
	for ($i = 1; $i <= $s_nato; $i++) {
		$_ = <STRUCT>;
		chomp;
		$at = substr $_,4,4;
		$atom = $at;
		# check if cubic
		$s_cubic[$i] = 0;
		$s_cubic[$i] = 1 if ($atom <0) ;

		$x = substr($_,12,11);
		$x =~ s/^ *//;
		$y = substr($_,25,11);
		$y =~ s/^ *//;
		$z = substr($_,38,11);
		$z =~ s/^ *//;

		$s_x[$i][1] = $x;
		$s_y[$i][1] = $y;
		$s_z[$i][1] = $z;

		#get multiplicity and isplit
		$isplit = 8;
		$_ = <STRUCT>;
		chomp;
		$multi = substr($_,15,2);
		$isplit = substr($_,35,2);

		$s_mult[$i] = $multi;
		$s_isplit[$i] = $isplit;

		if ($multi > 1) {
			#allow for multiplicity!
			for ($j = 2; $j <= $multi; $j++) {
				$_ = <STRUCT>;
				chomp;
				$atom = substr($_,5,4);
				$x = substr($_,12,11);
				$x =~ s/^ *//;
				$y = substr($_,25,11);
				$y =~ s/^ *//;
				$z = substr($_,38,11);
				$z =~ s/^ *//;

				$s_x[$i][$j] = $x;
				$s_y[$i][$j] = $y;
				$s_z[$i][$j] = $z;
			}
		}


		$nameadd = "";
		# get name and params
		$_ = <STRUCT>;
		chomp;
		$s_name[$i] = substr($_,0,2);
		$s_name[$i] =~ s/^ *//;
		$s_nameadd[$i] = substr($_,2,8);
		$s_nameadd[$i] =~ s/^ *//;
		$s_npt[$i] = substr($_,15,6);
		$s_npt[$i] =~ s/^ *//;
		$s_ro[$i] = substr($_,25,10);
		$s_rmt[$i] = substr($_,41,10);
		$s_rmt[$i] =~ s/^ *//;
		$s_zz[$i] = substr($_,56,5);
		$s_zz[$i] =~ s/^ *//;
				
		#read rotloc
		$_ = <STRUCT>;
		chomp;
		$s_rot11[$i] = substr($_,20,10);
		$s_rot12[$i] = substr($_,30,10);
		$s_rot13[$i] = substr($_,40,10);
		$_ = <STRUCT>;
		chomp;
		$s_rot21[$i] = substr($_,20,10);
		$s_rot22[$i] = substr($_,30,10);
		$s_rot23[$i] = substr($_,40,10);
		$_ = <STRUCT>;
		chomp;
		$s_rot31[$i] = substr($_,20,10);
		$s_rot32[$i] = substr($_,30,10);
		$s_rot33[$i] = substr($_,40,10);
	}

	# read nsym
	$_ = <STRUCT>;
	chomp;
	$s_nsym = substr($_,0,4);

	for ($i = 1; $i <= $s_nsym; $i++) {
		$_ = <STRUCT>;
		chomp;
		$s_s11[$i] = substr($_,0,2);
		$s_s12[$i] = substr($_,2,2);
		$s_s13[$i] = substr($_,4,2);
		$s_t1[$i] = substr($_,6,11);
		$_ = <STRUCT>;
		chomp;
		$s_s21[$i] = substr($_,0,2);
		$s_s22[$i] = substr($_,2,2);
		$s_s23[$i] = substr($_,4,2);
		$s_t2[$i] = substr($_,6,11);
		$_ = <STRUCT>;
		chomp;
		$s_s31[$i] = substr($_,0,2);
		$s_s32[$i] = substr($_,2,2);
		$s_s33[$i] = substr($_,4,2);
		$s_t3[$i] = substr($_,6,11);
		
		$_ = <STRUCT>;
		chomp;

	}
	close STRUCT;
}
# end StructRead

sub StructWrite {
	&SG2lat;
	$s_spacegr =~ s/ //;
	open (STRUCT, ">$DIR/$CASE.struct_i");
	print STRUCT sprintf "%-60s\n", $s_title;
	print STRUCT sprintf "%-4s%23s%3i%-30s\n", $s_lattice,"LATTICE,NONEQUIV.ATOMS:", $s_nato, $s_spacegr;
	print STRUCT sprintf "%-13s%-4s%-6s%-4s\n", "MODE OF CALC=", $s_rels," unit=",$unit;
	if ($unit =~ /ang/) {
		$s_a /= 0.529177;
		$s_b /= 0.529177;
		$s_c /= 0.529177;
	}
	print STRUCT sprintf "%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n", $s_a, $s_b, $s_c,$s_aa,$s_bb,$s_cc;
	for ($i = 1; $i <= $s_nato; $i++)
	{
	    $s_x[$i][1]= sprintf("%10.8f",&Calculette($s_x[$i][1], $i));
	    $s_y[$i][1]= sprintf("%10.8f",&Calculette($s_y[$i][1], $i));
	    $s_z[$i][1]= sprintf("%10.8f",&Calculette($s_z[$i][1], $i));
		if ($generateeq) {
			$sg = $s_spacegr;
			$sg =~ s/[0-9]*_//;
			if ( $s_lattice =~ /R/ ) {
				$OUT .= "<b>RHOMBOHEDRAL!!!!!</B>" if $debug;
				#spacegroup need hex-input!
				$in ="$s_x[$i][1] $s_y[$i][1] $s_z[$i][1]\n";
				$out = qx(/bin/echo "$in" | rhomb2hex | head -5 | tail -1);
				$OUT .= " xx $out  xx<br>" if $debug;
				$out .= " ";
				$out =~ s/  / /g;
				@line = split (/ /, $out);
				$count=0;
				foreach $j (@line) { 
					$hex[$count] = $j;
					$count++;
				}
				$OUT .= "$hex[1] --- $hex[2] --- $hex[3]<br>" if $debug;
				$s_x[$i][1]= $hex[1];
				$s_y[$i][1]= $hex[2];
				$s_z[$i][1]= $hex[3];
			}
			$in = "$s_a\n $s_b\n $s_c\n $s_aa\n $s_bb\n $s_cc\n $sg \n XXX\n $s_x[$i][1] \n $s_y[$i][1] \n $s_z[$i][1] \n";
                        $in =~ s/( [0-9]+)(e-[0-9]+)/$1.$2/ig; # this will change all xe-y to the fortran compliant format x.e-y
			$out = qx(/bin/echo "$in" | spacegroup);
			$OUT .= "$in \n $out" if $debug;
			@lines = split ('\n', $out);
			$count=0;
			$error=0;
			foreach $line (@lines) {
				if ($line =~ /Error - .*/ && $error==0) {
					$error=1;
					$tmp = $line;
					$tmp =~ s/^.* ://g;
					$tmp =~ s/Error/<b>Error<\/b>/g;
					$OUT .= "$tmp<br>";
				}
				$line =~ s/ Coordinates of the 1st atom : //;
				if ($line =~ / Xxx/) {
					$count++;
					$s_x[$i][$count] = substr($line,14,11);
					$s_y[$i][$count] = substr($line,26,11);
					$s_z[$i][$count] = substr($line,38,11);
					if ( $s_lattice =~ /R/ ) {
						#convert back to rhomb need hex-input!
						$in ="$s_x[$i][$count] $s_y[$i][$count] $s_z[$i][$count]\n";
						$out = qx(/bin/echo "$in" | hex2rhomb | head -5 | tail -1);
						$out .= " ";
						$out =~ s/  / /g;
						@line = split (/ /, $out);
						$c=0;
						foreach $item (@line) { 
							$rhomb[$c] = $item;
							$c++;
						}
						$OUT .= "$rhomb[1] --- $rhomb[2] --- $rhomb[3]<br>" if $debug;
						$s_x[$i][$count]= $rhomb[1];
						$s_y[$i][$count]= $rhomb[2];
						$s_z[$i][$count]= $rhomb[3];
					}
					$s_mult[$i]=$count;
				}
			}
		}

		$atomnum = $i;
		$atomnum = -1*$i if ($s_cubic[$i]);
		print STRUCT sprintf "%-4s%4i%-4s%10.8f%-3s%10.8f%-3s%10.8f\n","ATOM",$atomnum,": X=",$s_x[$i][1]," Y=",$s_y[$i][1]," Z=",$s_z[$i][1];
	
		print STRUCT sprintf "%-15s%2i%17s%2i\n", "          MULT=",$s_mult[$i],"ISPLIT=",$s_isplit[$i];

		if ($s_mult[$i] > 1) {
			#allow for multiplicity!
			for ($j = 2; $j <= $s_mult[$i]; $j++) {
				$s_x[$i][$j]= &Calculette($s_x[$i][$j], $i);
				$s_y[$i][$j]= &Calculette($s_y[$i][$j], $i);
				$s_z[$i][$j]= &Calculette($s_z[$i][$j], $i);
				print STRUCT sprintf "%-4s%4i%-4s%10.8f%-3s%10.8f%-3s%10.8f\n","ATOM",$atomnum,":X=",$s_x[$i][$j]," Y=",$s_y[$i][$j]," Z=",$s_z[$i][$j];
			}
		}
		$elem = $s_name[$i];
		$z = $s_zz[$i];
		$ro = $s_ro[$i];
		&Element() if ($s_zz[$i]==0);
		$s_name[$i] =$elem;
		$s_zz[$i] = $z;
		$s_ro[$i] = $ro;
		print STRUCT sprintf "%-2s%-8s%-5s%5i%-5s%10.8f%-5s%10.4f%-5s%5.1f\n", $s_name[$i],$s_nameadd[$i]," NPT=",$s_npt[$i],"  R0=",$s_ro[$i]," RMT=",$s_rmt[$i],"   Z:",$s_zz[$i];
		print STRUCT sprintf "%-20s%10.7f%10.7f%10.7f\n", "LOCAL ROT MATRIX:   ",$s_rot11[$i],$s_rot12[$i],$s_rot13[$i];
		print STRUCT sprintf "%20s%10.7f%10.7f%10.7f\n", "                    ",$s_rot21[$i],$s_rot22[$i],$s_rot23[$i];
		print STRUCT sprintf "%20s%10.7f%10.7f%10.7f\n", "                    ",$s_rot31[$i],$s_rot32[$i],$s_rot33[$i];

	}
	print STRUCT sprintf "%4i%-s\n",$s_nsym, "      NUMBER OF SYMMETRY OPERATIONS";
	for ($i = 1; $i <= $s_nsym; $i++) {
		print STRUCT sprintf "%2i%2i%2i%11.8f\n", $s_s11[$i],$s_s12[$i],$s_s13[$i],$s_t1[$i];
		print STRUCT sprintf "%2i%2i%2i%11.8f\n", $s_s21[$i],$s_s22[$i],$s_s23[$i],$s_t2[$i];
		print STRUCT sprintf "%2i%2i%2i%11.8f\n", $s_s31[$i],$s_s32[$i],$s_s33[$i],$s_t3[$i];
	print STRUCT sprintf "%8i\n", $i;
	}
	

	$OUT .= "struct_i written\n" if $debug;
}

sub Calculette {
	my ($myval, $myatom) = @_;
##	qx(echo "init:$myval/$myatom" >> $DIR/mistoutput);
	#	the genius calculette for doing wizard like
	#	input of positions like 1/2+1/6 !!!!
	#
	#	Taken from the original wienbox-calculette -
	#  this piece of code was written without the aid of grappa, honestly !

	#  (That might explain the bugs :-)
	
	#save the walrus
	$myval =~ s/e\-/EMINUS/;
	$myval =~ s/x/$s_x[$myatom][1]/;
	$myval =~ s/y/$s_y[$myatom][1]/;
	$myval =~ s/z/$s_z[$myatom][1]/;
	$myval =~ s/\+/*1.0+1.0*/;
	$myval =~ s/\-/+0.0-/;
	$myval =~ s/\-/*1.0-1.0*/;
	#restore the walrus
	$myval =~ s/EMINUS/e-/;
##	qx(echo "after regexp:$myval/$myatom" >> $DIR/mistoutput);
	$calcinput = "\$result = 1.0 * $myval * 1.0";
##	qx(echo "before eval:$calcinput" >> $DIR/mistoutput);
	eval $calcinput;
##	qx(echo "result:$result" >> $DIR/mistoutput);
	if ($result >=0) {
		$result=$result - int($result);
	} else {
		$result=$result +1 - int($result);
	}
##	qx(echo "result_afterwards:$result" >> $DIR/mistoutput);
	$result=0 if ($result==1) ;
	$OUT .= "$myval<br>$calcinput<br>$result<br>" if $debug;
	return $result;
}
sub InstGen {

# inst file generator

	$instfile = "$DIR/$CASE.inst";
	open(FILE,">$instfile");


	for ($i=1; $i <= $s_nato; $i++) {
		$name = $s_name[$i];
		$add = $s_nameadd[$i];
		if ($name =~ /H /)  {	
			print FILE "H $add\n";
			print FILE "1  \n";
			print FILE "1,-1,0.9  N\n";
			print FILE "1,-1,0.1  N\n";
		} elsif ($name =~ /D /) {
			print FILE "D $add\n";
			print FILE "1  \n";
			print FILE "1,-1,0.9  N\n";
			print FILE "1,-1,0.1  N\n";
		} elsif ($name =~ /He/) {
			print FILE "He $add\n";
			print FILE "1  \n";
			print FILE "1,-1,1.0  N\n";
			print FILE "1,-1,1.0  N\n";
		} elsif ($name =~ /Li/) {	
			print FILE "Li $add\n";
			print FILE "He 1  \n";
			print FILE "2,-1,1.0  N\n";
			print FILE "2,-1,0.0  N\n";
		} elsif ($name =~ /Be/) {	
			print FILE "Be $add\n";
			print FILE "He 1  \n";
			print FILE "2,-1,1.0  N\n";
			print FILE "2,-1,1.0  N\n";
		} elsif ($name =~ /B /) {	
			print FILE "B $add\n";
			print FILE "He 2  \n";
			print FILE "2,-1,1.0  N\n";
			print FILE "2,-1,1.0  N\n";
			print FILE "2, 1,1.0  N\n";
			print FILE "2, 1,0.0  N\n";
		} elsif ($name =~ /C /) {	
			print FILE "C $add\n";
			print FILE "He 3  \n";
			print FILE "2,-1,1.0  N\n";
			print FILE "2,-1,1.0  N\n";
			print FILE "2, 1,1.0  N\n";
			print FILE "2, 1,0.0  N\n";
			print FILE "2,-2,1.0  N\n";
			print FILE "2,-2,0.0  N\n";
		} elsif ($name =~ /N /) {	
			print FILE "N $add\n";
			print FILE "He 3  \n";
			print FILE "2,-1,1.0  N\n";
			print FILE "2,-1,1.0  N\n";
			print FILE "2, 1,1.0  N\n";
			print FILE "2, 1,0.0  N\n";
			print FILE "2,-2,2.0  N\n";
			print FILE "2,-2,0.0  N\n";
		} elsif ($name =~ /O /) {	
			print FILE "O $add\n";
			print FILE "He 3  \n";
			print FILE "2,-1,1.0  N\n";
			print FILE "2,-1,1.0  N\n";
			print FILE "2, 1,1.0  N\n";
			print FILE "2, 1,1.0  N\n";
			print FILE "2,-2,2.0  N\n";
			print FILE "2,-2,0.0  N\n";
		} elsif ($name =~ /F /) {	
			print FILE "F $add\n";
			print FILE "He 3  \n";
			print FILE "2,-1,1.0  N\n";
			print FILE "2,-1,1.0  N\n";
			print FILE "2, 1,1.0  N\n";
			print FILE "2, 1,1.0  N\n";
			print FILE "2,-2,2.0  N\n";
			print FILE "2,-2,1.0  N\n";
		} elsif ($name =~ /Ne/)  {	
			print FILE "Ne $add\n";
			print FILE "He 3  \n";
			print FILE "2,-1,1.0  N\n";
			print FILE "2,-1,1.0  N\n";
			print FILE "2, 1,1.0  N\n";
			print FILE "2, 1,1.0  N\n";
			print FILE "2,-2,2.0  N\n";
			print FILE "2,-2,2.0  N\n";
		} elsif ($name =~ /Na/) {	
			print FILE "Na $add\n";
			print FILE "Ne 1  \n";
			print FILE "3,-1,1.0  N\n";
			print FILE "3,-1,0.0  N\n";
		} elsif ($name =~ /Mg/) {	
			print FILE "Mg $add\n";
			print FILE "Ne 1  \n";
			print FILE "3,-1,1.0  N\n";
			print FILE "3,-1,1.0  N\n";
		} elsif ($name =~ /Al/) {	
			print FILE "Al $add\n";
			print FILE "Ne 2  \n";
			print FILE "3,-1,1.0  N\n";
			print FILE "3,-1,1.0  N\n";
			print FILE "3, 1,1.0  N\n";
			print FILE "3, 1,0.0  N\n";
		} elsif ($name =~ /Si/) {	
			print FILE "Si $add\n";
			print FILE "Ne 3  \n";
			print FILE "3,-1,1.0  N\n";
			print FILE "3,-1,1.0  N\n";
			print FILE "3, 1,1.0  N\n";
			print FILE "3, 1,0.0  N\n";
			print FILE "3,-2,1.0  N\n";
			print FILE "3,-2,0.0  N\n";
		} elsif ($name =~ /P /) {	
			print FILE "P $add\n";
			print FILE "Ne 3  \n";
			print FILE "3,-1,1.0  N\n";
			print FILE "3,-1,1.0  N\n";
			print FILE "3, 1,1.0  N\n";
			print FILE "3, 1,0.0  N\n";
			print FILE "3,-2,2.0  N\n";
			print FILE "3,-2,0.0  N\n";
		} elsif ($name =~ /S /) {	
			print FILE "S $add\n";
			print FILE "Ne 3  \n";
			print FILE "3,-1,1.0  N\n";
			print FILE "3,-1,1.0  N\n";
			print FILE "3, 1,1.0  N\n";
			print FILE "3, 1,1.0  N\n";
			print FILE "3,-2,2.0  N\n";
			print FILE "3,-2,0.0  N\n";
		} elsif ($name =~ /Cl/) {	
			print FILE "Cl $add\n";
			print FILE "Ne 3  \n";
			print FILE "3,-1,1.0  N\n";
			print FILE "3,-1,1.0  N\n";
			print FILE "3, 1,1.0  N\n";
			print FILE "3, 1,1.0  N\n";
			print FILE "3,-2,2.0  N\n";
			print FILE "3,-2,1.0  N\n";
		} elsif ($name =~ /Ar/)  {	
			print FILE "Ar $add\n";
			print FILE "Ne 3  \n";
			print FILE "3,-1,1.0  N\n";
			print FILE "3,-1,1.0  N\n";
			print FILE "3, 1,1.0  N\n";
			print FILE "3, 1,1.0  N\n";
			print FILE "3,-2,2.0  N\n";
			print FILE "3,-2,2.0  N\n";
		} elsif ($name =~ /K /)  {	
			print FILE "K $add\n";
			print FILE "Ar 1  \n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4,-1,0.0  N\n";
		} elsif ($name =~ /Ca/) {	
			print FILE "Ca $add\n";
			print FILE "Ar 1  \n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4,-1,1.0  N\n";
		} elsif ($name =~ /Sc/) {	
			print FILE "Sc $add\n";
			print FILE "Ar 2  \n";
			print FILE "3, 2,1.0  N\n";
			print FILE "3, 2,0.0  N\n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4,-1,1.0  N\n";
		} elsif ($name =~ /Ti/) {	
			print FILE "Ti $add\n";
			print FILE "Ar 2  \n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3, 2,0.0  N\n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4,-1,1.0  N\n";
		} elsif ($name =~ /V /) {	
			print FILE "V $add\n";
			print FILE "Ar 2  \n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3, 2,1.0  N\n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4,-1,1.0  N\n";
		} elsif ($name =~ /Cr/) {	
			print FILE "Cr $add\n";
			print FILE "Ar 3  \n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3,-3,1.0  N\n";
			print FILE "3,-3,0.0  N\n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4,-1,0.0  N\n";
		} elsif ($name =~ /Mn/) {	
			print FILE "Mn $add\n";
			print FILE "Ar 3  \n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3, 2,0.0  N\n";
			print FILE "3,-3,3.0  N\n";
			print FILE "3,-3,0.0  N\n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4,-1,1.0  N\n";
		} elsif ($name =~ /Fe/) {	
			print FILE "Fe $add\n";
			print FILE "Ar 3  \n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3,-3,2.5  N\n";
			print FILE "3,-3,0.0  N\n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4,-1,0.5  N\n";
		} elsif ($name =~ /Co/) {	
			print FILE "Co $add\n";
			print FILE "Ar 3  \n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3,-3,3.0  N\n";
			print FILE "3,-3,0.0  N\n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4,-1,1.0  N\n";
		} elsif ($name =~ /Ni/) {	
			print FILE "Ni $add\n";
			print FILE "Ar 3  \n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3,-3,3.0  N\n";
			print FILE "3,-3,1.0  N\n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4,-1,1.0  N\n";
		} elsif ($name =~ /Cu/) {	
			print FILE "Cu $add\n";
			print FILE "Ar 3  \n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3,-3,3.0  N\n";
			print FILE "3,-3,2.0  N\n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4,-1,1.0  N\n";
		} elsif ($name =~ /Zn/) {	
			print FILE "Zn $add\n";
			print FILE "Ar 3  \n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3,-3,3.0  N\n";
			print FILE "3,-3,3.0  N\n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4,-1,1.0  N\n";
		} elsif ($name =~ /Ga/) {	
			print FILE "Ga $add\n";
			print FILE "Ar 4  \n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3,-3,3.0  N\n";
			print FILE "3,-3,3.0  N\n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4, 1,1.0  N\n";
			print FILE "4, 1,0.0  N\n";
		} elsif ($name =~ /Ge/) {	
			print FILE "Ge $add\n";
			print FILE "Ar 5  \n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3,-3,3.0  N\n";
			print FILE "3,-3,3.0  N\n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4, 1,1.0  N\n";
			print FILE "4, 1,0.0  N\n";
			print FILE "4,-2,1.0  N\n";
			print FILE "4,-2,0.0  N\n";
		} elsif ($name =~ /As/) {	
			print FILE "As $add\n";
			print FILE "Ar 5  \n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3,-3,3.0  N\n";
			print FILE "3,-3,3.0  N\n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4, 1,1.0  N\n";
			print FILE "4, 1,0.0  N\n";
			print FILE "4,-2,2.0  N\n";
			print FILE "4,-2,0.0  N\n";
		} elsif ($name =~ /Se/) {	
			print FILE "Se $add\n";
			print FILE "Ar 5  \n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3,-3,3.0  N\n";
			print FILE "3,-3,3.0  N\n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4, 1,1.0  N\n";
			print FILE "4, 1,1.0  N\n";
			print FILE "4,-2,2.0  N\n";
			print FILE "4,-2,0.0  N\n";
		} elsif ($name =~ /Br/) {	
			print FILE "Br $add\n";
			print FILE "Ar 5  \n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3,-3,3.0  N\n";
			print FILE "3,-3,3.0  N\n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4, 1,1.0  N\n";
			print FILE "4, 1,1.0  N\n";
			print FILE "4,-2,2.0  N\n";
			print FILE "4,-2,1.0  N\n";
		} elsif ($name =~ /Kr/) {	
			print FILE "Kr $add\n";
			print FILE "Ar 5  \n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3, 2,2.0  N\n";
			print FILE "3,-3,3.0  N\n";
			print FILE "3,-3,3.0  N\n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4,-1,1.0  N\n";
			print FILE "4, 1,1.0  N\n";
			print FILE "4, 1,1.0  N\n";
			print FILE "4,-2,2.0  N\n";
			print FILE "4,-2,2.0  N\n";
		} elsif ($name =~ /Rb/)  {	
			print FILE "Rb $add\n";
			print FILE "Kr 1  \n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5,-1,0.0  N\n";
		} elsif ($name =~ /Sr/) {	
			print FILE "Sr $add\n";
			print FILE "Kr 1  \n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5,-1,1.0  N\n";
		} elsif ($name =~ /Y /) {	
			print FILE "Y $add\n";
			print FILE "Kr 2  \n";
			print FILE "4, 2,1.0  N\n";
			print FILE "4, 2,0.0  N\n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5,-1,1.0  N\n";
		} elsif ($name =~ /Zr/) {	
			print FILE "Zr $add\n";
			print FILE "Kr 2  \n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4, 2,0.0  N\n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5,-1,1.0  N\n";
		} elsif ($name =~ /Nb/) {	
			print FILE "Nb $add\n";
			print FILE "Kr 3  \n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4, 2,1.0  N\n";
			print FILE "4,-3,1.0  N\n";
			print FILE "4,-3,0.0  N\n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5,-1,0.0  N\n";
		} elsif ($name =~ /Mo/) {	
			print FILE "Mo $add\n";
			print FILE "Kr 3  \n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4, 2,1.0  N\n";
			print FILE "4,-3,2.0  N\n";
			print FILE "4,-3,0.0  N\n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5,-1,0.0  N\n";
		} elsif ($name =~ /Tc/) {	
			print FILE "Tc $add\n";
			print FILE "Kr 3  \n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4,-3,2.0  N\n";
			print FILE "4,-3,0.0  N\n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5,-1,0.0  N\n";
		} elsif ($name =~ /Ru/) {	
			print FILE "Ru $add\n";
			print FILE "Kr 3  \n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4,-3,2.0  N\n";
			print FILE "4,-3,1.0  N\n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5,-1,0.0  N\n";
		} elsif ($name =~ /Rh/) {	
			print FILE "Rh $add\n";
			print FILE "Kr 3  \n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4,-3,3.0  N\n";
			print FILE "4,-3,1.0  N\n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5,-1,0.0  N\n";
		} elsif ($name =~ /Pd/) {	
			print FILE "Pd $add\n";
			print FILE "Kr 3  \n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4,-3,3.0  N\n";
			print FILE "4,-3,2.0  N\n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5,-1,0.0  N\n";
		} elsif ($name =~ /Ag/) {	
			print FILE "Ag $add\n";
			print FILE "Kr 3  \n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4,-3,3.0  N\n";
			print FILE "4,-3,3.0  N\n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5,-1,0.0  N\n";
		} elsif ($name =~ /Cd/) {	
			print FILE "Cd $add\n";
			print FILE "Kr 3  \n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4,-3,3.0  N\n";
			print FILE "4,-3,3.0  N\n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5,-1,1.0  N\n";
		} elsif ($name =~ /In/) {	
			print FILE "In $add\n";
			print FILE "Kr 4  \n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4,-3,3.0  N\n";
			print FILE "4,-3,3.0  N\n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5, 1,1.0  N\n";
			print FILE "5, 1,0.0  N\n";
		} elsif ($name =~ /Sn/) {	
			print FILE "Sn $add\n";
			print FILE "Kr 5  \n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4,-3,3.0  N\n";
			print FILE "4,-3,3.0  N\n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5, 1,1.0  N\n";
			print FILE "5, 1,0.0  N\n";
			print FILE "5,-2,1.0  N\n";
			print FILE "5,-2,0.0  N\n";
		} elsif ($name =~ /Sb/) {	
			print FILE "Sb $add\n";
			print FILE "Kr 5  \n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4,-3,3.0  N\n";
			print FILE "4,-3,3.0  N\n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5, 1,1.0  N\n";
			print FILE "5, 1,0.0  N\n";
			print FILE "5,-2,2.0  N\n";
			print FILE "5,-2,0.0  N\n";
		} elsif ($name =~ /Te/) {	
			print FILE "Te $add\n";
			print FILE "Kr 5  \n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4,-3,3.0  N\n";
			print FILE "4,-3,3.0  N\n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5, 1,1.0  N\n";
			print FILE "5, 1,1.0  N\n";
			print FILE "5,-2,2.0  N\n";
			print FILE "5,-2,0.0  N\n";
		} elsif ($name =~ /I /) {	
			print FILE "I $add\n";
			print FILE "Kr 5  \n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4,-3,3.0  N\n";
			print FILE "4,-3,3.0  N\n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5, 1,1.0  N\n";
			print FILE "5, 1,1.0  N\n";
			print FILE "5,-2,2.0  N\n";
			print FILE "5,-2,1.0  N\n";
		} elsif ($name =~ /Xe/) {	
			print FILE "Xe $add\n";
			print FILE "Kr 5  \n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4, 2,2.0  N\n";
			print FILE "4,-3,3.0  N\n";
			print FILE "4,-3,3.0  N\n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5,-1,1.0  N\n";
			print FILE "5, 1,1.0  N\n";
			print FILE "5, 1,1.0  N\n";
			print FILE "5,-2,2.0  N\n";
			print FILE "5,-2,2.0  N\n";
		} elsif ($name =~ /Cs/) {	
			print FILE "Cs $add\n";
			print FILE "Xe 1  \n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,0.0  N\n";
		} elsif ($name =~ /Ba/) {	
			print FILE "Ba $add\n";
			print FILE "Xe 1  \n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
		} elsif ($name =~ /La/) {	
			print FILE "La $add\n";
			print FILE "Xe 2  \n";
			print FILE "5, 2,1.0  N\n";
			print FILE "5, 2,0.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
		} elsif ($name =~ /Ce/) {	
			print FILE "Ce $add\n";
			print FILE "Xe 4  \n";
			print FILE "4, 3,1.0  N\n";
			print FILE "4, 3,0.0  N\n";
			print FILE "4,-4,0.0  N\n";
			print FILE "4,-4,0.0  N\n";
			print FILE "5, 2,1.0  N\n";
			print FILE "5, 2,0.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
		} elsif ($name =~ /Pr/) {	
			print FILE "Pr $add\n";
			print FILE "Xe 4  \n";
			print FILE "4, 3,2.0  N\n";
			print FILE "4, 3,0.0  N\n";
			print FILE "4,-4,0.0  N\n";
			print FILE "4,-4,0.0  N\n";
      print FILE "5, 2,1.0  N\n";
      print FILE "5, 2,0.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
		} elsif ($name =~ /Nd/) {	
			print FILE "Nd $add\n";
			print FILE "Xe 4  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,0.0  N\n";
			print FILE "4,-4,0.0  N\n";
			print FILE "4,-4,0.0  N\n";
			print FILE "5, 2,1.0  N\n";
			print FILE "5, 2,0.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
		} elsif ($name =~ /Pm/) {	
			print FILE "Pm $add\n";
			print FILE "Xe 4  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,0.0  N\n";
			print FILE "4,-4,1.0  N\n";
			print FILE "4,-4,0.0  N\n";
      print FILE "5, 2,1.0  N\n";
      print FILE "5, 2,0.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
		} elsif ($name =~ /Sm/) {	
			print FILE "Sm $add\n";
			print FILE "Xe 4  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,0.0  N\n";
			print FILE "4,-4,2.0  N\n";
			print FILE "4,-4,0.0  N\n";
			print FILE "5, 2,1.0  N\n";
			print FILE "5, 2,0.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
		} elsif ($name =~ /Eu/) {	
			print FILE "Eu $add\n";
			print FILE "Xe 4  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,0.0  N\n";
			print FILE "4,-4,3.0  N\n";
			print FILE "4,-4,0.0  N\n";
			print FILE "5, 2,1.0  N\n";
			print FILE "5, 2,0.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
		} elsif ($name =~ /Gd/) {	
			print FILE "Gd $add\n";
			print FILE "Xe 4  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,0.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "4,-4,0.0  N\n";
			print FILE "5, 2,1.0  N\n";
			print FILE "5, 2,0.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
		} elsif ($name =~ /Tb/) {	
			print FILE "Tb $add\n";
			print FILE "Xe 4  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,1.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "4,-4,0.0  N\n";
      print FILE "5, 2,1.0  N\n";
      print FILE "5, 2,0.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
		} elsif ($name =~ /Dy/) {	
			print FILE "Dy $add\n";
			print FILE "Xe 4  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,2.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "4,-4,0.0  N\n";
      print FILE "5, 2,1.0  N\n";
      print FILE "5, 2,0.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
		} elsif ($name =~ /Ho/) {  
			print FILE "Ho $add\n";
      print FILE "Xe 4  \n";
      print FILE "4, 3,3.0  N\n";
      print FILE "4, 3,3.0  N\n";
      print FILE "4,-4,4.0  N\n";
      print FILE "4,-4,0.0  N\n";
      print FILE "6,-1,1.0  N\n";
      print FILE "6,-1,1.0  N\n";
      print FILE "5, 2,1.0  N\n";
      print FILE "5, 2,0.0  N\n";
		} elsif ($name =~ /Er/) {	
			print FILE "Er $add\n";
			print FILE "Xe 4  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "4,-4,1.0  N\n";
      print FILE "5, 2,1.0  N\n";
      print FILE "5, 2,0.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
		} elsif ($name =~ /Tm/) {	
			print FILE "Tm $add\n";
			print FILE "Xe 4  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "4,-4,2.0  N\n";
      print FILE "5, 2,1.0  N\n";
      print FILE "5, 2,0.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
		} elsif ($name =~ /Yb/) {	
			print FILE "Yb $add\n";
			print FILE "Xe 4  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "4,-4,3.0  N\n";
      print FILE "5, 2,1.0  N\n";
      print FILE "5, 2,0.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
		} elsif ($name =~ /Lu/) {	
			print FILE "Lu $add\n";
			print FILE "Xe 4  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "5, 2,1.0  N\n";
			print FILE "5, 2,0.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
		} elsif ($name =~ /Hf/) {	
			print FILE "Hf $add\n";
			print FILE "Xe 4  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5, 2,0.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
		} elsif ($name =~ /Ta/) {	
			print FILE "Ta $add\n";
			print FILE "Xe 4  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5, 2,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
		} elsif ($name =~ /W /) {	
			print FILE "W $add\n";
			print FILE "Xe 5  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5, 2,1.0  N\n";
			print FILE "5,-3,1.0  N\n";
			print FILE "5,-3,0.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
		} elsif ($name =~ /Re/) {	
			print FILE "Re $add\n";
			print FILE "Xe 5  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5,-3,1.0  N\n";
			print FILE "5,-3,0.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
		} elsif ($name =~ /Os/) {	
			print FILE "Os $add\n";
			print FILE "Xe 5  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5,-3,2.0  N\n";
			print FILE "5,-3,0.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
		} elsif ($name =~ /Ir/) {	
			print FILE "Ir $add\n";
			print FILE "Xe 5  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5,-3,2.0  N\n";
			print FILE "5,-3,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
		} elsif ($name =~ /Pt/) {	
			print FILE "Pt $add\n";
			print FILE "Xe 5  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5,-3,3.0  N\n";
			print FILE "5,-3,2.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,0.0  N\n";
		} elsif ($name =~ /Au/) {	
			print FILE "Au $add\n";
			print FILE "Xe 5  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5,-3,3.0  N\n";
			print FILE "5,-3,3.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,0.0  N\n";
		} elsif ($name =~ /Hg/) {	
			print FILE "Hg $add\n";
			print FILE "Xe 5  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5,-3,3.0  N\n";
			print FILE "5,-3,3.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
		} elsif ($name =~ /Tl/) {	
			print FILE "Tl $add\n";
			print FILE "Xe 6  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5,-3,3.0  N\n";
			print FILE "5,-3,3.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6, 1,1.0  N\n";
			print FILE "6, 1,0.0  N\n";
		} elsif ($name =~ /Pb/) {	
			print FILE "Pb $add\n";
			print FILE "Xe 7  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5,-3,3.0  N\n";
			print FILE "5,-3,3.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6, 1,1.0  N\n";
			print FILE "6, 1,0.0  N\n";
			print FILE "6,-2,1.0  N\n";
			print FILE "6,-2,0.0  N\n";
		} elsif ($name =~ /Bi/) {	
			print FILE "Bi $add\n";
			print FILE "Xe 7  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5,-3,3.0  N\n";
			print FILE "5,-3,3.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6, 1,1.0  N\n";
			print FILE "6, 1,1.0  N\n";
			print FILE "6,-2,1.0  N\n";
			print FILE "6,-2,0.0  N\n";
		} elsif ($name =~ /Po/) {	
			print FILE "Po $add\n";
			print FILE "Xe 7  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5,-3,3.0  N\n";
			print FILE "5,-3,3.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6, 1,1.0  N\n";
			print FILE "6, 1,1.0  N\n";
			print FILE "6,-2,2.0  N\n";
			print FILE "6,-2,0.0  N\n";
		} elsif ($name =~ /At/) {	
			print FILE "At $add\n";
			print FILE "Xe 7  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5,-3,3.0  N\n";
			print FILE "5,-3,3.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6, 1,1.0  N\n";
			print FILE "6, 1,1.0  N\n";
			print FILE "6,-2,2.0  N\n";
			print FILE "6,-2,1.0  N\n";
		} elsif ($name =~ /Rn/) {	
			print FILE "Rn $add\n";
			print FILE "Xe 7  \n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4, 3,3.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "4,-4,4.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5, 2,2.0  N\n";
			print FILE "5,-3,3.0  N\n";
			print FILE "5,-3,3.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6,-1,1.0  N\n";
			print FILE "6, 1,1.0  N\n";
			print FILE "6, 1,1.0  N\n";
			print FILE "6,-2,2.0  N\n";
			print FILE "6,-2,2.0  N\n";
		} elsif ($name =~ /Fr/) {	
			print FILE "Fr $add\n";
			print FILE "Rn 1  \n";
			print FILE "7,-1,1.0  N\n";
			print FILE "7,-1,0.0  N\n";
		} elsif ($name =~ /Ra/) {	
			print FILE "Ra $add\n";
			print FILE "Rn 1  \n";
			print FILE "7,-1,1.0  N\n";
			print FILE "7,-1,1.0  N\n";
		} elsif ($name =~ /Ac/) {	
			print FILE "Ac $add\n";
			print FILE "Rn 2  \n";
			print FILE "6, 2,1.0  N\n";
			print FILE "6, 2,0.0  N\n";
			print FILE "7,-1,1.0  N\n";
			print FILE "7,-1,1.0  N\n";
		} elsif ($name =~ /Th/) {	
			print FILE "Th $add\n";
			print FILE "Rn 2  \n";
			print FILE "6, 2,2.0  N\n";
			print FILE "6, 2,0.0  N\n";
			print FILE "7,-1,1.0  N\n";
			print FILE "7,-1,1.0  N\n";
		} elsif ($name =~ /Pa/) {	
			print FILE "Pa $add\n";
			print FILE "Rn 3  \n";
			print FILE "5, 3,2.0  N\n";
			print FILE "5, 3,0.0  N\n";
			print FILE "6, 2,1.0  N\n";
			print FILE "6, 2,0.0  N\n";
			print FILE "7,-1,1.0  N\n";
			print FILE "7,-1,1.0  N\n";
		} elsif ($name =~ /U /) {	
			print FILE "U $add\n";
			print FILE "Rn 3  \n";
			print FILE "5, 3,3.0  N\n";
			print FILE "5, 3,0.0  N\n";
			print FILE "6, 2,1.0  N\n";
			print FILE "6, 2,0.0  N\n";
			print FILE "7,-1,1.0  N\n";
			print FILE "7,-1,1.0  N\n";
		} elsif ($name =~ /Np/)  {	
			print FILE "Np $add\n";
			print FILE "Rn 4  \n";
			print FILE "5, 3,3.0  N\n";
			print FILE "5, 3,0.0  N\n";
			print FILE "5,-4,1.0  N\n";
			print FILE "5,-4,0.0  N\n";
			print FILE "6, 2,1.0  N\n";
			print FILE "6, 2,0.0  N\n";
			print FILE "7,-1,1.0  N\n";
			print FILE "7,-1,1.0  N\n";
		} elsif ($name =~ /Pu/) {	
			print FILE "Pu $add\n";
			print FILE "Rn 4  \n";
			print FILE "5, 3,3.0  N\n";
			print FILE "5, 3,0.0  N\n";
			print FILE "5,-4,2.0  N\n";
			print FILE "5,-4,0.0  N\n";
      print FILE "6, 2,1.0  N\n";
      print FILE "6, 2,0.0  N\n";
			print FILE "7,-1,1.0  N\n";
			print FILE "7,-1,1.0  N\n";
		} elsif ($name =~ /Am/) {	
			print FILE "Am $add\n";
			print FILE "Rn 4  \n";
			print FILE "5, 3,3.0  N\n";
			print FILE "5, 3,0.0  N\n";
			print FILE "5,-4,3.0  N\n";
			print FILE "5,-4,0.0  N\n";
      print FILE "6, 2,1.0  N\n";
      print FILE "6, 2,0.0  N\n";
			print FILE "7,-1,1.0  N\n";
			print FILE "7,-1,1.0  N\n";
		} elsif ($name =~ /Cm/) {	
			print FILE "Cm $add\n";
			print FILE "Rn 4  \n";
			print FILE "5, 3,3.0  N\n";
			print FILE "5, 3,0.0  N\n";
			print FILE "5,-4,4.0  N\n";
			print FILE "5,-4,0.0  N\n";
			print FILE "6, 2,1.0  N\n";
			print FILE "6, 2,0.0  N\n";
			print FILE "7,-1,1.0  N\n";
			print FILE "7,-1,1.0  N\n";
		} elsif ($name =~ /Bk/) {	
			print FILE "Bk $add\n";
			print FILE "Rn 3  \n";
			print FILE "5, 3,3.0  N\n";
			print FILE "5, 3,2.0  N\n";
			print FILE "5,-4,4.0  N\n";
			print FILE "5,-4,0.0  N\n";
			print FILE "7,-1,1.0  N\n";
			print FILE "7,-1,1.0  N\n";
		} elsif ($name =~ /Cf/) {	
			print FILE "Cf $add\n";
			print FILE "Rn 3  \n";
			print FILE "5, 3,3.0  N\n";
			print FILE "5, 3,3.0  N\n";
			print FILE "5,-4,4.0  N\n";
			print FILE "5,-4,0.0  N\n";
			print FILE "7,-1,1.0  N\n";
			print FILE "7,-1,1.0  N\n";
		} elsif ($name =~ /Es/) {	
			print FILE "Es $add\n";
			print FILE "Rn 3  \n";
			print FILE "5, 3,3.0  N\n";
			print FILE "5, 3,3.0  N\n";
			print FILE "5,-4,4.0  N\n";
			print FILE "5,-4,1.0  N\n";
			print FILE "7,-1,1.0  N\n";
			print FILE "7,-1,1.0  N\n";
		} elsif ($name =~ /Fm/) {	
			print FILE "Fm $add\n";
			print FILE "Rn 3  \n";
			print FILE "5, 3,3.0  N\n";
			print FILE "5, 3,3.0  N\n";
			print FILE "5,-4,4.0  N\n";
			print FILE "5,-4,2.0  N\n";
			print FILE "7,-1,1.0  N\n";
			print FILE "7,-1,1.0  N\n";
		} elsif ($name =~ /Md/) {	
			print FILE "Md $add\n";
			print FILE "Rn 3  \n";
			print FILE "5, 3,3.0  N\n";
			print FILE "5, 3,3.0  N\n";
			print FILE "5,-4,4.0  N\n";
			print FILE "5,-4,3.0  N\n";
			print FILE "7,-1,1.0  N\n";
			print FILE "7,-1,1.0  N\n";
		} elsif ($name =~ /No/) {	
			print FILE "No $add\n";
			print FILE "Rn 3  \n";
			print FILE "5, 3,3.0  N\n";
			print FILE "5, 3,3.0  N\n";
			print FILE "5,-4,4.0  N\n";
			print FILE "5,-4,4.0  N\n";
			print FILE "7,-1,1.0  N\n";
			print FILE "7,-1,1.0  N\n";
		} elsif ($name =~ /Lr/) {	
			print FILE "Lr $add\n";
			print FILE "Rn 4  \n";
			print FILE "5, 3,3.0  N\n";
			print FILE "5, 3,3.0  N\n";
			print FILE "5,-4,4.0  N\n";
			print FILE "5,-4,4.0  N\n";
			print FILE "6, 2,1.0  N\n";
			print FILE "6, 2,0.0  N\n";
			print FILE "7,-1,1.0  N\n";
			print FILE "7,-1,1.0  N\n";
		} else {
			    print FILE "WARNING: specified Element not in database!\n";
		}		
	}
	print FILE "****     End of Input\n";
	print FILE "****     End of Input\n";
	close(FILE);
}	      	

# check if inst inst inst_tmp file differ

#    catch { [exec diff $case.inst $case.inst_tmp | wc -l ] } diff

sub Element {
#	($elem,$z,$ro) = shift(@_);
	$elem .= " " if (length($elem)==1);
##	if ($z =~ //) {

		if($elem =~ /Ac/) {$z = 89.}
		if($elem =~ /Ag/) {$z = 47.}
		if($elem =~ /Al/) {$z = 13.}
		if($elem =~ /Am/) {$z = 95.}
		if($elem =~ /Ar/) {$z = 18.}
		if($elem =~ /As/) {$z = 33.}
		if($elem =~ /At/) {$z = 85.}
		if($elem =~ /Au/) {$z = 79.}
		if($elem =~ /B /) {$z =  5.}
		if($elem =~ /Ba/) {$z = 56.}
		if($elem =~ /Be/) {$z =  4.}
		if($elem =~ /Bi/) {$z = 83.}
		if($elem =~ /Bk/) {$z = 97.}
		if($elem =~ /Br/) {$z = 35.}
		if($elem =~ /C /) {$z =  6.}
		if($elem =~ /Ca/) {$z = 20.}
		if($elem =~ /Cd/) {$z = 48.}
		if($elem =~ /Ce/) {$z = 58.}
		if($elem =~ /Cf/) {$z = 98.}
		if($elem =~ /Cl/) {$z = 17.}
		if($elem =~ /Cm/) {$z = 96.}
		if($elem =~ /Co/) {$z = 27.}
		if($elem =~ /Cr/) {$z = 24.}  
		if($elem =~ /Cs/) {$z = 55.}
		if($elem =~ /Cu/) {$z = 29.}
		if($elem =~ /D /) {$z =  1.}
		if($elem =~ /Dy/) {$z = 66.}
		if($elem =~ /Er/) {$z = 68.}
		if($elem =~ /Es/) {$z = 99.}
		if($elem =~ /Eu/) {$z = 63.}
		if($elem =~ /F /) {$z =  9.}
		if($elem =~ /Fe/) {$z = 26.}
		if($elem =~ /Fm/) {$z = 100.}
		if($elem =~ /Fr/) {$z = 87.}
		if($elem =~ /Ga/) {$z = 31.}
		if($elem =~ /Gd/) {$z = 64.}
		if($elem =~ /Ge/) {$z = 32.}
		if($elem =~ /H /) {$z =  1.}
		if($elem =~ /He/) {$z =  2.}
		if($elem =~ /Hf/) {$z = 72.}
		if($elem =~ /Hg/) {$z = 80.}
		if($elem =~ /Ho/) {$z = 67.}
		if($elem =~ /I /) {$z = 53.}
		if($elem =~ /In/) {$z = 49.}
		if($elem =~ /Ir/) {$z = 77.}
		if($elem =~ /K /) {$z = 19.}
		if($elem =~ /Kr/) {$z = 36.}
		if($elem =~ /La/) {$z = 57.}
		if($elem =~ /Li/) {$z =  3.}
		if($elem =~ /Lr/) {$z = 103.}  
		if($elem =~ /Lu/) {$z = 71.}
		if($elem =~ /Md/) {$z = 101.}
		if($elem =~ /Mg/) {$z = 12.}
		if($elem =~ /Mn/) {$z = 25.}
		if($elem =~ /Mo/) {$z = 42.}
		if($elem =~ /N /) {$z =  7.}
		if($elem =~ /Na/) {$z = 11.}
		if($elem =~ /Nb/) {$z = 41.}
		if($elem =~ /Nd/) {$z = 60.}
		if($elem =~ /Ne/) {$z = 10.}
		if($elem =~ /Ni/) {$z = 28.}
		if($elem =~ /No/) {$z = 102.}
		if($elem =~ /Np/) {$z = 93.}
		if($elem =~ /O /) {$z =  8.}
		if($elem =~ /Os/) {$z = 76.}
		if($elem =~ /P /) {$z = 15.}
		if($elem =~ /Pa/) {$z = 91.}
		if($elem =~ /Pb/) {$z = 82.}
		if($elem =~ /Pd/) {$z = 46.}
		if($elem =~ /Pm/) {$z = 61.}
		if($elem =~ /Po/) {$z = 84.}
		if($elem =~ /Pr/) {$z = 59.}
		if($elem =~ /Pt/) {$z = 78.}
		if($elem =~ /Pu/) {$z = 94.}
		if($elem =~ /Ra/) {$z = 88.}
		if($elem =~ /Rb/) {$z = 37.}
		if($elem =~ /Re/) {$z = 75.}
		if($elem =~ /Rh/) {$z = 45.}
		if($elem =~ /Rn/) {$z = 86.}
		if($elem =~ /Ru/) {$z = 44.}
		if($elem =~ /S /) {$z = 16.}
		if($elem =~ /Sb/) {$z = 51.}
		if($elem =~ /Sc/) {$z = 21.}
		if($elem =~ /Se/) {$z = 34.}
		if($elem =~ /Si/) {$z = 14.}
		if($elem =~ /Sm/) {$z = 62.}
		if($elem =~ /Sn/) {$z = 50.}
		if($elem =~ /Sr/) {$z = 38.}
		if($elem =~ /Ta/) {$z = 73.}
		if($elem =~ /Tb/) {$z = 65.}
		if($elem =~ /Tc/) {$z = 43.}
		if($elem =~ /Te/) {$z = 52.}
		if($elem =~ /Th/) {$z = 90.}
		if($elem =~ /Ti/) {$z = 22.}
		if($elem =~ /Tl/) {$z = 81.}
		if($elem =~ /Tm/) {$z = 69.}
		if($elem =~ /U /) {$z = 92.}
		if($elem =~ /V /) {$z = 23.}
		if($elem =~ /W /) {$z = 74.}
		if($elem =~ /Xe/) {$z = 54.}
		if($elem =~ /Y/) {$z = 39.}
		if($elem =~ /Yb/) {$z = 70.}
		if($elem =~ /Zn/) {$z = 30.}
		if($elem =~ /Zr/) {$z = 40.}

		$ro = 0.000005;
		if ($z <= 71.) { $ro = 0.00001 }
		if ($z <= 36.) { $ro = 0.00005 }
		if ($z <= 18.) { $ro = 0.0001  }
##	}
}



#######################################
# SPACEGROUPS                         #
#######################################

@lattype = qw/P F B CXY CYZ CXZ R H 1_P1 2_P-1 3_P2 3_P2 3_P2 4_P21 4_P21 4_P21 5_B2 6_Pm 6_Pm 6_Pm 7_Pc 7_Pa 7_Pb 7_Pb 7_Pc 7_Pa 7_Pn 7_Pn 7_Pn 8_Bm 9_Bb 10_P2\/m 10_P2\/m 10_P2\/m 11_P21\/m 11_P21\/m 11_P21\/m 12_B2\/m 13_P2\/c 13_P2\/a 13_P2\/b 13_P2\/b 13_P2\/c 13_P2\/a 13_P2\/n 13_P2\/n 13_P2\/n 14_P21\/c 14_P21\/a 14_P21\/b 14_P21\/b 14_P21\/c 14_P21\/a 14_P21\/n 14_P21\/n 14_P21\/n 15_B2\/b 16_P222 17_P2221 17_P2122 17_P2212 18_P21212 18_P22121 18_P21221 19_P212121 20_C2221 20_A2122 20_B2212 21_C222 21_A222 21_B222 22_F222 23_I222 24_I212121 25_Pmm2 25_P2mm 25_Pm2m 26_Pmc21 26_P21ma 26_Pb21m 26_Pcm21 26_P21am 26_Pm21b 27_Pcc2 27_P2aa 27_Pb2b 28_Pma2 28_P2mb 28_Pc2m 28_Pbm2 28_P2cm 28_Pm2a 29_Pca21 29_P21ab 29_Pc21b 29_Pbc21 29_P21ca 29_Pb21a 30_Pnc2 30_P2na 30_Pb2n 30_Pcn2 30_P2an 30_Pn2b 31_Pmn21 31_P21mn 31_Pn21m 31_Pnm21 31_P21nm 31_Pm21n 32_Pba2 32_P2cb 32_Pc2a 33_Pna21 33_P21nb 33_Pc21n 33_Pbn21 33_P21cn 33_Pn21a 34_Pnn2 34_P2nn 34_Pn2n 35_Cmm2 35_A2mm 35_Bm2m 36_Cmc21 36_A21ma 36_Bb21m 36_Ccm21 36_A21am 36_Bm21b 37_Ccc2 37_A2aa 37_Bb2b 38_Amm2 38_B2mm 38_Cm2m 39_Abm2 39_B2cm 39_Cm2a 39_Bma2 39_C2mb 39_Ac2m 40_Ama2 40_B2mb 40_Cc2m 40_Bbm2 40_C2cm 40_Am2a 41_Aba2 41_B2cb 41_Cc2a 41_Bba2 41_C2cb 41_Ac2a 42_Fmm2 42_F2mm 42_Fm2m 43_Fdd2 43_F2dd 43_Fd2d 44_Imm2 44_I2mm 44_Im2m 45_Iba2 45_I2cb 45_Ic2a 46_Ima2 46_I2mb 46_Ic2m 46_Ibm2 46_I2cm 46_Im2a 47_Pmmm 48_Pnnn 49_Pccm 49_Pmaa 49_Pbmb 50_Pban 50_Pncb 50_Pcna 51_Pmma 51_Pbmm 51_Pmcm 51_Pmam 51_Pmmb 51_Pcmm 52_Pnna 52_Pbnn 52_Pncn 52_Pnan 52_Pnnb 52_Pcnn 53_Pmna 53_Pbmn 53_Pncm 53_Pman 53_Pnmb 53_Pcnm 54_Pcca 54_Pbaa 54_Pbcb 54_Pbab 54_Pccb 54_Pcaa 55_Pbam 55_Pmcb 55_Pcma 56_Pccn 56_Pnaa 56_Pbnb 57_Pbcm 57_Pmca 57_Pbma 57_Pcmb 57_Pcam 57_Pmab 58_Pnnm 58_Pmnn 58_Pnmn 59_Pmmn 59_Pnmm 59_Pmnm 60_Pbcn 60_Pnca 60_Pbna 60_Pcnb 60_Pcan 60_Pnab 61_Pbca 61_Pcab 62_Pnma 62_Pbnm 62_Pmcn 62_Pnam 62_Pmnb 62_Pcmn 63_Cmcm 63_Amma 63_Bbmm 63_Bmmb 63_Ccmm 63_Amam 64_Cmca 64_Abma 64_Bbcm 64_Bmab 64_Ccmb 64_Acam 65_Cmmm 65_Ammm 65_Bmmm 66_Cccm 66_Amaa 66_Bbmb 67_Cmma 67_Abmm 67_Bmcm 67_Bmam 67_Cmmb 67_Acmm 68_Ccca 68_Abaa 68_Bbcb 68_Bbab 68_Cccb 68_Acaa 69_Fmmm 70_Fddd 71_Immm 72_Ibam 72_Imcb 72_Icma 73_Ibca 73_Icab 74_Imma 74_Ibmm 74_Imcm 74_Imam 74_Immb 74_Icmm 75_P4 76_P41 77_P42 78_P43 79_I4 80_I41 81_P-4 82_I-4 83_P4\/m 84_P42\/m 85_P4\/n 86_P42\/n 87_I4\/m 88_I41\/a 89_P422 90_P4212 91_P4122 92_P41212 93_P4222 94_P42212 95_P4322 96_P43212 97_I422 98_I4122 99_P4mm 100_P4bm 101_P42cm 102_P42nm 103_P4cc 104_P4nc 105_P42mc 106_P42bc 107_I4mm 108_I4cm 109_I41md 110_I41cd 111_P-42m 112_P-42c 113_P-421m 114_P-421c 115_P-4m2 116_P-4c2 117_P-4b2 118_P-4n2 119_I-4m2 120_I-4c2 121_I-42m 122_I-42d 123_P4\/mmm 124_P4\/mcc 125_P4\/nbm 126_P4\/nnc 127_P4\/mbm 128_P4\/mnc 129_P4\/nmm 130_P4\/ncc 131_P42\/mmc 132_P42\/mcm 133_P42\/nbc 134_P42\/nnm 135_P42\/mbc 136_P42\/mnm 137_P42\/nmc 138_P42\/ncm 139_I4\/mmm 140_I4\/mcm 141_I41\/amd 142_I41\/acd 143_P3 144_P31 145_P32 146_R3 147_P-3 148_R-3 149_P312 150_P321 151_P3112 152_P3121 153_P3212 154_P3221 155_R32 156_P3m1 157_P31m 158_P3c1 159_P31c 160_R3m 161_R3c 162_P-31m 163_P-31c 164_P-3m1 165_P-3c1 166_R-3m 167_R-3c 168_P6 169_P61 170_P65 171_P62 172_P64 173_P63 174_P-6 175_P6\/m 176_P63\/m 177_P622 178_P6122 179_P6522 180_P6222 181_P6422 182_P6322 183_P6mm 184_P6cc 185_P63cm 186_P63mc 187_P-6m2 188_P-6c2 189_P-62m 190_P-62c 191_P6\/mmm 192_P6\/mcc 193_P63\/mcm 194_P63\/mmc 195_P23 196_F23 197_I23 198_P213 199_I213 200_Pm-3 201_Pn-3 202_Fm-3 203_Fd-3 204_Im-3 205_Pa-3 206_Ia-3 207_P432 208_P4232 209_F432 210_F4132 211_I432 212_P4332 213_P4132 214_I4132 215_P-43m 216_F-43m 217_I-43m 218_P-43n 219_F-43c 220_I-43d 221_Pm-3m 222_Pn-3n 223_Pm-3n 224_Pn-3m 225_Fm-3m 226_Fm-3c 227_Fd-3m 228_Fd-3c 229_Im-3m 230_Ia-3d/; 
#original list below
#@lattype = qw/P F B CXY CYZ CXZ R H 1_P1 2_P-1 3_P2 3_P2 3_P2 4_P21 4_P21 4_P21 5_C2 5_A2 5_B2 5_B2 5_C2 5_A2 6_Pm 6_Pm 6_Pm 7_Pc 7_Pa 7_Pb 7_Pb 7_Pc 7_Pa 7_Pn 7_Pn 7_Pn 8_Cm 8_Am 8_Bm 8_Bm 8_Cm 8_Am 9_Cc 9_Aa 9_Bb 9_Bb 9_Cc 9_Aa 10_P2\/m 10_P2\/m 10_P2\/m 11_P21\/m 11_P21\/m 11_P21\/m 12_C2\/m 12_A2\/m 12_B2\/m 12_B2\/m 12_C2\/m 12_A2\/m 13_P2\/c 13_P2\/a 13_P2\/b 13_P2\/b 13_P2\/c 13_P2\/a 13_P2\/n 13_P2\/n 13_P2\/n 14_P21\/c 14_P21\/a 14_P21\/b 14_P21\/b 14_P21\/c 14_P21\/a 14_P21\/n 14_P21\/n 14_P21\/n 15_C2\/c 15_A2\/a 15_B2\/b 15_B2\/b 15_C2\/c 15_A2\/a 16_P222 17_P2221 17_P2122 17_P2212 18_P21212 18_P22121 18_P21221 19_P212121 20_C2221 20_A2122 20_B2212 21_C222 21_A222 21_B222 22_F222 23_I222 24_I212121 25_Pmm2 25_P2mm 25_Pm2m 26_Pmc21 26_P21ma 26_Pb21m 26_Pcm21 26_P21am 26_Pm21b 27_Pcc2 27_P2aa 27_Pb2b 28_Pma2 28_P2mb 28_Pc2m 28_Pbm2 28_P2cm 28_Pm2a 29_Pca21 29_P21ab 29_Pc21b 29_Pbc21 29_P21ca 29_Pb21a 30_Pnc2 30_P2na 30_Pb2n 30_Pcn2 30_P2an 30_Pn2b 31_Pmn21 31_P21mn 31_Pn21m 31_Pnm21 31_P21nm 31_Pm21n 32_Pba2 32_P2cb 32_Pc2a 33_Pna21 33_P21nb 33_Pc21n 33_Pbn21 33_P21cn 33_Pn21a 34_Pnn2 34_P2nn 34_Pn2n 35_Cmm2 35_A2mm 35_Bm2m 36_Cmc21 36_A21ma 36_Bb21m 36_Ccm21 36_A21am 36_Bm21b 37_Ccc2 37_A2aa 37_Bb2b 38_Amm2 38_B2mm 38_Cm2m 39_Abm2 39_B2cm 39_Cm2a 39_Bma2 39_C2mb 39_Ac2m 40_Ama2 40_B2mb 40_Cc2m 40_Bbm2 40_C2cm 40_Am2a 41_Aba2 41_B2cb 41_Cc2a 41_Bba2 41_C2cb 41_Ac2a 42_Fmm2 42_F2mm 42_Fm2m 43_Fdd2 43_F2dd 43_Fd2d 44_Imm2 44_I2mm 44_Im2m 45_Iba2 45_I2cb 45_Ic2a 46_Ima2 46_I2mb 46_Ic2m 46_Ibm2 46_I2cm 46_Im2a 47_Pmmm 48_Pnnn 49_Pccm 49_Pmaa 49_Pbmb 50_Pban 50_Pncb 50_Pcna 51_Pmma 51_Pbmm 51_Pmcm 51_Pmam 51_Pmmb 51_Pcmm 52_Pnna 52_Pbnn 52_Pncn 52_Pnan 52_Pnnb 52_Pcnn 53_Pmna 53_Pbmn 53_Pncm 53_Pman 53_Pnmb 53_Pcnm 54_Pcca 54_Pbaa 54_Pbcb 54_Pbab 54_Pccb 54_Pcaa 55_Pbam 55_Pmcb 55_Pcma 56_Pccn 56_Pnaa 56_Pbnb 57_Pbcm 57_Pmca 57_Pbma 57_Pcmb 57_Pcam 57_Pmab 58_Pnnm 58_Pmnn 58_Pnmn 59_Pmmn 59_Pnmm 59_Pmnm 60_Pbcn 60_Pnca 60_Pbna 60_Pcnb 60_Pcan 60_Pnab 61_Pbca 61_Pcab 62_Pnma 62_Pbnm 62_Pmcn 62_Pnam 62_Pmnb 62_Pcmn 63_Cmcm 63_Amma 63_Bbmm 63_Bmmb 63_Ccmm 63_Amam 64_Cmca 64_Abma 64_Bbcm 64_Bmab 64_Ccmb 64_Acam 65_Cmmm 65_Ammm 65_Bmmm 66_Cccm 66_Amaa 66_Bbmb 67_Cmma 67_Abmm 67_Bmcm 67_Bmam 67_Cmmb 67_Acmm 68_Ccca 68_Abaa 68_Bbcb 68_Bbab 68_Cccb 68_Acaa 69_Fmmm 70_Fddd 71_Immm 72_Ibam 72_Imcb 72_Icma 73_Ibca 73_Icab 74_Imma 74_Ibmm 74_Imcm 74_Imam 74_Immb 74_Icmm 75_P4 76_P41 77_P42 78_P43 79_I4 80_I41 81_P-4 82_I-4 83_P4\/m 84_P42\/m 85_P4\/n 86_P42\/n 87_I4\/m 88_I41\/a 89_P422 90_P4212 91_P4122 92_P41212 93_P4222 94_P42212 95_P4322 96_P43212 97_I422 98_I4122 99_P4mm 100_P4bm 101_P42cm 102_P42nm 103_P4cc 104_P4nc 105_P42mc 106_P42bc 107_I4mm 108_I4cm 109_I41md 110_I41cd 111_P-42m 112_P-42c 113_P-421m 114_P-421c 115_P-4m2 116_P-4c2 117_P-4b2 118_P-4n2 119_I-4m2 120_I-4c2 121_I-42m 122_I-42d 123_P4\/mmm 124_P4\/mcc 125_P4\/nbm 126_P4\/nnc 127_P4\/mbm 128_P4\/mnc 129_P4\/nmm 130_P4\/ncc 131_P42\/mmc 132_P42\/mcm 133_P42\/nbc 134_P42\/nnm 135_P42\/mbc 136_P42\/mnm 137_P42\/nmc 138_P42\/ncm 139_I4\/mmm 140_I4\/mcm 141_I41\/amd 142_I41\/acd 143_P3 144_P31 145_P32 146_R3 147_P-3 148_R-3 149_P312 150_P321 151_P3112 152_P3121 153_P3212 154_P3221 155_R32 156_P3m1 157_P31m 158_P3c1 159_P31c 160_R3m 161_R3c 162_P-31m 163_P-31c 164_P-3m1 165_P-3c1 166_R-3m 167_R-3c 168_P6 169_P61 170_P65 171_P62 172_P64 173_P63 174_P-6 175_P6\/m 176_P63\/m 177_P622 178_P6122 179_P6522 180_P6222 181_P6422 182_P6322 183_P6mm 184_P6cc 185_P63cm 186_P63mc 187_P-6m2 188_P-6c2 189_P-62m 190_P-62c 191_P6\/mmm 192_P6\/mcc 193_P63\/mcm 194_P63\/mmc 195_P23 196_F23 197_I23 198_P213 199_I213 200_Pm-3 201_Pn-3 202_Fm-3 203_Fd-3 204_Im-3 205_Pa-3 206_Ia-3 207_P432 208_P4232 209_F432 210_F4132 211_I432 212_P4332 213_P4132 214_I4132 215_P-43m 216_F-43m 217_I-43m 218_P-43n 219_F-43c 220_I-43d 221_Pm-3m 222_Pn-3n 223_Pm-3n 224_Pn-3m 225_Fm-3m 226_Fm-3c 227_Fd-3m 228_Fd-3c 229_Im-3m 230_Ia-3d/; 



sub SG2lat {
	@lattic = qw/P F B CXY CYZ CXZ R H/;
	my $mylat = $s_lattice;

	if($mylat =~ /P1/) {$s_lattice = "P"};
	if($mylat =~ /P-1/) {$s_lattice = "P"};
	if($mylat =~ /P2/) {$s_lattice = "P"};
	if($mylat =~ /P2/) {$s_lattice = "P"};
	if($mylat =~ /P2/) {$s_lattice = "P"};
	if($mylat =~ /P21/) {$s_lattice = "P"};
	if($mylat =~ /P21/) {$s_lattice = "P"};
	if($mylat =~ /P21/) {$s_lattice = "P"};
	if($mylat =~ /C2/) {$s_lattice = "CXY"};
	if($mylat =~ /A2/) {$s_lattice = "CYZ"};
	if($mylat =~ /B2/) {$s_lattice = "CXZ"};
	if($mylat =~ /B2/) {$s_lattice = "CXZ"};
	if($mylat =~ /C2/) {$s_lattice = "CXY"};
	if($mylat =~ /A2/) {$s_lattice = "CYZ"};
	if($mylat =~ /Pm/) {$s_lattice = "P"};
	if($mylat =~ /Pm/) {$s_lattice = "P"};
	if($mylat =~ /Pm/) {$s_lattice = "P"};
	if($mylat =~ /Pc/) {$s_lattice = "P"};
	if($mylat =~ /Pa/) {$s_lattice = "P"};
	if($mylat =~ /Pb/) {$s_lattice = "P"};
	if($mylat =~ /Pb/) {$s_lattice = "P"};
	if($mylat =~ /Pc/) {$s_lattice = "P"};
	if($mylat =~ /Pa/) {$s_lattice = "P"};
	if($mylat =~ /Pn/) {$s_lattice = "P"};
	if($mylat =~ /Pn/) {$s_lattice = "P"};
	if($mylat =~ /Pn/) {$s_lattice = "P"};
	if($mylat =~ /Cm/) {$s_lattice = "CXY"};
	if($mylat =~ /Am/) {$s_lattice = "CYZ"};
	if($mylat =~ /Bm/) {$s_lattice = "CXZ"};
	if($mylat =~ /Bm/) {$s_lattice = "CXZ"};
	if($mylat =~ /Cm/) {$s_lattice = "CXY"};
	if($mylat =~ /Am/) {$s_lattice = "CYZ"};
	if($mylat =~ /Cc/) {$s_lattice = "CXY"};
	if($mylat =~ /Aa/) {$s_lattice = "CYZ"};
	if($mylat =~ /Bb/) {$s_lattice = "CXZ"};
	if($mylat =~ /Bb/) {$s_lattice = "CXZ"};
	if($mylat =~ /Cc/) {$s_lattice = "CXY"};
	if($mylat =~ /Aa/) {$s_lattice = "CYZ"};
	if($mylat =~ /P2\/m/) {$s_lattice = "P"};
	if($mylat =~ /P2\/m/) {$s_lattice = "P"};
	if($mylat =~ /P2\/m/) {$s_lattice = "P"};
	if($mylat =~ /P21\/m/) {$s_lattice = "P"};
	if($mylat =~ /P21\/m/) {$s_lattice = "P"};
	if($mylat =~ /P21\/m/) {$s_lattice = "P"};
	if($mylat =~ /C2\/m/) {$s_lattice = "CXY"};
	if($mylat =~ /A2\/m/) {$s_lattice = "CYZ"};
	if($mylat =~ /B2\/m/) {$s_lattice = "CXZ"};
	if($mylat =~ /B2\/m/) {$s_lattice = "CXZ"};
	if($mylat =~ /C2\/m/) {$s_lattice = "CXY"};
	if($mylat =~ /A2\/m/) {$s_lattice = "CYZ"};
	if($mylat =~ /P2\/c/) {$s_lattice = "P"};
	if($mylat =~ /P2\/a/) {$s_lattice = "P"};
	if($mylat =~ /P2\/b/) {$s_lattice = "P"};
	if($mylat =~ /P2\/b/) {$s_lattice = "P"};
	if($mylat =~ /P2\/c/) {$s_lattice = "P"};
	if($mylat =~ /P2\/a/) {$s_lattice = "P"};
	if($mylat =~ /P2\/n/) {$s_lattice = "P"};
	if($mylat =~ /P2\/n/) {$s_lattice = "P"};
	if($mylat =~ /P2\/n/) {$s_lattice = "P"};
	if($mylat =~ /P21\/c/) {$s_lattice = "P"};
	if($mylat =~ /P21\/a/) {$s_lattice = "P"};
	if($mylat =~ /P21\/b/) {$s_lattice = "P"};
	if($mylat =~ /P21\/b/) {$s_lattice = "P"};
	if($mylat =~ /P21\/c/) {$s_lattice = "P"};
	if($mylat =~ /P21\/a/) {$s_lattice = "P"};
	if($mylat =~ /P21\/n/) {$s_lattice = "P"};
	if($mylat =~ /P21\/n/) {$s_lattice = "P"};
	if($mylat =~ /P21\/n/) {$s_lattice = "P"};
	if($mylat =~ /C2\/c/) {$s_lattice = "CXY"};
	if($mylat =~ /A2\/a/) {$s_lattice = "CYZ"};
	if($mylat =~ /B2\/b/) {$s_lattice = "CXZ"};
	if($mylat =~ /B2\/b/) {$s_lattice = "CXZ"};
	if($mylat =~ /C2\/c/) {$s_lattice = "CXY"};
	if($mylat =~ /A2\/a/) {$s_lattice = "CYZ"};
	if($mylat =~ /P222/) {$s_lattice = "P"};
	if($mylat =~ /P2221/) {$s_lattice = "P"};
	if($mylat =~ /P2122/) {$s_lattice = "P"};
	if($mylat =~ /P2212/) {$s_lattice = "P"};
	if($mylat =~ /P21212/) {$s_lattice = "P"};
	if($mylat =~ /P22121/) {$s_lattice = "P"};
	if($mylat =~ /P21221/) {$s_lattice = "P"};
	if($mylat =~ /P212121/) {$s_lattice = "P"};
	if($mylat =~ /C2221/) {$s_lattice = "CXY"};
	if($mylat =~ /A2122/) {$s_lattice = "CYZ"};
	if($mylat =~ /B2212/) {$s_lattice = "CXZ"};
	if($mylat =~ /C222/) {$s_lattice = "CXY"};
	if($mylat =~ /A222/) {$s_lattice = "CYZ"};
	if($mylat =~ /B222/) {$s_lattice = "CXZ"};
	if($mylat =~ /F222/) {$s_lattice = "F"};
	if($mylat =~ /I222/) {$s_lattice = "B"};
	if($mylat =~ /I212121/) {$s_lattice = "B"};
	if($mylat =~ /Pmm2/) {$s_lattice = "P"};
	if($mylat =~ /P2mm/) {$s_lattice = "P"};
	if($mylat =~ /Pm2m/) {$s_lattice = "P"};
	if($mylat =~ /Pmc21/) {$s_lattice = "P"};
	if($mylat =~ /P21ma/) {$s_lattice = "P"};
	if($mylat =~ /Pb21m/) {$s_lattice = "P"};
	if($mylat =~ /Pcm21/) {$s_lattice = "P"};
	if($mylat =~ /P21am/) {$s_lattice = "P"};
	if($mylat =~ /Pm21b/) {$s_lattice = "P"};
	if($mylat =~ /Pcc2/) {$s_lattice = "P"};
	if($mylat =~ /P2aa/) {$s_lattice = "P"};
	if($mylat =~ /Pb2b/) {$s_lattice = "P"};
	if($mylat =~ /Pma2/) {$s_lattice = "P"};
	if($mylat =~ /P2mb/) {$s_lattice = "P"};
	if($mylat =~ /Pc2m/) {$s_lattice = "P"};
	if($mylat =~ /Pbm2/) {$s_lattice = "P"};
	if($mylat =~ /P2cm/) {$s_lattice = "P"};
	if($mylat =~ /Pm2a/) {$s_lattice = "P"};
	if($mylat =~ /Pca21/) {$s_lattice = "P"};
	if($mylat =~ /P21ab/) {$s_lattice = "P"};
	if($mylat =~ /Pc21b/) {$s_lattice = "P"};
	if($mylat =~ /Pbc21/) {$s_lattice = "P"};
	if($mylat =~ /P21ca/) {$s_lattice = "P"};
	if($mylat =~ /Pb21a/) {$s_lattice = "P"};
	if($mylat =~ /Pnc2/) {$s_lattice = "P"};
	if($mylat =~ /P2na/) {$s_lattice = "P"};
	if($mylat =~ /Pb2n/) {$s_lattice = "P"};
	if($mylat =~ /Pcn2/) {$s_lattice = "P"};
	if($mylat =~ /P2an/) {$s_lattice = "P"};
	if($mylat =~ /Pn2b/) {$s_lattice = "P"};
	if($mylat =~ /Pmn21/) {$s_lattice = "P"};
	if($mylat =~ /P21mn/) {$s_lattice = "P"};
	if($mylat =~ /Pn21m/) {$s_lattice = "P"};
	if($mylat =~ /Pnm21/) {$s_lattice = "P"};
	if($mylat =~ /P21nm/) {$s_lattice = "P"};
	if($mylat =~ /Pm21n/) {$s_lattice = "P"};
	if($mylat =~ /Pba2/) {$s_lattice = "P"};
	if($mylat =~ /P2cb/) {$s_lattice = "P"};
	if($mylat =~ /Pc2a/) {$s_lattice = "P"};
	if($mylat =~ /Pna21/) {$s_lattice = "P"};
	if($mylat =~ /P21nb/) {$s_lattice = "P"};
	if($mylat =~ /Pc21n/) {$s_lattice = "P"};
	if($mylat =~ /Pbn21/) {$s_lattice = "P"};
	if($mylat =~ /P21cn/) {$s_lattice = "P"};
	if($mylat =~ /Pn21a/) {$s_lattice = "P"};
	if($mylat =~ /Pnn2/) {$s_lattice = "P"};
	if($mylat =~ /P2nn/) {$s_lattice = "P"};
	if($mylat =~ /Pn2n/) {$s_lattice = "P"};
	if($mylat =~ /Cmm2/) {$s_lattice = "CXY"};
	if($mylat =~ /A2mm/) {$s_lattice = "CYZ"};
	if($mylat =~ /Bm2m/) {$s_lattice = "CXZ"};
	if($mylat =~ /Cmc21/) {$s_lattice = "CXY"};
	if($mylat =~ /A21ma/) {$s_lattice = "CYZ"};
	if($mylat =~ /Bb21m/) {$s_lattice = "CXZ"};
	if($mylat =~ /Ccm21/) {$s_lattice = "CXY"};
	if($mylat =~ /A21am/) {$s_lattice = "CYZ"};
	if($mylat =~ /Bm21b/) {$s_lattice = "CXZ"};
	if($mylat =~ /Ccc2/) {$s_lattice = "CXY"};
	if($mylat =~ /A2aa/) {$s_lattice = "CYZ"};
	if($mylat =~ /Bb2b/) {$s_lattice = "CXZ"};
	if($mylat =~ /Amm2/) {$s_lattice = "CYZ"};
	if($mylat =~ /B2mm/) {$s_lattice = "CXZ"};
	if($mylat =~ /Cm2m/) {$s_lattice = "CXY"};
	if($mylat =~ /Abm2/) {$s_lattice = "CYZ"};
	if($mylat =~ /B2cm/) {$s_lattice = "CXZ"};
	if($mylat =~ /Cm2a/) {$s_lattice = "CXY"};
	if($mylat =~ /Bma2/) {$s_lattice = "CXZ"};
	if($mylat =~ /C2mb/) {$s_lattice = "CXY"};
	if($mylat =~ /Ac2m/) {$s_lattice = "CYZ"};
	if($mylat =~ /Ama2/) {$s_lattice = "CYZ"};
	if($mylat =~ /B2mb/) {$s_lattice = "CXZ"};
	if($mylat =~ /Cc2m/) {$s_lattice = "CXY"};
	if($mylat =~ /Bbm2/) {$s_lattice = "CXZ"};
	if($mylat =~ /C2cm/) {$s_lattice = "CXY"};
	if($mylat =~ /Am2a/) {$s_lattice = "CYZ"};
	if($mylat =~ /Aba2/) {$s_lattice = "CYZ"};
	if($mylat =~ /B2cb/) {$s_lattice = "CXZ"};
	if($mylat =~ /Cc2a/) {$s_lattice = "CXY"};
	if($mylat =~ /Bba2/) {$s_lattice = "CXZ"};
	if($mylat =~ /C2cb/) {$s_lattice = "CXY"};
	if($mylat =~ /Ac2a/) {$s_lattice = "CYZ"};
	if($mylat =~ /Fmm2/) {$s_lattice = "F"};
	if($mylat =~ /F2mm/) {$s_lattice = "F"};
	if($mylat =~ /Fm2m/) {$s_lattice = "F"};
	if($mylat =~ /Fdd2/) {$s_lattice = "F"};
	if($mylat =~ /F2dd/) {$s_lattice = "F"};
	if($mylat =~ /Fd2d/) {$s_lattice = "F"};
	if($mylat =~ /Imm2/) {$s_lattice = "B"};
	if($mylat =~ /I2mm/) {$s_lattice = "B"};
	if($mylat =~ /Im2m/) {$s_lattice = "B"};
	if($mylat =~ /Iba2/) {$s_lattice = "B"};
	if($mylat =~ /I2cb/) {$s_lattice = "B"};
	if($mylat =~ /Ic2a/) {$s_lattice = "B"};
	if($mylat =~ /Ima2/) {$s_lattice = "B"};
	if($mylat =~ /I2mb/) {$s_lattice = "B"};
	if($mylat =~ /Ic2m/) {$s_lattice = "B"};
	if($mylat =~ /Ibm2/) {$s_lattice = "B"};
	if($mylat =~ /I2cm/) {$s_lattice = "B"};
	if($mylat =~ /Im2a/) {$s_lattice = "B"};
	if($mylat =~ /Pmmm/) {$s_lattice = "P"};
	if($mylat =~ /Pnnn/) {$s_lattice = "P"};
	if($mylat =~ /Pccm/) {$s_lattice = "P"};
	if($mylat =~ /Pmaa/) {$s_lattice = "P"};
	if($mylat =~ /Pbmb/) {$s_lattice = "P"};
	if($mylat =~ /Pban/) {$s_lattice = "P"};
	if($mylat =~ /Pncb/) {$s_lattice = "P"};
	if($mylat =~ /Pcna/) {$s_lattice = "P"};
	if($mylat =~ /Pmma/) {$s_lattice = "P"};
	if($mylat =~ /Pbmm/) {$s_lattice = "P"};
	if($mylat =~ /Pmcm/) {$s_lattice = "P"};
	if($mylat =~ /Pmam/) {$s_lattice = "P"};
	if($mylat =~ /Pmmb/) {$s_lattice = "P"};
	if($mylat =~ /Pcmm/) {$s_lattice = "P"};
	if($mylat =~ /Pnna/) {$s_lattice = "P"};
	if($mylat =~ /Pbnn/) {$s_lattice = "P"};
	if($mylat =~ /Pncn/) {$s_lattice = "P"};
	if($mylat =~ /Pnan/) {$s_lattice = "P"};
	if($mylat =~ /Pnnb/) {$s_lattice = "P"};
	if($mylat =~ /Pcnn/) {$s_lattice = "P"};
	if($mylat =~ /Pmna/) {$s_lattice = "P"};
	if($mylat =~ /Pbmn/) {$s_lattice = "P"};
	if($mylat =~ /Pncm/) {$s_lattice = "P"};
	if($mylat =~ /Pman/) {$s_lattice = "P"};
	if($mylat =~ /Pnmb/) {$s_lattice = "P"};
	if($mylat =~ /Pcnm/) {$s_lattice = "P"};
	if($mylat =~ /Pcca/) {$s_lattice = "P"};
	if($mylat =~ /Pbaa/) {$s_lattice = "P"};
	if($mylat =~ /Pbcb/) {$s_lattice = "P"};
	if($mylat =~ /Pbab/) {$s_lattice = "P"};
	if($mylat =~ /Pccb/) {$s_lattice = "P"};
	if($mylat =~ /Pcaa/) {$s_lattice = "P"};
	if($mylat =~ /Pbam/) {$s_lattice = "P"};
	if($mylat =~ /Pmcb/) {$s_lattice = "P"};
	if($mylat =~ /Pcma/) {$s_lattice = "P"};
	if($mylat =~ /Pccn/) {$s_lattice = "P"};
	if($mylat =~ /Pnaa/) {$s_lattice = "P"};
	if($mylat =~ /Pbnb/) {$s_lattice = "P"};
	if($mylat =~ /Pbcm/) {$s_lattice = "P"};
	if($mylat =~ /Pmca/) {$s_lattice = "P"};
	if($mylat =~ /Pbma/) {$s_lattice = "P"};
	if($mylat =~ /Pcmb/) {$s_lattice = "P"};
	if($mylat =~ /Pcam/) {$s_lattice = "P"};
	if($mylat =~ /Pmab/) {$s_lattice = "P"};
	if($mylat =~ /Pnnm/) {$s_lattice = "P"};
	if($mylat =~ /Pmnn/) {$s_lattice = "P"};
	if($mylat =~ /Pnmn/) {$s_lattice = "P"};
	if($mylat =~ /Pmmn/) {$s_lattice = "P"};
	if($mylat =~ /Pnmm/) {$s_lattice = "P"};
	if($mylat =~ /Pmnm/) {$s_lattice = "P"};
	if($mylat =~ /Pbcn/) {$s_lattice = "P"};
	if($mylat =~ /Pnca/) {$s_lattice = "P"};
	if($mylat =~ /Pbna/) {$s_lattice = "P"};
	if($mylat =~ /Pcnb/) {$s_lattice = "P"};
	if($mylat =~ /Pcan/) {$s_lattice = "P"};
	if($mylat =~ /Pnab/) {$s_lattice = "P"};
	if($mylat =~ /Pbca/) {$s_lattice = "P"};
	if($mylat =~ /Pcab/) {$s_lattice = "P"};
	if($mylat =~ /Pnma/) {$s_lattice = "P"};
	if($mylat =~ /Pbnm/) {$s_lattice = "P"};
	if($mylat =~ /Pmcn/) {$s_lattice = "P"};
	if($mylat =~ /Pnam/) {$s_lattice = "P"};
	if($mylat =~ /Pmnb/) {$s_lattice = "P"};
	if($mylat =~ /Pcmn/) {$s_lattice = "P"};
	if($mylat =~ /Cmcm/) {$s_lattice = "CXY"};
	if($mylat =~ /Amma/) {$s_lattice = "CYZ"};
	if($mylat =~ /Bbmm/) {$s_lattice = "CXZ"};
	if($mylat =~ /Bmmb/) {$s_lattice = "CXZ"};
	if($mylat =~ /Ccmm/) {$s_lattice = "CXY"};
	if($mylat =~ /Amam/) {$s_lattice = "CYZ"};
	if($mylat =~ /Cmca/) {$s_lattice = "CXY"};
	if($mylat =~ /Abma/) {$s_lattice = "CYZ"};
	if($mylat =~ /Bbcm/) {$s_lattice = "CXZ"};
	if($mylat =~ /Bmab/) {$s_lattice = "CXZ"};
	if($mylat =~ /Ccmb/) {$s_lattice = "CXY"};
	if($mylat =~ /Acam/) {$s_lattice = "CYZ"};
	if($mylat =~ /Cmmm/) {$s_lattice = "CXY"};
	if($mylat =~ /Ammm/) {$s_lattice = "CYZ"};
	if($mylat =~ /Bmmm/) {$s_lattice = "CXZ"};
	if($mylat =~ /Cccm/) {$s_lattice = "CXY"};
	if($mylat =~ /Amaa/) {$s_lattice = "CYZ"};
	if($mylat =~ /Bbmb/) {$s_lattice = "CXZ"};
	if($mylat =~ /Cmma/) {$s_lattice = "CXY"};
	if($mylat =~ /Abmm/) {$s_lattice = "CYZ"};
	if($mylat =~ /Bmcm/) {$s_lattice = "CXZ"};
	if($mylat =~ /Bmam/) {$s_lattice = "CXZ"};
	if($mylat =~ /Cmmb/) {$s_lattice = "CXY"};
	if($mylat =~ /Acmm/) {$s_lattice = "CYZ"};
	if($mylat =~ /Ccca/) {$s_lattice = "CXY"};
	if($mylat =~ /Abaa/) {$s_lattice = "CYZ"};
	if($mylat =~ /Bbcb/) {$s_lattice = "CXZ"};
	if($mylat =~ /Bbab/) {$s_lattice = "CXZ"};
	if($mylat =~ /Cccb/) {$s_lattice = "CXY"};
	if($mylat =~ /Acaa/) {$s_lattice = "CYZ"};
	if($mylat =~ /Fmmm/) {$s_lattice = "F"};
	if($mylat =~ /Fddd/) {$s_lattice = "F"};
	if($mylat =~ /Immm/) {$s_lattice = "B"};
	if($mylat =~ /Ibam/) {$s_lattice = "B"};
	if($mylat =~ /Imcb/) {$s_lattice = "B"};
	if($mylat =~ /Icma/) {$s_lattice = "B"};
	if($mylat =~ /Ibca/) {$s_lattice = "B"};
	if($mylat =~ /Icab/) {$s_lattice = "B"};
	if($mylat =~ /Imma/) {$s_lattice = "B"};
	if($mylat =~ /Ibmm/) {$s_lattice = "B"};
	if($mylat =~ /Imcm/) {$s_lattice = "B"};
	if($mylat =~ /Imam/) {$s_lattice = "B"};
	if($mylat =~ /Immb/) {$s_lattice = "B"};
	if($mylat =~ /Icmm/) {$s_lattice = "B"};
	if($mylat =~ /P4/) {$s_lattice = "P"};
	if($mylat =~ /P41/) {$s_lattice = "P"};
	if($mylat =~ /P42/) {$s_lattice = "P"};
	if($mylat =~ /P43/) {$s_lattice = "P"};
	if($mylat =~ /I4/) {$s_lattice = "B"};
	if($mylat =~ /I41/) {$s_lattice = "B"};
	if($mylat =~ /P-4/) {$s_lattice = "P"};
	if($mylat =~ /I-4/) {$s_lattice = "B"};
	if($mylat =~ /P4\/m/) {$s_lattice = "P"};
	if($mylat =~ /P42\/m/) {$s_lattice = "P"};
	if($mylat =~ /P4\/n/) {$s_lattice = "P"};
	if($mylat =~ /P42\/n/) {$s_lattice = "P"};
	if($mylat =~ /I4\/m/) {$s_lattice = "B"};
	if($mylat =~ /I41\/a/) {$s_lattice = "B"};
	if($mylat =~ /P422/) {$s_lattice = "P"};
	if($mylat =~ /P4212/) {$s_lattice = "P"};
	if($mylat =~ /P4122/) {$s_lattice = "P"};
	if($mylat =~ /P41212/) {$s_lattice = "P"};
	if($mylat =~ /P4222/) {$s_lattice = "P"};
	if($mylat =~ /P42212/) {$s_lattice = "P"};
	if($mylat =~ /P4322/) {$s_lattice = "P"};
	if($mylat =~ /P43212/) {$s_lattice = "P"};
	if($mylat =~ /I422/) {$s_lattice = "B"};
	if($mylat =~ /I4122/) {$s_lattice = "B"};
	if($mylat =~ /P4mm/) {$s_lattice = "P"};
	if($mylat =~ /P4bm/) {$s_lattice = "P"};
	if($mylat =~ /P42cm/) {$s_lattice = "P"};
	if($mylat =~ /P42nm/) {$s_lattice = "P"};
	if($mylat =~ /P4cc/) {$s_lattice = "P"};
	if($mylat =~ /P4nc/) {$s_lattice = "P"};
	if($mylat =~ /P42mc/) {$s_lattice = "P"};
	if($mylat =~ /P42bc/) {$s_lattice = "P"};
	if($mylat =~ /I4mm/) {$s_lattice = "B"};
	if($mylat =~ /I4cm/) {$s_lattice = "B"};
	if($mylat =~ /I41md/) {$s_lattice = "B"};
	if($mylat =~ /I41cd/) {$s_lattice = "B"};
	if($mylat =~ /P-42m/) {$s_lattice = "P"};
	if($mylat =~ /P-42c/) {$s_lattice = "P"};
	if($mylat =~ /P-421m/) {$s_lattice = "P"};
	if($mylat =~ /P-421c/) {$s_lattice = "P"};
	if($mylat =~ /P-4m2/) {$s_lattice = "P"};
	if($mylat =~ /P-4c2/) {$s_lattice = "P"};
	if($mylat =~ /P-4b2/) {$s_lattice = "P"};
	if($mylat =~ /P-4n2/) {$s_lattice = "P"};
	if($mylat =~ /I-4m2/) {$s_lattice = "B"};
	if($mylat =~ /I-4c2/) {$s_lattice = "B"};
	if($mylat =~ /I-42m/) {$s_lattice = "B"};
	if($mylat =~ /I-42d/) {$s_lattice = "B"};
	if($mylat =~ /P4\/mmm/) {$s_lattice = "P"};
	if($mylat =~ /P4\/mcc/) {$s_lattice = "P"};
	if($mylat =~ /P4\/nbm/) {$s_lattice = "P"};
	if($mylat =~ /P4\/nnc/) {$s_lattice = "P"};
	if($mylat =~ /P4\/mbm/) {$s_lattice = "P"};
	if($mylat =~ /P4\/mnc/) {$s_lattice = "P"};
	if($mylat =~ /P4\/nmm/) {$s_lattice = "P"};
	if($mylat =~ /P4\/ncc/) {$s_lattice = "P"};
	if($mylat =~ /P42\/mmc/) {$s_lattice = "P"};
	if($mylat =~ /P42\/mcm/) {$s_lattice = "P"};
	if($mylat =~ /P42\/nbc/) {$s_lattice = "P"};
	if($mylat =~ /P42\/nnm/) {$s_lattice = "P"};
	if($mylat =~ /P42\/mbc/) {$s_lattice = "P"};
	if($mylat =~ /P42\/mnm/) {$s_lattice = "P"};
	if($mylat =~ /P42\/nmc/) {$s_lattice = "P"};
	if($mylat =~ /P42\/ncm/) {$s_lattice = "P"};
	if($mylat =~ /I4\/mmm/) {$s_lattice = "B"};
	if($mylat =~ /I4\/mcm/) {$s_lattice = "B"};
	if($mylat =~ /I41\/amd/) {$s_lattice = "B"};
	if($mylat =~ /I41\/acd/) {$s_lattice = "B"};
	if($mylat =~ /P3/) {$s_lattice = "H"};
	if($mylat =~ /P31/) {$s_lattice = "H"};
	if($mylat =~ /P32/) {$s_lattice = "H"};
	if($mylat =~ /R3/) {$s_lattice = "R"};
	if($mylat =~ /P-3/) {$s_lattice = "H"};
	if($mylat =~ /R-3/) {$s_lattice = "R"};
	if($mylat =~ /P312/) {$s_lattice = "H"};
	if($mylat =~ /P321/) {$s_lattice = "H"};
	if($mylat =~ /P3112/) {$s_lattice = "H"};
	if($mylat =~ /P3121/) {$s_lattice = "H"};
	if($mylat =~ /P3212/) {$s_lattice = "H"};
	if($mylat =~ /P3221/) {$s_lattice = "H"};
	if($mylat =~ /R32/) {$s_lattice = "R"};
	if($mylat =~ /P3m1/) {$s_lattice = "H"};
	if($mylat =~ /P31m/) {$s_lattice = "H"};
	if($mylat =~ /P3c1/) {$s_lattice = "H"};
	if($mylat =~ /P31c/) {$s_lattice = "H"};
	if($mylat =~ /R3m/) {$s_lattice = "R"};
	if($mylat =~ /R3c/) {$s_lattice = "R"};
	if($mylat =~ /P-31m/) {$s_lattice = "H"};
	if($mylat =~ /P-31c/) {$s_lattice = "H"};
	if($mylat =~ /P-3m1/) {$s_lattice = "H"};
	if($mylat =~ /P-3c1/) {$s_lattice = "H"};
	if($mylat =~ /R-3m/) {$s_lattice = "R"};
	if($mylat =~ /R-3c/) {$s_lattice = "R"};
	if($mylat =~ /P6/) {$s_lattice = "H"};
	if($mylat =~ /P61/) {$s_lattice = "H"};
	if($mylat =~ /P65/) {$s_lattice = "H"};
	if($mylat =~ /P62/) {$s_lattice = "H"};
	if($mylat =~ /P64/) {$s_lattice = "H"};
	if($mylat =~ /P63/) {$s_lattice = "H"};
	if($mylat =~ /P-6/) {$s_lattice = "H"};
	if($mylat =~ /P6\/m/) {$s_lattice = "H"};
	if($mylat =~ /P63\/m/) {$s_lattice = "H"};
	if($mylat =~ /P622/) {$s_lattice = "H"};
	if($mylat =~ /P6122/) {$s_lattice = "H"};
	if($mylat =~ /P6522/) {$s_lattice = "H"};
	if($mylat =~ /P6222/) {$s_lattice = "H"};
	if($mylat =~ /P6422/) {$s_lattice = "H"};
	if($mylat =~ /P6322/) {$s_lattice = "H"};
	if($mylat =~ /P6mm/) {$s_lattice = "H"};
	if($mylat =~ /P6cc/) {$s_lattice = "H"};
	if($mylat =~ /P63cm/) {$s_lattice = "H"};
	if($mylat =~ /P63mc/) {$s_lattice = "H"};
	if($mylat =~ /P-6m2/) {$s_lattice = "H"};
	if($mylat =~ /P-6c2/) {$s_lattice = "H"};
	if($mylat =~ /P-62m/) {$s_lattice = "H"};
	if($mylat =~ /P-62c/) {$s_lattice = "H"};
	if($mylat =~ /P6\/mmm/) {$s_lattice = "H"};
	if($mylat =~ /P6\/mcc/) {$s_lattice = "H"};
	if($mylat =~ /P63\/mcm/) {$s_lattice = "H"};
	if($mylat =~ /P63\/mmc/) {$s_lattice = "H"};
	if($mylat =~ /P23/) {$s_lattice = "P"};
	if($mylat =~ /F23/) {$s_lattice = "F"};
	if($mylat =~ /I23/) {$s_lattice = "B"};
	if($mylat =~ /P213/) {$s_lattice = "P"};
	if($mylat =~ /I213/) {$s_lattice = "B"};
	if($mylat =~ /Pm-3/) {$s_lattice = "P"};
	if($mylat =~ /Pn-3/) {$s_lattice = "P"};
	if($mylat =~ /Fm-3/) {$s_lattice = "F"};
	if($mylat =~ /Fd-3/) {$s_lattice = "F"};
	if($mylat =~ /Im-3/) {$s_lattice = "B"};
	if($mylat =~ /Pa-3/) {$s_lattice = "P"};
	if($mylat =~ /Ia-3/) {$s_lattice = "B"};
	if($mylat =~ /P432/) {$s_lattice = "P"};
	if($mylat =~ /P4232/) {$s_lattice = "P"};
	if($mylat =~ /F432/) {$s_lattice = "F"};
	if($mylat =~ /F4132/) {$s_lattice = "F"};
	if($mylat =~ /I432/) {$s_lattice = "B"};
	if($mylat =~ /P4332/) {$s_lattice = "P"};
	if($mylat =~ /P4132/) {$s_lattice = "P"};
	if($mylat =~ /I4132/) {$s_lattice = "B"};
	if($mylat =~ /P-43m/) {$s_lattice = "P"};
	if($mylat =~ /F-43m/) {$s_lattice = "F"};
	if($mylat =~ /I-43m/) {$s_lattice = "B"};
	if($mylat =~ /P-43n/) {$s_lattice = "P"};
	if($mylat =~ /F-43c/) {$s_lattice = "F"};
	if($mylat =~ /I-43d/) {$s_lattice = "B"};
	if($mylat =~ /Pm-3m/) {$s_lattice = "P"};
	if($mylat =~ /Pn-3n/) {$s_lattice = "P"};
	if($mylat =~ /Pm-3n/) {$s_lattice = "P"};
	if($mylat =~ /Pn-3m/) {$s_lattice = "P"};
	if($mylat =~ /Fm-3m/) {$s_lattice = "F"};
	if($mylat =~ /Fm-3c/) {$s_lattice = "F"};
	if($mylat =~ /Fd-3m/) {$s_lattice = "F"};
	if($mylat =~ /Fd-3c/) {$s_lattice = "F"};
	if($mylat =~ /Im-3m/) {$s_lattice = "B"};
	if($mylat =~ /Ia-3d/) {$s_lattice = "B"};

	my $mytest = 0;
	foreach $v (@lattic) {
		if ($mylat =~ /^$v/) { $mytest = 1;}
	}
	if ($mytest) {
		$s_spacegr = "";
	} else {
		$s_spacegr = "$mylat";
	}
}

1;
