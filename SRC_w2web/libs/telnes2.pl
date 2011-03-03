$t_output = 0;
$t_atoms  = 0;
$t_det    = 0;
$t_modus  = 0;
$t_srule    = 0;
$t_lrule    = 0;
$t_split    = 0;
$t_branch   = 0;
$t_init     = 0;
$t_nrel     = 0;
$t_orient   = 0;
$t_qgrid    = 0;
$t_noheader = 0;
$t_xqtl     = 0;

sub InnesRead {
    my $i;
    my $j;
#    $debug = 1;

    $file = "$DIR/$CASE.innes";
    
    if ( -e $file ) {
	# good, we have a innes-file
    } else {
	# we need a template
	&RequiredFile("innes");
	exit;
    }


    open (INNES, $file);
    # read title
    $_ = <INNES>;
    chomp;
    $t_title = $_;

    $t_title =~ s/^ *//;
    $t_title =~ s/ *$//;

    &w("title: $t_title");

    # read atom
    $_ = <INNES>;
    chomp;
    @dum = split;
    $t_atom = @dum[0];
    &w("atom:  $t_atom");

    # read N + L
    $_ = <INNES>;
    chomp;
    @dum = split;
    $t_n=@dum[0];
    $t_l=@dum[1];
    &w("n= $t_n l=$t_l");

    # read E_loss
    $_ = <INNES>;
    chomp;
    @dum = split;
    $t_ene = @dum[0];
    &w("E_loss= $t_ene");

    # read beam energy
    $_ = <INNES>;
    chomp;
    @dum = split;
    $t_beam = @dum[0];
    &w("E_beam= $t_beam");

    # energy grid
    $_ = <INNES>;
    chomp;
    @dum = split;
    $t_egrid1= @dum[0];
    $t_egrid2= @dum[1];
    $t_egrid3= @dum[2];
    &w("E_grid:  $t_egrid1 to $t_egrid2 in steps of $t_egrid3");

    # read coll. and conv. semiangles
    $_ = <INNES>;
    chomp;
    @dum = split;
    $t_coll= @dum[0];
    $t_conv= @dum[1];
    &w("coll= $t_coll, conv=$t_conv");

    # read NR and NT
    $_ = <INNES>;
    chomp;
    @dum = split;
    $t_nr= @dum[0];
    $t_nt= @dum[1];
    &w("NR= $t_nr, NT=$t_nt");

    # read spectr. broadening
    $_ = <INNES>;
    chomp;
    @dum = split;
    $t_spec= @dum[0];
    &w("spectr.broad=$t_spec");


    # this was obligatory input, now we parse for optional things

    while(<INNES>) {
      chop;
      # only first 4 letters are significant !
      

      if ( /^OUTP/ ) {
	  &w("found: OUTPUT");
	  $t_output=1;
	  $_ = <INNES>;
	  chomp;
	  @dum = split;
	  $t_output1= @dum[0];
	  &w("-> $t_output1");
      }
      if ( /^ATOM/ ) {
	  &w("found: ATOMS");
	  $t_atoms=1;
	  $_ = <INNES>;
	  chomp;
	  @dum = split;
	  $t_atoms1= @dum[0];
	  $t_atoms2= @dum[1];
	  &w("-> $t_atoms1, $t_atoms2");
      }
      if ( /^DETE/ ) {
	  &w("found: DETECTOR POSITION");
	  $t_det=1;
	  $_ = <INNES>;
	  chomp;
	  @dum = split;
	  $t_det1= @dum[0];
	  $t_det2= @dum[1];
	  &w("-> $t_det1, $t_det2");
      }
      if ( /^MODU/ ) {
	  &w("found: MODUS");
	  $t_modus=1;
	  $_ = <INNES>;
	  chomp;
	  @dum = split;
	  $t_modus1= @dum[0];
	  &w("-> $t_modus1");
      }
      if ( /^SELE/ ) {
	  &w("found: SELECTION RULE");
	  $t_srule=1;
	  $_ = <INNES>;
	  chomp;
	  @dum = split;
	  $t_srule1= @dum[0];
	  &w("-> $t_srule1");
      }
      if ( /^LSEL/ ) {
	  &w("found: LSELECTION RULE");
	  $t_lrule=1;
	  $_ = <INNES>;
	  chomp;
	  @dum = split;
	  $t_lrule1= @dum[0];
	  &w("-> $t_lrule1");
      }
      if ( /^SPLI/ ) {
	  &w("found: SPLIT");
	  $t_split=1;
	  $_ = <INNES>;
	  chomp;
	  @dum = split;
	  $t_split1= @dum[0];
	  &w("-> $t_split1");
      }
      if ( /^BRAN/ ) {
	  &w("found: BRANCHING RATIO");
	  $t_branch=1;
	  $_ = <INNES>;
	  chomp;
	  @dum = split;
	  $t_branch1= @dum[0];
	  &w("-> $t_branch1");
      }
      if ( /^INIT/ ) {
	  &w("found: INITIALIZATION");
	  $t_init=1;
	  $_ = <INNES>;
	  chomp;
	  @dum = split;
	  $t_init1= @dum[0];
	  $t_init2= @dum[1];
	  # one more line
	  $_ = <INNES>;
	  
	  chomp;
	  @dum = split;
	  $t_init3= @dum[0];
	  $t_init4= @dum[1];

	  &w("-> $t_init1 $t_init2 $t_init3 $t_init4");
      }
      if ( /^NONR/ ) {
	  &w("found: NONRELATIVISTIC");
	  $t_nrel=1;
	  &w("-> (no value)");
      }

      if ( /^ORIE/ ) {
	  &w("found: ORIENTATION SENSITIVE");
	  $t_orient=1;
	  $_ = <INNES>;
	  chomp;
	  @dum = split;
	  $t_orient1= @dum[0];
	  $t_orient2= @dum[1];
	  $t_orient3= @dum[2];

	  &w("-> $t_orient1, $t_orient2, $t_orient3");
      }
      if ( /^QGRI/ ) {
	  &w("found: QGRID");
	  $t_qgrid=1;
	  $_ = <INNES>;
	  chomp;
	  @dum = split;
	  $t_qgrid1 = @dum[0];
	  &w("-> $t_qgrid1");
	  if($t_qgrid1 =~ /U/) {
	      #we are done
	  } else {
	      #extra parameter here
	      $_ = <INNES>;
	      chomp;
	      @dum = split;
	      $t_qgrid2 = @dum[0];
	      &w("->$t_qgrid2");
	  }
      }
      if ( /^NOHE/ ) {
	  &w("found: NOHEADERS");
	  $t_noheader=1;
	  &w("-> (no value)");
      }
      if ( /^XQTL/ ) {
	  &w("found: XQTL");
	  $t_xqtl=1;
	  &w("-> (no value)");
      }
  }
}


sub w {
  my ($myval) = @_;
  $OUT .= "$myval<br>" if $debug;
}

    
sub InnesWrite {
    #$debug = 1;
    
    open (INNES, ">$DIR/$CASE.innes");
    
    print INNES sprintf "%-60s\n", $t_title;
    print INNES sprintf "%i\n", $t_atom;
    print INNES sprintf "%i %i\n", $t_n, $t_l;
    print INNES sprintf "%.2f\n", $t_ene;
    print INNES sprintf "%i\n", $t_beam;
    print INNES sprintf "%.4f %.4f %.4f\n",$t_egrid1,$t_egrid2,$t_egrid3;
    print INNES sprintf "%.2f %.2f\n", $t_coll, $t_conv;
    print INNES sprintf "%i %i\n", $t_nr, $t_nt;
    print INNES sprintf "%.2f\n", $t_spec;
    # this was obligatory

    if ($t_output) {
	print INNES "OUTPUT\n";
	print INNES sprintf "%i\n",$t_output1;
    }
    if ($t_atoms) {
	print INNES "ATOMS\n";
	print INNES sprintf "%i %i\n",$t_atoms1,$t_atoms2;
    }
    if ($t_det) {
	print INNES "DETECTOR POSITION\n";
	print INNES sprintf "%.3f %.3f\n",$t_det1,$t_det2;
    }
    if ($t_modus) {
	print INNES "MODUS\n";
	print INNES sprintf "%s\n",$t_modus1;
    }
    if ($t_srule) {
	print INNES "SELECTION RULE\n";
	print INNES sprintf "%s\n",$t_srule1;
    }
    if ($t_lrule) {
	print INNES "LSELECTION RULE\n";
	print INNES sprintf "%s\n",$t_lrule1;
    }

    if ($t_split) {
	print INNES "SPLIT\n";
	print INNES sprintf "%.3f\n",$t_split1;
    }
    if ($t_branch) {
	print INNES "BRANCHING RATIO\n";
	print INNES sprintf "%.4f\n",$t_branch1;
    }
    if ($t_init) {
	print INNES "INITIALIZATION\n";
	print INNES sprintf "%s %s\n%s %s\n",$t_init1,$t_init2,$t_init3,$t_init4;
    }
    if ($t_nrel) {
	print INNES "NONRELATIVISTIC\n";
    }
   if ($t_orient) {
	print INNES "ORIENTATION SENSITIVE\n";
	print INNES sprintf "%.2f %.2f %.2f\n",$t_orient1,$t_orient2,$t_orient3;
    }
   if ($t_qgrid) {
	print INNES "QGRID\n";
	print INNES sprintf "%s\n",$t_qgrid1;
	if ($t_qgrid =~ /U/) {
	    #done
	} else {
	    print INNES sprintf "%.5f\n",$t_qgrid2;
	}
    }
    if ($t_noheader) {
	print INNES "NOHEADERS\n";
    }
    if ($t_xqtl) {
	print INNES "XQTL\n";
    }

    print INNES "END\n";
    &w("innes_i written");
}

1;
