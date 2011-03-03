  function sr=rescale_c(s,znew)

#   usage:  sr=rescale_c(s,znew)
#
#        rescales atomic postions according to new c lattice parameter
#        for "vacuum" at z=0.5   (use rescale_c_2 for vacuum at z=0)
#
#        s      input structure
#        znew   new c lattice parameter (becomes a(3))
#
#   example: s2=rescale_c(s,30.0)
#

     sr=s;
     for i=1:s.nat
         if (s.pos(i,3) < 0.5) 
            sr.pos(i,3)=sr.pos(i,3)*sr.a(3)/znew;
          else
            sr.pos(i,3)=1-((1-sr.pos(i,3))*sr.a(3)/znew); 
          endif
     end
     sr.a(3)=znew;
     tm=[1,0,0;
         0,1,0;
         0,0,znew/s.a(3)];
     sr.lat2car= tm*sr.lat2car;
     sr.brlat= tm*sr.brlat;

 end
