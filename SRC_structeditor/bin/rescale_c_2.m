  function sr=resclale_c_2(s,znew)

     sr=s;
     for i=1:s.nat
         if (s.pos(i,3) < 1.5) 
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
