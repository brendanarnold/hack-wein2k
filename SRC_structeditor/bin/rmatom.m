function sr=rmatom(s,ind)

#   usage:  sr=rmatom(s,ind)
#
#       removes an atom from the structure s and returns all as a structure sr
#
#        sr      output structure
#        s       matlab structure used to store input structural information
#        ind     a index of a removed atom (position in the list),
#
#  example:
#  
#         sr=rmatom(s,3)        
#

      for i=ind:s.nat-1
          s.aname(i,:)=s.aname(i+1,:);
          s.jrj(i)=s.jrj(i+1);
          s.pos(i,:)=s.pos(i+1,:);
          s.r0(i)=s.r0(i+1);
          s.rmt(i)=s.rmt(i+1);
          s.zz(i)=s.zz(i+1); 
     end 

     sr.nat=s.nat-1;
     
     sr.aname(1:sr.nat,:)=s.aname(1:sr.nat,:);
     sr.jrj(1:sr.nat,1)=s.jrj(1:sr.nat);
     sr.pos(1:sr.nat,:)=s.pos(1:sr.nat,:);
     sr.r0(1:sr.nat,1)=s.r0(1:sr.nat);
     sr.rmt(1:sr.nat,1)=s.rmt(1:sr.nat);
     sr.zz(1:sr.nat,1)=s.zz(1:sr.nat);  

     sr.lattic=s.lattic;
     sr.a=s.a;
     sr.alpha=s.alpha;
     sr.lat2car=s.lat2car;
     sr.brlat=s.brlat;
     
end
