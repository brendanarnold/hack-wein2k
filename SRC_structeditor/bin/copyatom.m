 function  sr=copyatom(s,ind1,ind2,pos)

#    usage:  sr=copyatom(s,ind1,ind2,pos)
#
#        creates a copy of atom ind1 and puts it in the structure with
#        index ind2
#
#        s     input structure
#        ind1  index of oryginal atom
#        ind2  index of new atom (if = 0, then adds it at the end)
#        pos   position of new atom
#       
#    example: sr=copyatom(s,1,0,[0.1 0.1 0.1}) 
#
#
#
    aname(:)=s.aname(ind1,:);
    sr=addatom(s,aname,pos,ind2);
    if (ind2 == 0)
       ind2=sr.nat;
    endif
    sr.rmt(ind2,1)=s.rmt(ind1,1);  
    sr.jrj(ind2,1)=s.jrj(ind1,1);
    sr.r0(ind2,1)=s.r0(ind1,1);

 end
