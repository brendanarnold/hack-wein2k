 function  sr=copyatom2(s,s1,ind1,ind2,pos)

#    usage:  sr=copyatom2(s,s1,ind1,ind2,pos)
#
#        creates a copy of atom ind1 and puts it in the structure with
#        index ind2
#
#        s     input structure containing the atom
#        s1    input structure intowhich atom is copied
#        ind1  index of oryginal atom
#        ind2  index of new atom (if = 0, then adds it at the end)
#        pos   position of new atom
#       
#    example: sr=copyatom(s,1,0,[0.1 0.1 0.1}) 
#
#
#
    aname(:)=s.aname(ind1,:);
    rmt=s.rmt(ind1);
    jrj=s.jrj(ind1);
    r0=s.r0(ind1); 
    sr=addatom(s1,aname,pos,ind2);
    if (ind2 == 0)
       ind2=sr.nat;
    endif
    sr.rmt(ind2,1)=rmt;  
    sr.jrj(ind2,1)=jrj;
    sr.r0(ind2,1)=r0;

 end
