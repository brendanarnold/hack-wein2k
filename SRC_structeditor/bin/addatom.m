function sr=addatom(s,atom,pos,ind)

#  usage:  sr=addatom(s,atom,pos,ind) 
#
#        adds an atom to the structure s and returns all as structure sr
#
#        sr      output structure
#        s       matlab structure used to store structural information
#        atom    identyfies atom, can be string variable containing
#                atomic symbol or scalar with atomic number
#        pos     position in unit cell
#        ind     a index of new atom (position in the list),
#                if 0 then atom is added to the end of the list
#
#  example:
#  
#         sr=addatom(s,"N",[0,0,0.1],3)        
#

    if (isnumeric(atom)) 
       elmzz=atom;
       elem="            "
       elem=getaname(elmzz);      
    else
       elem=[atom," "];
       elmzz=getazz(elem);  
    endif
 
    if (ind==0)
       ind=s.nat+1 ;
    endif
 
    if ((ind>0)&(ind<=s.nat)) 

       for i=s.nat:-1:ind

           s.aname(i+1,:)=s.aname(i,:);
           s.jrj(i+1,1)=s.jrj(i,1);
           s.pos(i+1,:)=s.pos(i,:);
           s.r0(i+1,1)=s.r0(i,1);
           s.rmt(i+1,1)=s.rmt(i,1);
           s.zz(i+1,1)=s.zz(i,1);
       end	        
    endif 

    elem1=[elem,"        "];
    s.aname(ind,1:10)=elem1(1:10);
    s.jrj(ind,1)=781;
    s.pos(ind,1:3)=pos(1:3);
    s.r0(ind,1)=getar0(elmzz);
    s.rmt(ind,1)=2.0;
    s.zz(ind,1)=elmzz;
    s.nat=s.nat+1;
         
    sr=s;

end	
