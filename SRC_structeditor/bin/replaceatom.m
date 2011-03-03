 function sr=replaceatom(s,ind,atom)

#  usage:  sr=replaceatom(s,ind,atom) 
#
#       replaces an atom ind1 with other atom
#
#       s     input structure
#       ind   index of replaced atom
#       atom  new atom name 
#       sr    output structure 
#
# example: sr=replaceatom (s,3,"Ge")   
#

     sr=s;
     pos(:)=sr.pos(ind,:); 
     sr=rmatom(sr,ind);
     sr=addatom(sr,atom,pos,ind);

  end 
