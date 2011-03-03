function sr=addeqatom(s,atom,pos,ind,spgn)

#  usage:  sr=addeqatom(s,atom,pos,ind)
#
#        adds an atom and all equivalent to the structure s and returns all as structure sr
#
#        sr      output structure
#        s       matlab structure used to store structural information
#        atom    identyfies atom, can be string variable containing
#                atomic symbol or scalar with atomic number
#        pos     position in unit cell
#        ind     a index of new atom (position in the list),
#                if 0 then atom is added to the end of the list
#        spgn    spacegroup name, if not specified sym. operations
#                are generated         
# 
#  example:
#
#         sr=addeqatom(s,"N",[0,0,0.1],3,"P63mc")
#

     sr=addatom(s,atom,pos,ind);

     if (ind ==0 )
        ind=sr.nat;
     end

     if (nargin == 5)
 
	sr=smultatom(sr,ind,spgn);

     else

        sr=smultatom(sr,ind);

     endif 

 end
