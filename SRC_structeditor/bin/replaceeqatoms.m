 function sr=replaceeqatoms(s,ind,atom,spgn)

#  usage:  sr=replaceeqatoms(s,ind,atom,spgn)
#
#          replaces atom ind and all eqivalent with other atoms
#          
#        s        input structure  
#        ind      index of removed atoms
#        atom     name of new atoms
#        spgn     space group name, if not present sym. operations
#                 are genereted
#        sr       output structure
#
#  example: sr=replaceeqatoms (s,3,"Ge","P63mc") 
#   


     if (nargin == 4)

        lea=showequivalent(s,ind,spgn);

     else

        lea=showequivalent(s,ind);

     endif

        nle=columns(lea);

        sr=s;
        for i=nle:-1:1
            pos(i,:)=s.pos(lea(i),:);
            sr=rmatom(sr,lea(i));
        end
        for i=1:nle
            sr=addatom(sr,atom,pos(i,:),ind);
        end


 end
