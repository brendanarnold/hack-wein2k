  function sr=rmeqatoms(s,ind,spgn)

#  usage:  sr=rmeqatoms(s,ind,spgn)
#
#       removes atom ind and all equivalent
#
#       s     input structure 
#       ind   index of an atom
#       spgn  space group name
#       sr    output structure 
#
#  example:  sr=rmeqatoms(s,5,'P63mc')
#

   if (nargin == 3)

     lea=showequivalent(s,ind,spgn);

   else

     lea=showequivalent(s,ind);

   end

     nle=columns(lea);
      
     sr=s;
     for i=nle:-1:1
         sr=rmatom(sr,lea(i));
     end

  end
