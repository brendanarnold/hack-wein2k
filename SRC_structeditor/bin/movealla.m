function sr=movealla(s,vec)

#  usage:  sr=movealla(s,vec)
#
#        moves all atoms with vector vec
#
#       s       input structure
#       vec     vector
#       sr      output structure  
#
#  example:
#
#       sr=movealla(s,[0 0 0.5])
#

      sr=s;
      svec=[0.00001 0.00001 0.00001];
      for i=1:s.nat
          sr.pos(i,1:3)=mod((sr.pos(i,1:3)+vec+svec),1)-svec;
      end

end
