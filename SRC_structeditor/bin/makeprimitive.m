 function sr=makeprimitive(s)

#  usage:  sr=makeprimitive(s)
#
#         converts structure to the primitive form (fcc to rhomboherral)
#
#         sr    output structure
#         s     input structure  
#
#    example : sr=makeprimitive(s)
#

         sr=s; 
        
         sr.a(1)=sqrt(sr.brlat(1,1:3)*sr.brlat(1,1:3)');
         sr.a(2)=sqrt(sr.brlat(2,1:3)*sr.brlat(2,1:3)');
         sr.a(3)=sqrt(sr.brlat(3,1:3)*sr.brlat(3,1:3)');

         sr.alpha(1)=180.0d0*acos((sr.brlat(2,1:3)*sr.brlat(3,1:3)')/(sr.a(2)*sr.a(3)))/pi;
         sr.alpha(2)=180.0d0*acos((sr.brlat(1,1:3)*sr.brlat(3,1:3)')/(sr.a(1)*sr.a(3)))/pi;
         sr.alpha(3)=180.0d0*acos((sr.brlat(2,1:3)*sr.brlat(1,1:3)')/(sr.a(2)*sr.a(1)))/pi;
         tra(:,:)=inv(sr.brlat);
         tra=sr.lat2car*tra;
         for i=1:sr.nat
             sr.pos(i,1:3)=sr.pos(i,1:3)*tra;
         end
         sr.lat2car=sr.brlat;
         sr.lattic="P";

 end
