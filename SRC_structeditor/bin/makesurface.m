 function sr=makesurface(s,n,ind,depth,vac)

# usage:  sr=makesurface(s,n,ind,depth,vac)
#
#      creates surface for a given unitcells
#
#      s       input structure   
#      n       normal vectror (in lattice coordinates)
#      ind     an index of an atom which should be in (0 0 0) 
#      depth   thickness of the material  
#      vac     thickness of the vacum layer
#
# evample: sr=makesurface(s,[0 0 1],1,30.0,20.0)
#
    
   imax=3;
   jest=0;
   for i1=0:imax
       if (jest)
          break;
       endif 
       j1=ji(i1); 
       for i2=0:imax
           if (jest)
              break;
           endif 
           j2=ji(i2); 
           for i3=0:imax
               if (jest)
                  break;
               endif 
               j3=ji(i3);
               if (!(i1==0 & i2==0 & i3==0))
                  np(:)=[j1,j2,j3];
                  if (n*np'==0)
                     np1=np;
                     jest=1;
                  end
               end       
           end
       end
   end
               
   dnp1=sqrt(np1*np1');
    
   jest=0;
   for i1=0:imax
       if (jest)
          break;
       endif 
       j1=ji(i1); 
       for i2=0:imax
           if (jest)
              break;
           endif 
           j2=ji(i2); 
           for i3=0:imax
               if (jest)
                  break;
               endif 
               j3=ji(i3);
               if (!(i1==0 & i2==0 & i3==0))
                  np(:)=[j1,j2,j3];
                  dnp=sqrt(np*np');
                  spp1=abs(np1*np'/(dnp*dnp1));
                  if (n*np'==0 & spp1 < 0.9999)
                     np2=np;
                     jest=1;
                  end
               end       
           end
       end
   end

   a(1,1:3)=np1(1:3); 
   a(2,1:3)=np2(1:3); 
   a(3,1:3)=n(1:3);

    printf("Lattice vectors: \n")
    a

    brlat=a*s.brlat;
    c=sqrt(brlat(3,1:3)*brlat(3,1:3)');
    a(3,1:3)=a(3,1:3)*depth/c;
 
   if (ind == 0 | ind > s.nat)   
      sr=makesupercell(s,a);
   else
      vec(:)=-s.pos(ind,:);
      s=movealla(s,vec);
      sr=makesupercell(s,a);
   endif

   sr.a(3)=sr.a(3)+vac;
   tm=[1,0,0;
       0,1,0;
       0,0,sr.a(3)/(sr.a(3)-vac)];
   sr.lat2car= tm*sr.lat2car;
   sr.brlat= tm*sr.brlat;

   for i=1:sr.nat
       sr.pos(i,3)=sr.pos(i,3)*(sr.a(3)-vac)/sr.a(3);
   end

 end


 function ji=ji(i)

       if (mod(i,2) == 0)
          ji=i/2;
       else
          ji=-(i+1)/2;
       end

 end
