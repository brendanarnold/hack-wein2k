 function sr=makesupercell(s,a)

#   usage:  sr=makesupercell(s,a)
#
#      creates supercell based on structure s
#
#       s     input structure 
#       a     new brave lattice in the basis of old lattice vectors
#             (rowwise), works only with non-centered WIEN structures
#       
#   example:  sr=makesupercell(s,[1 0 0; 0 1 0; 0 0 2]) 
#

     if (s.lattic == "H" | s.lattic == "P")

      a1max=max(a(:,1))+1;
      a2max=max(a(:,2))+1;
      a3max=max(a(:,3))+1;
   
      a1min=min(a(:,1));
      a2min=min(a(:,2));
      a3min=min(a(:,3));
   
      if (a1min > 0) 
         a1min=0;
      endif
      if (a2min > 0) 
         a2min=0;
      endif
      if (a3min > 0) 
         a3min=0;
      endif
   
      if (a1max < 0) 
         a1max=0;
      endif
      if (a2max < 0) 
         a2max=0;
      endif
      if (a3max < 0) 
         a3max=0;
      endif
   
      inva=inv(a);
   
      sr.nat=0;

      for i=1:s.nat

          for i1=a1min:a1max
              for i2=a2min:a2max
                  for i3=a3min:a3max

                         poso(:)=s.pos(i,:)+[i1 i2 i3];              
                         posn(:)=poso*inva;

                         if ((posn <= 0.99999) & (posn >= -0.0001) )
                            jest=0;  
                            for j=1:sr.nat

                                vec=sr.pos(j,:)-posn(:)';
                                dd=sqrt(vec*vec');
                                
                                if (dd < 0.00001) 
                                   jest=1;
                                   break
                                endif
                            end
                            if (!jest)
                               sr=copyatom2(s,sr,i,0,posn);
                            endif
                         endif
   
                  end
              end
          end

      end

      sr.brlat=a*s.brlat;
      sr.lat2car=sr.brlat;

      sr.a(1)=sqrt(sr.brlat(1,1:3)*sr.brlat(1,1:3)');
      sr.a(2)=sqrt(sr.brlat(2,1:3)*sr.brlat(2,1:3)');
      sr.a(3)=sqrt(sr.brlat(3,1:3)*sr.brlat(3,1:3)');

      sr.alpha(1)=180.0d0*acos((sr.brlat(2,1:3)*sr.brlat(3,1:3)')/(sr.a(2)*sr.a(3)))/pi;
      sr.alpha(2)=180.0d0*acos((sr.brlat(1,1:3)*sr.brlat(3,1:3)')/(sr.a(1)*sr.a(3)))/pi;
      sr.alpha(3)=180.0d0*acos((sr.brlat(2,1:3)*sr.brlat(1,1:3)')/(sr.a(2)*sr.a(1)))/pi;

      sr.lattic='P';

     else
        printf("makesupercell:  remove contering with mace primitive or makeconventional\n")
        sr=s;        
     endif
     
 end

         
