 function sr=makeconventional(s)

#    usage:  sr=makeconventional(s)
#
#         convestrs structure into the conventional form (removes centering)
#            
#         sr    output structure
#         s     input structure  
#
#    example : sr=makeconventional(s)
#

         sr=s;
         if (s.lattic == 'P')
            ntr=1;
            tr(1,1:3)=[0.0,0.0,0.0];
            sr.lattic="P";
         endif 
         if (s.lattic == 'H')
            ntr=1;
            tr(1,1:3)=[0.0,0.0,0.0];
            sr.lattic="H";
         endif 
         if (s.lattic == 'F') 
            ntr=4;
            tr(1,1:3)=[0.0,0.0,0.0];
            tr(2,1:3)=[0.5,0.5,0.0];           
            tr(3,1:3)=[0.5,0.0,0.5];           
            tr(4,1:3)=[0.0,0.5,0.5];           
            sr.lattic="P";
         endif
         if (s.lattic == 'B')
            ntr=2;
            tr(1,1:3)=[0.0,0.0,0.0];
            tr(2,1:3)=[0.5,0.5,0.5];            
            sr.lattic="P";
         endif
         if (s.lattic == 'R') 
            ntr=3;
            tr(1,1:3)=[0.0,0.0d0,0.0d0];
            tr(2,1:3)=[1.0/3.0,2.0/3.0,2.0/3.0];            
            tr(3,1:3)=[2.0/3.0,1.0/3.0,1.0/3.0];

            hex2ort=[ 0 1 0 ; sqrt(3)/2 -0.5 0; 0 0 1 ]; 
            ort2rho=[ 1/sqrt(3) 1/sqrt(3) -2/sqrt(3) ; -1 1 0 ; 1 1 1 ];
            hex2rho=hex2ort*ort2rho;
            rho2hex=inv(hex2rho);
            for i=1:nat
                sr.pos(i,:)=s.pos(i,:)*rho2hex;
            end
            sr.lattic="H";
         endif
         if (s.lattic == 'CXY') 
            ntr=2;
            tr(1,1:3)=[0.0,0.0,0.0];
            tr(2,1:3)=[0.5,0.5,0.0];            
            sr.lattic="P";
         endif
         if (s.lattic == 'CXZ')
            ntr=2;
            tr(1,1:3)=[0.0,0.0,0.0];
            tr(2,1:3)=[0.5,0.0,0.5];            
            sr.lattic="P";
         endif
         if (s.lattic == 'CYZ') 
            ntr=2;
            tr(1,1:3)=[0.0,0.0,0.0];
            tr(2,1:3)=[0.0,0.5,0.5];            
            sr.lattic="P";
         endif

         ii=s.nat;
         for i=1:s.nat
             for k=2:ntr          
               ii=ii+1;
               pos(1:3)=s.pos(i,1:3)+tr(k,1:3);
               for j=1:3      
                   if (pos(j) < 0) 
                      pos(j)=pos(j)+1;
                   endif
                   if (pos(j) > 1) 
                      pos(j)=pos(j)-1;
                   endif
               end
               sr=copyatom(sr,i,ii,pos);

            end
         end

         sr.brlat=sr.lat2car;

 end
