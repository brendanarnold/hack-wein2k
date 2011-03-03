 function sr=sshift(s,ind,vec,spgn)

#   usage:  sr=sshift(s,ind,vec,spgn)
#
#        shifts atom ind with vector vec and all equivalent with
#        corresponding equivalents vectors
#
#        s      input structure
#        ind    index of an atom  
#        vec    displacement vector
#        spgn   space group name (optional), if not present
#               sym. op. are generated
#
#   example: s2=sshift (s1,5,[ 0.1 0.1 0.1])  
#

     pos(:)=s.pos(ind,:);
     pos(:)=pos(:)+vec(:);
     atom(:)=s.aname(ind,:);

     if (nargin == 4)
        s1=addeqatom(s,atom,pos,0,spgn);
        lea=showequivalent(s,ind,spgn);
     else
        s1=addeqatom(s,atom,pos,0);
        lea=showequivalent(s,ind);
     endif

     nle=columns(lea);

     sr=s;
    
     dpos(:)=vec*s.lat2car;
     dis0=sqrt(dpos*dpos');

     if (s1.nat != s.nat+nle)
        printf("sshift: displacement is incompatible with space group\n")
     else

        ibb=0; 
        for i=1:nle
            j=lea(i);
            ib=1;
            for k=s.nat+1:s.nat+nle
               
                dis=a2adist(s1,j,k)-dis0;

                if (abs(dis) < 0.0001)                    
                   sr.pos(j,:)=s1.pos(k,:);
                   ib=0; 
                   break;
                endif 

            end

            if (ib)
               ibb=1;
               break;
            endif 

        end

        if (ibb)
           k=s.nat;
           for i=1:nle
               j=lea(i);
               k=k+1 ;
               sr.pos(j,:)=s1.pos(k,:);
           end
           printf("sshift: equivalent atoms may be interchanged\n");
        endif
 
     endif

     
 end

