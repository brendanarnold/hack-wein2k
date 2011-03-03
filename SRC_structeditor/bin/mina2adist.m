function dis=mina2adist(s,ind1,ind2)

#   usage:   dis=mina2adist(s,ind1,ind2)
#
#         calculates minimum distance between atoms ind1 and ind2 
#         taking into accound periodic images         
# 
#         dis       result
#         s         input structure
#         ind1      index of first atom
#         ind2      index of second atom
#
#   example:    
#
#       dis=mina2adist(s,1,2)  
#


    dis=1e10;
    for i1=-1:1
    for i2=-1:1
    for i3=-1:1

        dpos(:)=s.pos(ind2,:)-s.pos(ind1,:)+[i1,i2,i3];
        dposc(:)=s.lat2car'*dpos(:);
        d=sqrt(dposc*dposc');
        dis=min(dis,d);

    end
    end
    end



end
