function dis=a2adist(s,ind1,ind2)

#   usage:   dis=a2adist(s,ind1,ind2)
#
#         calculates  distance between atoms ind1 and ind2,
#         periodic images are not taken into account
#
#         dis       result
#         s         input structure
#         ind1      index of first atom
#         ind2      index of second atom
#
#   example:
#
#       dis=a2adist(s,1,2)
#


        dpos(:)=s.pos(ind2,:)-s.pos(ind1,:);
        dposc(:)=s.lat2car'*dpos(:);
        dis=sqrt(dposc*dposc');

end
