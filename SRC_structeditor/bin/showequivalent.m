  function lea=showequivalent(s,ind,spgn)

#   usage:  lea=showequivalent(s,ind,spgn) 
#
#         outputs list of equivalent atoms with atom ind
#
#         s      input structure  
#         ind    index of an atom
#         spgn   space group name, if not specified sym. operations 
#                are denerated
#         lea    output vector containing indexes of eqivalent atoms
#     
#   example:  el=showequivalent (s2,5,'P63mc')
#

   if (nargin == 3)

      s1.a=s.a;
      s1.alpha=s.alpha; 	
      s1.pos(:)=s.pos(ind,:);
      s1.spgn=spgn;
      s1.aname(:)=s.aname(ind,:);

      save struct.octave s1;

      system("structgen struct.octave");
      s2=loadstruct("struct.octave.struct");

      ilea=0;
      for j=1:s.nat
          for i=1:s2.nat
              dpos(:)=s2.pos(i,:)-s.pos(j,:);
              dposc(:)=s.lat2car'*dpos(:);
              dis=sqrt(dposc*dposc');
              if (dis < 0.0001)
                 ilea=ilea+1 ;               
                 lea(ilea)=j ;
              endif
          end
      end	
           
      system("rm struct.octave");
      system("rm struct.octave.struct");

   else

      printf("showequivalent:  no spgroup specified, sym. op. are generated\n")
      pos(1:3)=s.pos(ind,1:3);
      aname(:)=s.aname(ind,:);
      savestruct(s,"tmp.struct",1);
      system("readwrite  tmp.struct read_sym  sym ");
      load struct.read.sym
      system("rm struct.read.sym");

      nnew=0;
      for i=1:sym.nsym
          osym(1:3,1:3)=sym.rot(1:3,1:3,i);
          posadd(1:3)=osym(1:3,1:3)*pos(1:3)';
          posadd(1:3)=posadd(1:3)+sym.tr(1:3,i)';

          for j=1:3
             if (posadd(1) >= 1)
                posadd(1)=posadd(1)-1;
             endif
             if (posadd(2) >= 1)
                posadd(2)=posadd(2)-1;
             endif
             if (posadd(3) >= 1)
                posadd(3)=posadd(3)-1;
             endif
             if (posadd(1) < 0)
                posadd(1)=posadd(1)+1;
             endif
             if (posadd(2) < 0)
                posadd(2)=posadd(2)+1;
             endif
             if (posadd(3) < 0)
                posadd(3)=posadd(3)+1;
             endif
          end
          new=1;
          for  j=1:nnew
            dpos(:)=posnew(j,:)-posadd(:)';
             dis=sqrt(dpos*dpos');
             if (dis < 0.0001)
                new=0;
             endif
          end
          if (new)
             nnew=nnew+1;
             posnew(nnew,:)=posadd(1:3) ;
          endif
       end

       ilea=0;
       for j=1:s.nat
           for i=1:nnew
              dpos(:)=posnew(i,:)-s.pos(j,:);
              dposc(:)=s.lat2car'*dpos(:);
              dis=sqrt(dposc*dposc');
              if (dis < 0.0001)
                 ilea=ilea+1 ;
                 lea(ilea)=j ;
              endif
          end
      end

    endif

  end
