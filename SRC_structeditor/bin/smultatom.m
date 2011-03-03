 function sr=smultatom(s,ind,spgn)

#  usage:  sr=smultatom(s,ind,spgn)
#
#       creates symmetry equivalent positions for a given atom
#
#       s      input structure 
#       ind    index of a multiplied atom 
#       spgn   space group symbol, if not specified sym. operations
#              are generated 
#       sr     output structure
#
#  example:  s2=smultatom (s1,5,"P63mc")
#
     
   if (nargin == 3)

      sr=s; 

      s1.a=s.a;
      s1.alpha=s.alpha; 	
      s1.pos(:)=s.pos(ind,:);
      s1.spgn=spgn;
      s1.aname(:)=s.aname(ind,:);

      save struct.octave s1;

      system("structgen struct.octave");
      s2=loadstruct("struct.octave.struct");
      sr=rmatom(sr,ind);
      ind1=ind;
      for i=1:s2.nat
          pos(:)=s2.pos(i,:);
          aname(:)=s.aname(ind,:);
          sr=addatom(sr,aname,pos,ind1);
          ind1=ind1+1;
      end
	            
      system("rm struct.octave");
      system("rm struct.octave.struct");

   else

     sr=s;
     printf("smultatom:  no spgroup specified, sym. op. are generated\n")
     pos(1:3)=sr.pos(ind,1:3);
     aname(:)=sr.aname(ind,:);
     sr=rmatom(sr,ind);
     savestruct(sr,"tmp.struct",1);
     system("readwrite  tmp.struct read_sym  sym ");
     load struct.read.sym
     system("rm struct.read.sym");

     nnew=0;

     for i=1:sym.nsym
         osym(1:3,1:3)=sym.rot(1:3,1:3,i); 
         posadd(1:3)=osym(1:3,1:3)*pos';
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
            sr=addatom(sr,aname,posadd,ind); 
         endif 
      end
 
   endif

 end
