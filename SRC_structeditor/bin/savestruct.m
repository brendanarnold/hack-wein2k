  function savestruct(s,filename,dosym)

#   usage:  savestruct(s,filename,dosym)
#          
#           saves crystal structure stored in variable s as
#           Wien2k structfile.
# 
#           s              is a structure
#           filename       is a text variable
#           dosym          1 - with symmetry, 0 - no symmetry
#
#   Example:
#
#         savestruct(s,"GaN.struct",0)
#
             
      if (nargin == 2)
         dosym=1;
      endif

      format long
      save  struct.write s;        
      system(["readwrite ",filename," write "," s ",  int2str(dosym)]);
      system("rm struct.write");
     
      if (dosym)
         system(["cp ", filename," ", "tmpsave.struct"]);
         system("xncm ncmsymmetry -f tmpsave -nomm>& error.ncmsym ");
         system(["cp ", "tmpsave.struct_ncmsym ", filename]);
         system("rm -f tmpsave.* ncmsymmetry.error ncmsymmetry.def ");
      endif	 

  end	
