function showstruct(s)

#    usage:  showstruct(s)
#
#            displays structure stored in a variable s  
#     

      save  struct.write s;
      system(["readwrite "," s.struct "," write "," s " ]);
      system("struct2dx -nomm -new -f s >& error.struct2dx", 1 , "async");
      system("rm -f s.outputstr2mol struct2mol.def struct2mol.error");
              
end
