  function s = loadstruct(filename)

#   usage:  s = loadstruct(filename)
#
#           reads Wien2k structfile (filename) and returns 
#           apopriate structure variable (s)
#
#           filename       is a text variable
#           s              is a structure
#
#   Example:
#
#         s = loadstruct("GaN.struct")
#

      format long

      system(["readwrite ",filename," read "," s " ]);
      load struct.read;        
      system("rm struct.read");

  end	
