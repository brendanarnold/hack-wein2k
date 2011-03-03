function r0=getar0(zz)

#     usage:  r0=getar0(zz)  
#    
#          calculates r0 from atomic number zz
#
#     example:
#
#           r0=getar0(5)
#
#

     r0=5.0e-6 ;
     if (zz<72) 
        r0=1.0e-5;
     endif	
     if (zz<37) 
        r0=5.0e-5;
     endif
     if (zz<19) 
        r0=1.0e-4;
     endif
end
