function elmzz=getazz(elem)

#      usage:  elmzz=getazz(elem)  
#
#             converts atomic name into atomic number
#
#      example:
#
#            elmzz=getazz("Ga") 
#
#

     elmzz=1;
     if (elem(1:2) == "Ac")
         elmzz= 89 ;
     endif
     if (elem(1:2) == "Ag")
         elmzz= 47 ;
     endif
     if (elem(1:2) == "Al")
         elmzz= 13 ;
     endif
     if (elem(1:2) == "Am")
         elmzz= 95 ;
     endif
     if (elem(1:2) == "Ar")
         elmzz= 18 ;
     endif
     if (elem(1:2) == "As")
         elmzz= 33 ;
     endif
     if (elem(1:2) == "At")
         elmzz= 85 ;
     endif
     if (elem(1:2) == "Au")
         elmzz= 79 ;
     endif
     if (elem(1:1) == "B")
          elmzz=  5 ;
     endif
     if (elem(1:2) == "Ba")
         elmzz= 56 ;
     endif
     if (elem(1:2) == "Be")
         elmzz=  4 ;
     endif
     if (elem(1:2) == "Bi")
         elmzz= 83 ;
     endif
     if (elem(1:2) == "Bk")
         elmzz= 97 ;
     endif
     if (elem(1:2) == "Br")
         elmzz= 35 ;
     endif
     if (elem(1:1) == "C")
          elmzz=  6 ;
     endif
     if (elem(1:2) == "Ca")
         elmzz= 20 ;
     endif
     if (elem(1:2) == "Cd")
         elmzz= 48 ;
     endif
     if (elem(1:2) == "Ce")
         elmzz= 58 ;
     endif
     if (elem(1:2) == "Cf")
         elmzz= 98 ;
     endif
     if (elem(1:2) == "Cl")
         elmzz= 17 ;
     endif
     if (elem(1:2) == "Cm")
         elmzz= 96 ;
     endif
     if (elem(1:2) == "Co")
         elmzz= 27 ;
     endif
     if (elem(1:2) == "Cr")
         elmzz= 24 ;
     endif
     if (elem(1:2) == "Cs")
         elmzz= 55 ;
     endif
     if (elem(1:2) == "Cu")
         elmzz= 29 ;
     endif
     if (elem(1:2) == "Dy")
         elmzz= 66 ;
     endif
     if (elem(1:2) == "Er")
         elmzz= 68 ;
     endif
     if (elem(1:2) == "Es")
         elmzz= 99 ;
     endif
     if (elem(1:2) == "Eu")
         elmzz= 63 ;
     endif
     if (elem(1:1)=="F")
          elmzz=  9 ;
     endif
     if (elem(1:2)=="Fe")
         elmzz= 26 ;
     endif
     if (elem(1:2)=="Fm")
         elmzz=100 ;
     endif
     if (elem(1:2)=="Fr")
         elmzz= 87 ;
     endif
     if (elem(1:2)=="Ga")
         elmzz= 31 ;
     endif
     if (elem(1:2)=="Gd")
         elmzz= 64 ;
     endif
     if (elem(1:2)=="Ge")
         elmzz= 32 ;
     endif
     if (elem(1:1)=="H")
          elmzz=  1 ;
     endif
     if (elem(1:2)=="He")
         elmzz=  2 ;
     endif
     if (elem(1:2)=="Hf")
         elmzz= 72 ;
     endif
     if (elem(1:2)=="Hg")
         elmzz= 80 ;
     endif
     if (elem(1:2)=="Ho")
         elmzz= 67 ;
     endif
     if (elem(1:1)=="I")
          elmzz= 53 ;
     endif
     if (elem(1:2)=="In")
         elmzz= 49 ;
     endif
     if (elem(1:2)=="Ir")
         elmzz= 77 ;
     endif
     if (elem(1:1)=="K")
          elmzz= 19 ;
     endif
     if (elem(1:2)=="Kr")
         elmzz= 36 ;
     endif
     if (elem(1:2)=="La")
         elmzz= 57 ;
     endif
     if (elem(1:2)=="Li")
         elmzz=  3 ;
     endif
     if (elem(1:2)=="Lr")
         elmzz=103 ;
     endif
     if (elem(1:2)=="Lu")
         elmzz= 71 ;
     endif
     if (elem(1:2)=="Md")
         elmzz=101 ;
     endif
     if (elem(1:2)=="Mg")
         elmzz= 12 ;
     endif
     if (elem(1:2)=="Mn")
         elmzz= 25 ;
     endif
     if (elem(1:2)=="Mo")
         elmzz= 42 ;
     endif
     if (elem(1:1)=="N")
          elmzz=  7 ;
     endif
     if (elem(1:2)=="Na")
         elmzz= 11 ;
     endif
     if (elem(1:2)=="Nb")
         elmzz= 41 ;
     endif
     if (elem(1:2)=="Nd")
         elmzz= 60 ;
     endif
     if (elem(1:2)=="Ne")
         elmzz= 10 ;
     endif
     if (elem(1:2)=="Ni")
         elmzz= 28 ;
     endif
     if (elem(1:2)=="No")
         elmzz=102 ;
     endif
     if (elem(1:2)=="Np")
         elmzz= 93 ;
     endif
     if (elem(1:1)=="O")
          elmzz=  8 ;
     endif
     if (elem(1:2)=="Os")
         elmzz= 76 ;
     endif
     if (elem(1:1)=="P")
          elmzz= 15 ;
     endif
     if (elem(1:2)=="Pa")
         elmzz= 91 ;
     endif
     if (elem(1:2)=="Pb")
         elmzz= 82 ;
     endif
     if (elem(1:2)=="Pd")
         elmzz= 46 ;
     endif
     if (elem(1:2)=="Pm")
         elmzz= 61 ;
     endif
     if (elem(1:2)=="Po")
         elmzz= 84 ;
     endif
     if (elem(1:2)=="Pr")
         elmzz= 59 ;
     endif
     if (elem(1:2)=="Pt")
         elmzz= 78 ;
     endif
     if (elem(1:2)=="Pu")
         elmzz= 94 ;
     endif
     if (elem(1:2)=="Ra")
         elmzz= 88 ;
     endif
     if (elem(1:2)=="Rb")
         elmzz= 37 ;
     endif
     if (elem(1:2)=="Re")
         elmzz= 75 ;
     endif
     if (elem(1:2)=="Rh")
         elmzz= 45 ;
     endif
     if (elem(1:2)=="Rn")
         elmzz= 86 ;
     endif
     if (elem(1:2)=="Ru")
         elmzz= 44 ;
     endif
     if (elem(1:1)=="S")
          elmzz= 16 ;
     endif
     if (elem(1:2)=="Sb")
         elmzz= 51 ;
     endif
     if (elem(1:2)=="Sc")
         elmzz= 21 ;
     endif
     if (elem(1:2)=="Se")
         elmzz= 34 ;
     endif
     if (elem(1:2)=="Si")
         elmzz= 14 ;
     endif
     if (elem(1:2)=="Sm")
         elmzz= 62 ;
     endif
     if (elem(1:2)=="Sn")
         elmzz= 50 ;
     endif
     if (elem(1:2)=="Sr")
         elmzz= 38 ;
     endif
     if (elem(1:2)=="Ta")
         elmzz= 73 ;
     endif
     if (elem(1:2)=="Tb")
         elmzz= 65 ;
     endif
     if (elem(1:2)=="Tc")
         elmzz= 43 ;
     endif
     if (elem(1:2)=="Te")
         elmzz= 52 ;
     endif
     if (elem(1:2)=="Th")
         elmzz= 90 ;
     endif
     if (elem(1:2)=="Ti")
         elmzz= 22 ;
     endif
     if (elem(1:2)=="Tl")
         elmzz= 81 ;
     endif
     if (elem(1:2)=="Tm")
         elmzz= 69 ;
     endif
     if (elem(1:1)=="U")
          elmzz= 92 ;
     endif
     if (elem(1:1)=="V")
          elmzz= 23 ;
     endif
     if (elem(1:1)=="W")
          elmzz= 74 ;
     endif
     if (elem(1:2)=="Xe")
         elmzz= 54 ;
     endif
     if (elem(1:1)=="Y")
          elmzz= 39 ;
     endif
     if (elem(1:2)=="Yb")
         elmzz= 70 ;
     endif
     if (elem(1:2)=="Zn")
         elmzz= 30 ;
     endif
     if (elem(1:2)=="Zr")
         elmzz= 40 ;
     endif

end
