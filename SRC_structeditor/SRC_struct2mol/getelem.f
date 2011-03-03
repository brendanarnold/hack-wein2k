 subroutine getelem(nspec,specname,specrad,specrgb,speccr,nbo,indb,bonddis)
   
   use struct, only: nat,aname,a2b

   integer     nspec
   character*6 specname(*)
   real*8      specrad(*)
   real*8      speccr(*)
   real*8      specrgb(3,*)
   integer     nbo
   integer     indb(2,*) 
   real*8      bonddis(*)  

   character elem*2,elemname*20
   integer elemzz
   real*8  elemw,elemcr,elemvr,elemrgb(3)

   character*10 tmpname(0:1000)
   real*8       tmpvr(0:1000)
   real*8       tmpcr(0:1000)
   real*8       tmprgb(0:1000,1:3)
   integer i,tmp(0:1000),melem

   tmp(1:1000)=0
     melem=0

   do i=1,nat
      
      elem='  '
      elem=trim(aname(i)(1:2))
      if (elem(2:2).eq.'0') elem=elem(1:1)//' '
      if (elem(2:2).eq.'1') elem=elem(1:1)//' '
      if (elem(2:2).eq.'2') elem=elem(1:1)//' '
      if (elem(2:2).eq.'3') elem=elem(1:1)//' '
      if (elem(2:2).eq.'4') elem=elem(1:1)//' '
      if (elem(2:2).eq.'5') elem=elem(1:1)//' '
      if (elem(2:2).eq.'6') elem=elem(1:1)//' '
      if (elem(2:2).eq.'7') elem=elem(1:1)//' '
      if (elem(2:2).eq.'8') elem=elem(1:1)//' '
      if (elem(2:2).eq.'9') elem=elem(1:1)//' '
      write(*,*) elem
     if (elem.eq.'X') then
        elemname='Dummy'
        elemzz=   0
        elemw=  0.000000
        elemcr=  0.100000
        elemvr=  1.000000
        elemrgb(1:3)=(/  0.000000,  0.000000,  0.984970/)
     else if (elem.eq.'H') then
        elemname='Hydrogen'
        elemzz=   1
        elemw=  1.007900
        elemcr=  0.330000
        elemvr=  1.170000
        elemrgb(1:3)=(/  0.984970,  0.984970,  0.984970/)
     else if (elem.eq.'He') then
        elemname='Helium'
        elemzz=   2
        elemw=  4.002600
        elemcr=  0.000000
        elemvr=  0.400000
        elemrgb(1:3)=(/  0.984970,  0.984970,  0.984970/)
     else if (elem.eq.'Li') then
        elemname='Lithium'
        elemzz=   3
        elemw=  6.941000
        elemcr=  0.680000
        elemvr=  1.010000
        elemrgb(1:3)=(/  0.450891,  0.790561,  0.910799/)
     else if (elem.eq.'Be') then
        elemname='Beryllium'
        elemzz=   4
        elemw=  9.012200
        elemcr=  0.350000
        elemvr=  0.520000
        elemrgb(1:3)=(/  0.075148,  0.536560,  0.760502/)
     else if (elem.eq.'B') then
        elemname='Boron'
        elemzz=   5
        elemw= 10.811000
        elemcr=  0.820000
        elemvr=  1.700000
        elemrgb(1:3)=(/  0.984970,  0.611708,  0.966409/)
     else if (elem.eq.'C') then
        elemname='Carbon'
        elemzz=   6
        elemw= 12.011000
        elemcr=  0.770000
        elemvr=  1.750000
        elemrgb(1:3)=(/  0.533554,  0.910348,  0.595176/)
     else if (elem.eq.'N') then
        elemname='Nitrogen'
        elemzz=   7
        elemw= 14.006700
        elemcr=  0.750000
        elemvr=  1.550000
        elemrgb(1:3)=(/  0.804088,  0.455399,  0.984970/)
     else if (elem.eq.'O') then
        elemname='Oxygen'
        elemzz=   8
        elemw= 15.999400
        elemcr=  0.730000
        elemvr=  1.400000
        elemrgb(1:3)=(/  0.984970,  0.300594,  0.240475/)
     else if (elem.eq.'F') then
        elemname='Fluorine'
        elemzz=   9
        elemw= 18.998400
        elemcr=  0.720000
        elemvr=  1.300000
        elemrgb(1:3)=(/  0.339671,  0.925829,  0.804088/)
     else if (elem.eq.'Ne') then
        elemname='Neon'
        elemzz=  10
        elemw= 20.179701
        elemcr=  0.000000
        elemvr=  0.700000
        elemrgb(1:3)=(/  0.984970,  0.984970,  0.984970/)
     else if (elem.eq.'Na') then
        elemname='Sodium'
        elemzz=  11
        elemw= 22.989799
        elemcr=  0.970000
        elemvr=  1.450000
        elemrgb(1:3)=(/  0.450891,  0.790561,  0.910799/)
     else if (elem.eq.'Mg') then
        elemname='Magnesium'
        elemzz=  12
        elemw= 24.305000
        elemcr=  1.100000
        elemvr=  1.640000
        elemrgb(1:3)=(/  0.075148,  0.536560,  0.760502/)
     else if (elem.eq.'Al') then
        elemname='Aluminum'
        elemzz=  13
        elemw= 26.981501
        elemcr=  1.180000
        elemvr=  2.010000
        elemrgb(1:3)=(/  0.339671,  0.925829,  0.804088/)
     else if (elem.eq.'Si') then
        elemname='Silicon'
        elemzz=  14
        elemw= 28.085501
        elemcr=  1.110000
        elemvr=  2.000000
        elemrgb(1:3)=(/  0.984970,  0.611708,  0.966409/)
     else if (elem.eq.'P') then
        elemname='Phosphorus'
        elemzz=  15
        elemw= 30.973801
        elemcr=  1.060000
        elemvr=  1.900000
        elemrgb(1:3)=(/  0.984970,  0.855940,  0.736454/)
     else if (elem.eq.'S') then
        elemname='Sulfur'
        elemzz=  16
        elemw= 32.066002
        elemcr=  1.020000
        elemvr=  1.800000
        elemrgb(1:3)=(/  0.981438,  0.984970,  0.573382/)
     else if (elem.eq.'Cl') then
        elemname='Chlorine'
        elemzz=  17
        elemw= 35.452702
        elemcr=  0.990000
        elemvr=  1.770000
        elemrgb(1:3)=(/  0.339671,  0.925829,  0.804088/)
     else if (elem.eq.'Ar') then
        elemname='Argon'
        elemzz=  18
        elemw= 39.948002
        elemcr=  0.000000
        elemvr=  0.970000
        elemrgb(1:3)=(/  0.984970,  0.984970,  0.984970/)
     else if (elem.eq.'K') then
        elemname='Potassium'
        elemzz=  19
        elemw= 39.098301
        elemcr=  1.330000
        elemvr=  1.990000
        elemrgb(1:3)=(/  0.450891,  0.790561,  0.910799/)
     else if (elem.eq.'Ca') then
        elemname='Calcium'
        elemzz=  20
        elemw= 40.077999
        elemcr=  0.990000
        elemvr=  1.480000
        elemrgb(1:3)=(/  0.075148,  0.536560,  0.760502/)
     else if (elem.eq.'Sc') then
        elemname='Scandium'
        elemzz=  21
        elemw= 44.955898
        elemcr=  1.440000
        elemvr=  2.150000
        elemrgb(1:3)=(/  0.635756,  0.916811,  0.730443/)
     else if (elem.eq.'Ti') then
        elemname='Titanium'
        elemzz=  22
        elemw= 47.880001
        elemcr=  1.470000
        elemvr=  2.190000
        elemrgb(1:3)=(/  0.635756,  0.916811,  0.730443/)
     else if (elem.eq.'V') then
        elemname='Vanadium'
        elemzz=  23
        elemw= 50.941502
        elemcr=  1.330000
        elemvr=  1.990000
        elemrgb(1:3)=(/  0.635756,  0.916811,  0.730443/)
     else if (elem.eq.'Cr') then
        elemname='Chromium'
        elemzz=  24
        elemw= 51.996101
        elemcr=  1.350000
        elemvr=  2.010000
        elemrgb(1:3)=(/  0.635756,  0.916811,  0.730443/)
     else if (elem.eq.'Mn') then
        elemname='Manganese'
        elemzz=  25
        elemw= 54.938000
        elemcr=  1.350000
        elemvr=  2.010000
        elemrgb(1:3)=(/  0.635756,  0.916811,  0.730443/)
     else if (elem.eq.'Fe') then
        elemname='Iron'
        elemzz=  26
        elemw= 55.847000
        elemcr=  1.340000
        elemvr=  2.000000
        elemrgb(1:3)=(/  0.635756,  0.916811,  0.730443/)
     else if (elem.eq.'Co') then
        elemname='Cobalt'
        elemzz=  27
        elemw= 58.933201
        elemcr=  1.330000
        elemvr=  1.990000
        elemrgb(1:3)=(/  0.635756,  0.916811,  0.730443/)
     else if (elem.eq.'Ni') then
        elemname='Nickel'
        elemzz=  28
        elemw= 58.693401
        elemcr=  1.500000
        elemvr=  1.810000
        elemrgb(1:3)=(/  0.635756,  0.916811,  0.730443/)
     else if (elem.eq.'Cu') then
        elemname='Copper'
        elemzz=  29
        elemw= 63.546001
        elemcr=  1.520000
        elemvr=  1.540000
        elemrgb(1:3)=(/  0.635756,  0.916811,  0.730443/)
     else if (elem.eq.'Zn') then
        elemname='Zinc'
        elemzz=  30
        elemw= 65.389999
        elemcr=  1.450000
        elemvr=  2.160000
        elemrgb(1:3)=(/  0.635756,  0.916811,  0.730443/)
     else if (elem.eq.'Ga') then
        elemname='Gallium'
        elemzz=  31
        elemw= 69.723000
        elemcr=  1.220000
        elemvr=  1.820000
        elemrgb(1:3)=(/  0.339671,  0.925829,  0.804088/)
     else if (elem.eq.'Ge') then
        elemname='Germanium'
        elemzz=  32
        elemw= 72.610001
        elemcr=  1.170000
        elemvr=  1.750000
        elemrgb(1:3)=(/  0.339671,  0.925829,  0.804088/)
     else if (elem.eq.'As') then
        elemname='Arsenic'
        elemzz=  33
        elemw= 74.921600
        elemcr=  1.210000
        elemvr=  2.200000
        elemrgb(1:3)=(/  0.984970,  0.611708,  0.966409/)
     else if (elem.eq.'Se') then
        elemname='Selenium'
        elemzz=  34
        elemw= 78.959999
        elemcr=  1.220000
        elemvr=  2.000000
        elemrgb(1:3)=(/  0.984970,  0.794469,  0.899977/)
     else if (elem.eq.'Br') then
        elemname='Bromine'
        elemzz=  35
        elemw= 79.903999
        elemcr=  1.210000
        elemvr=  1.950000
        elemrgb(1:3)=(/  0.880739,  0.670324,  0.483956/)
     else if (elem.eq.'Kr') then
        elemname='Krypton'
        elemzz=  36
        elemw= 83.800003
        elemcr=  1.890000
        elemvr=  2.820000
        elemrgb(1:3)=(/  0.984970,  0.984970,  0.984970/)
     else if (elem.eq.'Rb') then
        elemname='Rubidium'
        elemzz=  37
        elemw= 85.467796
        elemcr=  1.470000
        elemvr=  2.190000
        elemrgb(1:3)=(/  0.450891,  0.790561,  0.910799/)
     else if (elem.eq.'Sr') then
        elemname='Strontium'
        elemzz=  38
        elemw= 87.620003
        elemcr=  1.120000
        elemvr=  1.670000
        elemrgb(1:3)=(/  0.075148,  0.536560,  0.760502/)
     else if (elem.eq.'Y') then
        elemname='Yttrium'
        elemzz=  39
        elemw= 88.905899
        elemcr=  1.780000
        elemvr=  2.660000
        elemrgb(1:3)=(/  0.106711,  0.790561,  0.813106/)
     else if (elem.eq.'Zr') then
        elemname='Zirconium'
        elemzz=  40
        elemw= 91.223999
        elemcr=  1.560000
        elemvr=  2.330000
        elemrgb(1:3)=(/  0.106711,  0.790561,  0.813106/)
     else if (elem.eq.'Nb') then
        elemname='Niobium'
        elemzz=  41
        elemw= 92.906403
        elemcr=  1.480000
        elemvr=  2.210000
        elemrgb(1:3)=(/  0.106711,  0.790561,  0.813106/)
     else if (elem.eq.'Mo') then
        elemname='Molybdenum'
        elemzz=  42
        elemw= 95.940002
        elemcr=  1.470000
        elemvr=  2.190000
        elemrgb(1:3)=(/  0.106711,  0.790561,  0.813106/)
     else if (elem.eq.'Tc') then
        elemname='Technetium'
        elemzz=  43
        elemw= 97.907204
        elemcr=  1.350000
        elemvr=  2.010000
        elemrgb(1:3)=(/  0.798076,  0.798076,  0.798076/)
     else if (elem.eq.'Ru') then
        elemname='Ruthenium'
        elemzz=  44
        elemw=101.070000
        elemcr=  1.400000
        elemvr=  2.090000
        elemrgb(1:3)=(/  0.106711,  0.790561,  0.813106/)
     else if (elem.eq.'Rh') then
        elemname='Rhodium'
        elemzz=  45
        elemw=102.905502
        elemcr=  1.450000
        elemvr=  2.160000
        elemrgb(1:3)=(/  0.106711,  0.790561,  0.813106/)
     else if (elem.eq.'Pd') then
        elemname='Palladium'
        elemzz=  46
        elemw=106.419998
        elemcr=  1.500000
        elemvr=  2.240000
        elemrgb(1:3)=(/  0.106711,  0.790561,  0.813106/)
     else if (elem.eq.'Ag') then
        elemname='Silver'
        elemzz=  47
        elemw=107.868202
        elemcr=  1.590000
        elemvr=  2.370000
        elemrgb(1:3)=(/  0.106711,  0.790561,  0.813106/)
     else if (elem.eq.'Cd') then
        elemname='Cadmium'
        elemzz=  48
        elemw=112.411003
        elemcr=  1.690000
        elemvr=  2.520000
        elemrgb(1:3)=(/  0.106711,  0.790561,  0.813106/)
     else if (elem.eq.'In') then
        elemname='Indium'
        elemzz=  49
        elemw=114.818001
        elemcr=  1.630000
        elemvr=  2.430000
        elemrgb(1:3)=(/  0.339671,  0.925829,  0.804088/)
     else if (elem.eq.'Sn') then
        elemname='Tin'
        elemzz=  50
        elemw=118.709999
        elemcr=  1.460000
        elemvr=  2.180000
        elemrgb(1:3)=(/  0.339671,  0.925829,  0.804088/)
     else if (elem.eq.'Sb') then
        elemname='Antimony'
        elemzz=  51
        elemw=121.760002
        elemcr=  1.460000
        elemvr=  2.170000
        elemrgb(1:3)=(/  0.339671,  0.925829,  0.804088/)
     else if (elem.eq.'Te') then
        elemname='Tellurium'
        elemzz=  52
        elemw=127.599998
        elemcr=  1.470000
        elemvr=  2.200000
        elemrgb(1:3)=(/  0.984970,  0.611708,  0.966409/)
     else if (elem.eq.'I') then
        elemname='Iodine'
        elemzz=  53
        elemw=126.904503
        elemcr=  1.330000
        elemvr=  2.100000
        elemrgb(1:3)=(/  0.930337,  0.275043,  0.984970/)
     else if (elem.eq.'Xe') then
        elemname='Xenon'
        elemzz=  54
        elemw=131.289993
        elemcr=  0.000000
        elemvr=  1.100000
        elemrgb(1:3)=(/  0.984970,  0.984970,  0.984970/)
     else if (elem.eq.'Cs') then
        elemname='Cesium'
        elemzz=  55
        elemw=132.905396
        elemcr=  1.670000
        elemvr=  2.490000
        elemrgb(1:3)=(/  0.450891,  0.790561,  0.910799/)
     else if (elem.eq.'Ba') then
        elemname='Barium'
        elemzz=  56
        elemw=137.326996
        elemcr=  1.340000
        elemvr=  2.000000
        elemrgb(1:3)=(/  0.075148,  0.536560,  0.760502/)
     else if (elem.eq.'La') then
        elemname='Lanthanum'
        elemzz=  57
        elemw=138.905502
        elemcr=  1.870000
        elemvr=  2.790000
        elemrgb(1:3)=(/  0.054107,  0.813106,  0.610205/)
     else if (elem.eq.'Ce') then
        elemname='Cerium'
        elemzz=  58
        elemw=140.115005
        elemcr=  1.830000
        elemvr=  2.730000
        elemrgb(1:3)=(/  0.054107,  0.813106,  0.610205/)
     else if (elem.eq.'Pr') then
        elemname='Praseodymium'
        elemzz=  59
        elemw=140.907593
        elemcr=  1.820000
        elemvr=  2.720000
        elemrgb(1:3)=(/  0.054107,  0.813106,  0.610205/)
     else if (elem.eq.'Nd') then
        elemname='Neodymium'
        elemzz=  60
        elemw=144.240005
        elemcr=  1.810000
        elemvr=  2.700000
        elemrgb(1:3)=(/  0.054107,  0.813106,  0.610205/)
     else if (elem.eq.'Pm') then
        elemname='Promethium'
        elemzz=  61
        elemw=144.912704
        elemcr=  1.800000
        elemvr=  2.690000
        elemrgb(1:3)=(/  0.798076,  0.798076,  0.798076/)
     else if (elem.eq.'Sm') then
        elemname='Samarium'
        elemzz=  62
        elemw=150.360001
        elemcr=  1.800000
        elemvr=  2.690000
        elemrgb(1:3)=(/  0.054107,  0.813106,  0.610205/)
     else if (elem.eq.'Eu') then
        elemname='Europium'
        elemzz=  63
        elemw=151.964996
        elemcr=  1.990000
        elemvr=  2.970000
        elemrgb(1:3)=(/  0.054107,  0.813106,  0.610205/)
     else if (elem.eq.'Gd') then
        elemname='Gadolinium'
        elemzz=  64
        elemw=157.250000
        elemcr=  1.790000
        elemvr=  2.670000
        elemrgb(1:3)=(/  0.054107,  0.813106,  0.610205/)
     else if (elem.eq.'Tb') then
        elemname='Terbium'
        elemzz=  65
        elemw=158.925293
        elemcr=  1.760000
        elemvr=  2.630000
        elemrgb(1:3)=(/  0.054107,  0.813106,  0.610205/)
     else if (elem.eq.'Dy') then
        elemname='Dysprosium'
        elemzz=  66
        elemw=162.500000
        elemcr=  1.750000
        elemvr=  2.610000
        elemrgb(1:3)=(/  0.054107,  0.813106,  0.610205/)
     else if (elem.eq.'Ho') then
        elemname='Holmium'
        elemzz=  67
        elemw=164.930298
        elemcr=  1.740000
        elemvr=  2.600000
        elemrgb(1:3)=(/  0.054107,  0.813106,  0.610205/)
     else if (elem.eq.'Er') then
        elemname='Erbium'
        elemzz=  68
        elemw=167.259995
        elemcr=  1.730000
        elemvr=  2.580000
        elemrgb(1:3)=(/  0.054107,  0.813106,  0.610205/)
     else if (elem.eq.'Tm') then
        elemname='Thulium'
        elemzz=  69
        elemw=168.934204
        elemcr=  1.720000
        elemvr=  2.570000
        elemrgb(1:3)=(/  0.054107,  0.813106,  0.610205/)
     else if (elem.eq.'Yb') then
        elemname='Ytterbium'
        elemzz=  70
        elemw=173.039993
        elemcr=  1.940000
        elemvr=  2.900000
        elemrgb(1:3)=(/  0.054107,  0.813106,  0.610205/)
     else if (elem.eq.'Lu') then
        elemname='Lutetium'
        elemzz=  71
        elemw=174.966995
        elemcr=  1.720000
        elemvr=  2.570000
        elemrgb(1:3)=(/  0.054107,  0.813106,  0.610205/)
     else if (elem.eq.'Hf') then
        elemname='Hafnium'
        elemzz=  72
        elemw=178.490005
        elemcr=  1.570000
        elemvr=  2.340000
        elemrgb(1:3)=(/  0.099196,  0.526039,  0.661306/)
     else if (elem.eq.'Ta') then
        elemname='Tantalum'
        elemzz=  73
        elemw=180.947906
        elemcr=  1.430000
        elemvr=  2.130000
        elemrgb(1:3)=(/  0.099196,  0.526039,  0.661306/)
     else if (elem.eq.'W') then
        elemname='Tungsten'
        elemzz=  74
        elemw=183.839996
        elemcr=  1.370000
        elemvr=  2.040000
        elemrgb(1:3)=(/  0.099196,  0.526039,  0.661306/)
     else if (elem.eq.'Re') then
        elemname='Rhenium'
        elemzz=  75
        elemw=186.207001
        elemcr=  1.350000
        elemvr=  2.010000
        elemrgb(1:3)=(/  0.099196,  0.526039,  0.661306/)
     else if (elem.eq.'Os') then
        elemname='Osmium'
        elemzz=  76
        elemw=190.229996
        elemcr=  1.370000
        elemvr=  2.040000
        elemrgb(1:3)=(/  0.099196,  0.526039,  0.661306/)
     else if (elem.eq.'Ir') then
        elemname='Iridium'
        elemzz=  77
        elemw=192.220001
        elemcr=  1.320000
        elemvr=  1.970000
        elemrgb(1:3)=(/  0.099196,  0.526039,  0.661306/)
     else if (elem.eq.'Pt') then
        elemname='Platinum'
        elemzz=  78
        elemw=195.080002
        elemcr=  1.500000
        elemvr=  1.970000
        elemrgb(1:3)=(/  0.099196,  0.526039,  0.661306/)
     else if (elem.eq.'Au') then
        elemname='Gold'
        elemzz=  79
        elemw=196.966507
        elemcr=  1.500000
        elemvr=  1.850000
        elemrgb(1:3)=(/  0.099196,  0.526039,  0.661306/)
     else if (elem.eq.'Hg') then
        elemname='Mercury'
        elemzz=  80
        elemw=200.589996
        elemcr=  1.700000
        elemvr=  1.900000
        elemrgb(1:3)=(/  0.099196,  0.526039,  0.661306/)
     else if (elem.eq.'Tl') then
        elemname='Thallium'
        elemzz=  81
        elemw=204.383301
        elemcr=  1.550000
        elemvr=  2.310000
        elemrgb(1:3)=(/  0.339671,  0.925829,  0.804088/)
     else if (elem.eq.'Pb') then
        elemname='Lead'
        elemzz=  82
        elemw=207.199997
        elemcr=  1.540000
        elemvr=  2.300000
        elemrgb(1:3)=(/  0.339671,  0.925829,  0.804088/)
     else if (elem.eq.'Bi') then
        elemname='Bismuth'
        elemzz=  83
        elemw=208.980392
        elemcr=  1.540000
        elemvr=  2.300000
        elemrgb(1:3)=(/  0.339671,  0.925829,  0.804088/)
     else if (elem.eq.'Po') then
        elemname='Polonium'
        elemzz=  84
        elemw=208.982407
        elemcr=  1.680000
        elemvr=  2.510000
        elemrgb(1:3)=(/  0.339671,  0.925829,  0.804088/)
     else if (elem.eq.'At') then
        elemname='Astatine'
        elemzz=  85
        elemw=209.987106
        elemcr=  0.000000
        elemvr=  1.270000
        elemrgb(1:3)=(/  0.984970,  0.611708,  0.966409/)
     else if (elem.eq.'Rn') then
        elemname='Radon'
        elemzz=  86
        elemw=222.017593
        elemcr=  0.000000
        elemvr=  1.200000
        elemrgb(1:3)=(/  0.984970,  0.984970,  0.984970/)
     else if (elem.eq.'Fr') then
        elemname='Francium'
        elemzz=  87
        elemw=223.019699
        elemcr=  0.000000
        elemvr=  3.000000
        elemrgb(1:3)=(/  0.450891,  0.790561,  0.910799/)
     else if (elem.eq.'Ra') then
        elemname='Radium'
        elemzz=  88
        elemw=226.025406
        elemcr=  1.900000
        elemvr=  2.840000
        elemrgb(1:3)=(/  0.075148,  0.536560,  0.760502/)
     else if (elem.eq.'Ac') then
        elemname='Actinium'
        elemzz=  89
        elemw=227.027802
        elemcr=  1.880000
        elemvr=  2.810000
        elemrgb(1:3)=(/  0.011122,  0.626738,  0.535057/)
     else if (elem.eq.'Th') then
        elemname='Thorium'
        elemzz=  90
        elemw=232.038101
        elemcr=  1.790000
        elemvr=  2.670000
        elemrgb(1:3)=(/  0.011122,  0.626738,  0.535057/)
     else if (elem.eq.'Pa') then
        elemname='Protactinium'
        elemzz=  91
        elemw=231.035904
        elemcr=  1.610000
        elemvr=  2.400000
        elemrgb(1:3)=(/  0.011122,  0.626738,  0.535057/)
     else if (elem.eq.'U') then
        elemname='Uranium'
        elemzz=  92
        elemw=238.028900
        elemcr=  1.580000
        elemvr=  2.360000
        elemrgb(1:3)=(/  0.011122,  0.626738,  0.535057/)
     else if (elem.eq.'Np') then
        elemname='Neptunium'
        elemzz=  93
        elemw=237.048004
        elemcr=  1.550000
        elemvr=  2.310000
        elemrgb(1:3)=(/  0.798076,  0.798076,  0.798076/)
     else if (elem.eq.'Pu') then
        elemname='Plutonium'
        elemzz=  94
        elemw=244.064194
        elemcr=  1.530000
        elemvr=  2.280000
        elemrgb(1:3)=(/  0.798076,  0.798076,  0.798076/)
     else if (elem.eq.'Am') then
        elemname='Americium'
        elemzz=  95
        elemw=243.061401
        elemcr=  1.510000
        elemvr=  2.250000
        elemrgb(1:3)=(/  0.798076,  0.798076,  0.798076/)
     else if (elem.eq.'Cm') then
        elemname='Curium'
        elemzz=  96
        elemw=247.070297
        elemcr=  0.000000
        elemvr=  2.000000
        elemrgb(1:3)=(/  0.798076,  0.798076,  0.798076/)
     else if (elem.eq.'Bk') then
        elemname='Berkelium'
        elemzz=  97
        elemw=247.070297
        elemcr=  0.000000
        elemvr=  2.000000
        elemrgb(1:3)=(/  0.798076,  0.798076,  0.798076/)
     else if (elem.eq.'Cf') then
        elemname='Californium'
        elemzz=  98
        elemw=251.079605
        elemcr=  0.000000
        elemvr=  2.000000
        elemrgb(1:3)=(/  0.798076,  0.798076,  0.798076/)
     else if (elem.eq.'Es') then
        elemname='Einsteinium'
        elemzz=  99
        elemw=252.082993
        elemcr=  0.000000
        elemvr=  2.000000
        elemrgb(1:3)=(/  0.798076,  0.798076,  0.798076/)
     else if (elem.eq.'Fm') then
        elemname='Fermium'
        elemzz= 100
        elemw=257.095093
        elemcr=  0.000000
        elemvr=  2.000000
        elemrgb(1:3)=(/  0.798076,  0.798076,  0.798076/)
     else if (elem.eq.'Md') then
        elemname='Mendelevium'
        elemzz= 101
        elemw=258.100006
        elemcr=  0.000000
        elemvr=  2.000000
        elemrgb(1:3)=(/  0.798076,  0.798076,  0.798076/)
     else if (elem.eq.'No') then
        elemname='Nobelium'
        elemzz= 102
        elemw=259.100891
        elemcr=  0.000000
        elemvr=  2.000000
        elemrgb(1:3)=(/  0.798076,  0.798076,  0.798076/)
     else if (elem.eq.'Lr') then
        elemname='Lawrencium'
        elemzz= 103
        elemw=262.109985
        elemcr=  0.000000
        elemvr=  2.000000
        elemrgb(1:3)=(/  0.798076,  0.798076,  0.798076/)
     else if (elem.eq.'Rf') then
        elemname='Rutherfordium'
        elemzz= 104
        elemw=261.000000
        elemcr=  0.000000
        elemvr=  2.000000
        elemrgb(1:3)=(/  0.798076,  0.798076,  0.798076/)
     else if (elem.eq.'Db') then
        elemname='Dubnium'
        elemzz= 105
        elemw=262.000000
        elemcr=  0.000000
        elemvr=  2.000000
        elemrgb(1:3)=(/  0.798076,  0.798076,  0.798076/)
     else if (elem.eq.'Sg') then
        elemname='Seaborgium'
        elemzz= 106
        elemw=266.000000
        elemcr=  0.000000
        elemvr=  2.000000
        elemrgb(1:3)=(/  0.798076,  0.798076,  0.798076/)
     else if (elem.eq.'Bh') then
        elemname='Bohrium'
        elemzz= 107
        elemw=264.000000
        elemcr=  0.000000
        elemvr=  2.000000
        elemrgb(1:3)=(/  0.798076,  0.798076,  0.798076/)
     else if (elem.eq.'Hs') then
        elemname='Hassium'
        elemzz= 108
        elemw=269.000000
        elemcr=  0.000000
        elemvr=  2.000000
        elemrgb(1:3)=(/  0.798076,  0.798076,  0.798076/)
     else if (elem.eq.'Mt') then
        elemname='Meitnerium'
        elemzz= 109
        elemw=268.000000
        elemcr=  0.000000
        elemvr=  2.000000
        elemrgb(1:3)=(/  0.798076,  0.798076,  0.798076/)
     else
        elemzz=   0
        elemw=  0.000000
        elemcr=  0.100000
        elemvr=  1.000000
        elemrgb(1:3)=(/  0.000000,  0.000000,  0.984970/)        
     endif

     melem=max(melem,elemzz)
     tmp(elemzz)=tmp(elemzz)+1
     tmpname(elemzz)=elem
     tmpvr(elemzz)=elemvr
     tmpcr(elemzz)=elemcr
     tmprgb(elemzz,1:3)=elemrgb(1:3)     
  enddo

  nspec=0
  do i=0,melem
     if (tmp(i).gt.0) then
        nspec=nspec+1
        specname(nspec)='      '
        specname(nspec)=''''//trim(tmpname(i))//'*'''
        specrad(nspec)=tmpvr(i)*a2b
        specrgb(1:3,nspec)=tmprgb(i,1:3)
        speccr(nspec)=tmpcr(i)*a2b        
     endif
  end do

  nbo=0
  do i=1,nspec
     do j=i,nspec
        nbo=nbo+1
        indb(1,nbo)=i
        indb(2,nbo)=j
        bonddis(nbo)=(speccr(i)+speccr(j))*1.05d0
     enddo
  enddo

end subroutine getelem
