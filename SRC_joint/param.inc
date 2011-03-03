!b     
!     .........PARAMETER INCLUDE FILE.............
!     
!      PARAMETER(NKPT=2000)
!      PARAMETER(NUME=300)  
!     
!      PARAMETER(INUME=5000)
!      PARAMETER(INUMEDEN=50)
!      PARAMETER(MET=10000)
      PARAMETER(MG0=9)
!      PARAMETER(NSW=9) 
!..new struct + latgen 
!      PARAMETER          (NATO=   18)
!      PARAMETER          (NDIF=   64)
      PARAMETER(NSYM= 48)                                              
!     NATO        MAX NUMBER OF INEQ: ATOMS
!
!     NDIF        MAX NUMBER OF ATOMPOSITION IN CELL
!
!     NSYM        MAX NUMBER OF SYM-OPERATION
!
!     NUME        MAX NUMBER OF EIGENVALUES PER K - POINT 
!
!     INUME       MAX NUMBER OF BANDCOMBINATIONS nume*(nume+1)/2 
!
!     INUMEden    MAX NUMBER OF BANDCOMBINATIONS for band analysis     
!
!     MET         MAXIMUM FOR ENERGYMESH
!
!     MG          MAXIMUM NUMBER OF SPECTRA 
!                (polarizations: see input case.inop) 
!
!     NSW         NUMBER OF SWITCHCODES   fixed!
!
!....special: not used in standart application
      PARAMETER   (NEKPT=5)   
      PARAMETER   (NUMTETR= 6)
!     NEKPT       NUMBER OF NOT EQUIVALENT K_POINTS   fixed!
!
!     NUMTETR     NUMBER OF TETRAEDERS        fixed!
