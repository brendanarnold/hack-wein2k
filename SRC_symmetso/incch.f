      subroutine incch(ifr,ifw,nat,ieq)
! this program serves to read input files in FLAPW WIEN program
! then to construct the input files 
! when s-o coupling and magnetization are present
       implicit real*8(a-h,o-z)
       character b
        dimension b(100,60)
        integer ieq(*)
        do 80 jatom=1,nat
        read(ifr,*)n,shift
         do j=1,n
         read(ifr,1980)(b(j,i),i=1,37)
         enddo
2000   format(i2,f5.2,'     NUMBER OF ORBITALS (EXCLUDING SPIN)')
        do 81 jeq=1,ieq(jatom)
        write(ifw,2000)n,shift
         do j=1,n
         write(ifw,1980)(b(j,i),i=1,37)
         enddo
 81     continue
 80     continue
        n=0
        write(ifw,100)n
100     format(i2)
 1980   format(80a1)
        return
        end
