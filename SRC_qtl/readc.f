        SUBROUTINE READC(MODUS,iso)
        USE param
        USE case
        implicit real*8 (a-h,o-z)

!  variable BLEND detects end of the block for SUMA calculation
!  (i.g. mark blocks belonging to J=1/2 or 3/2), ISUM(0:ndim2,nato)
!  controls the devision of CF into blocks: ISUM(0)...number of blocks
!  ISUM(n)...last line of block n

	character*4 MODUS
	character*200 aline
	complex*16 imag,czero,sum,tt(ndim2,ndim2)
	dimension a(ndim2),b(ndim2) 

      DATA IMAG/(0D0,1D0)/,CZERO/(0D0,0D0)/

       do 200 i=1,natom
	nb=0
	ifile=50+i
	n=2*lcase(i)+1

        do 150 j=1,iso*n
        read(ifile,5)aline
	write(90,5)aline(2:200)
 5      format(200A)
	close(90)

        if (aline(1:1).eq.'*') then 
        nb=nb+1
        isum(nb,i)=j
        end if

 	read(90,*)(a(m),b(m),m=1,iso*n)
	rewind(90)

	do m=1,2*n
	cf(j,m,i)=a(m)+imag*b(m)
	end do

 150   continue
	isum(0,i)=nb


 200   continue

	if (MODUS.eq.'SUMA') then
	write(6,*)'BLOCKING SUMMARY:'
	do i=1,natom
	write(6,*)'atom:',i,'   number of blocks:',isum(0,i)
	write(6,*)'end lines of blocks:'
	write(6,*)(isum(j,i),j=1,isum(0,i))
	end do
	end if

        write(6,*)'reading transformation Y(l,m)f(s)->D(l,j)'
	do ia=1,natom
!..... unitarity test
        n=iso*(2*lcase(ia)+1)
        do i=1,n
        do j=1,n
        sum=czero
        do k=1,n
        sum=sum+cf(i,k,1)*conjg(cf(j,k,1))
        end do
        tt(i,j)=sum
        end do
        end do

	write(6,*)'ATOM NUM:',IATOM(ia)
        write(6,*)'unitarity test:'
        write(6,*)'REAL:'
        do i=1,n
        write(6,357)(dble(tt(i,j)),j=1,n)
        end do
        write(6,*)'IMAG:'
        do i=1,n
        write(6,357)(aimag(tt(i,j)),j=1,n)
        end do
	end do
 357  format(14f10.6)
	end 
