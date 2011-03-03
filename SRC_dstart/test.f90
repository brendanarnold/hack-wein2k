	module testallocate
	double precision, pointer :: testfeld(:)
	double precision, pointer :: testfeld2(:)
	contains
	subroutine doreallocate(tf, newdimension)
	double precision, pointer :: hilfsfeld(:)
	double precision, pointer :: tf(:)
	allocate(hilfsfeld(newdimension))
	hilfsfeld=tf
	deallocate(tf)
	allocate(tf(newdimension))
	tf=hilfsfeld
	end subroutine doreallocate
	end module testallocate

	program testprogram
	use testallocate

	allocate(testfeld(5))
	allocate(testfeld2(5))
	do i=1,5
	  testfeld(i)=i
	  testfeld2(i)=2*i
	end do
	call doreallocate(testfeld,8)
	call doreallocate(testfeld2,3)
	write(*,*) testfeld
	write(*,*)
	write(*,*) testfeld2
	end
