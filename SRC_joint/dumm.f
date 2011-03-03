	subroutine dumm(n,an)
        character*2 an
        write(22,'(i2)') n
	rewind 22
	read(22,'(a2)') an
        rewind 22
	return 
	end
