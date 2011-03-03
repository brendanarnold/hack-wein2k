	subroutine cputim(cp)
	real*8 cp
	integer tms (4)
	integer utime,stime
	equivalence (tms(1),utime)
	equivalence (tms(2),stime)
	integer cutime,cstime
	equivalence (tms(3),cutime)
	equivalence (tms(4),cstime)

	integer HZ
	parameter (HZ = 100)

	call times (tms)
	cp = dble (utime+stime) / HZ
	return
	end

