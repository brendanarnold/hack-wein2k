        subroutine filename(string,length,fname)	
	character*80 fname,string
        do i=1,80
        if(string(i:i).eq.'.') then
        length=i-1
        goto 10
        endif
        enddo
 10     continue
        fname=string(1:length)
        return
        end
