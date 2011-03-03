


      PROGRAM setelast
      CHARACTER*1 choi,lat,cho
      INTEGER numb,i
      CHARACTER*79 row
      REAL*8 var(100),a0


      OPEN(unit=10,file='init.struct',status='old',err=8000)
      READ(10,1000) row
      READ(10,1010) lat
      READ(10,1000) row
      READ(10,1020) a0
      CLOSE(10)

 5    WRITE(*,*) 'Do you want to create input files for EOS'
      WRITE(*,*) 'calculations (i.e. simply vary volume) (y/n) ?'
      READ(*,"(A1)",err=5) cho

      IF (cho.eq.'y') then      


 10   WRITE(*,*) ' '
      WRITE(*,*) '**** Setup for EOS calculations'
      WRITE(*,*) ' '
      WRITE(*,*) ' '
      WRITE(*,*) 'Do you want to use default parameters ?(y/n)'
      WRITE(*,*) '     (-10%,-5%,0,5%,10%)'
      READ(*,"(A1)",err=10) choi
      
      IF (choi.eq.'n') then
         WRITE(*,*) 'Number of structure changes(max=99)?'
         READ(*,"(I2)",err=10) numb
         DO i=1,numb
            WRITE(*,*) 'Enter value ',i,' in %'
            READ(*,*) var(i)
         ENDDO

      ELSEIF (choi.eq.'y') then
         numb=5
         var(1)=-10.d0
         var(2)=-5.d0
         var(3)=0.d0
         var(4)=5.d0
         var(5)=10.d0
      ELSE 
         GOTO 10
      ENDIF

      CALL creatfeos(numb,var,a0)

      ENDIF

 49   WRITE(*,*) 'Do you want to create input files for tetragonal'
      WRITE(*,*) 'strain calculations (i.e. vary c/a ratio) (y/n) ?'
      READ(*,"(A1)",err=49) cho

      IF (cho.eq.'y') then      



 50   WRITE(*,*) ' '
      WRITE(*,*) '**** Setup for tetragonal calculations'
      WRITE(*,*) ' '
      WRITE(*,*) ' '
      WRITE(*,*) 'Do you want to use default parameters ?(y/n)'
      WRITE(*,*) '     (-4%,-2%,0,2%,4%)'
      READ(*,"(A1)",err=50) choi
      
      IF (choi.eq.'n') then
         WRITE(*,*) 'Number of structure changes(max=99)?'
         READ(*,"(I2)",err=50) numb
         WRITE(*,*) ' '
         WRITE(*,*) '******** 0.0 value MUST be calculated for analysis'
         WRITE(*,*) ' '
         DO i=1,numb
            WRITE(*,*) 'Enter value ',i,' in %'
            READ(*,*) var(i)
         ENDDO

      ELSEIF (choi.eq.'y') then
         numb=5
         var(1)=-4.d0
         var(2)=-2.d0
         var(3)=0.d0
         var(4)=2.d0
         var(5)=4.d0
      ELSE 
         GOTO 50
      ENDIF

      CALL creatftet(numb,var,a0)

      ENDIF

 89   WRITE(*,*) 'Do you want to create input files for rhombohedral'
      WRITE(*,*) 'strain calculations (i.e. vary a+b+c length) (y/n) ?'
      READ(*,"(A1)",err=89) cho

      IF (cho.eq.'y') then      

 90   WRITE(*,*) ' '
      WRITE(*,*) '**** Setup for rhombohedral calculations'
      WRITE(*,*) ' '
      WRITE(*,*) ' '
      WRITE(*,*) 'Do you want to use default parameters ?(y/n)'
      WRITE(*,*) '     (-4%,-2%,0,2%,4%)'
      READ(*,"(A1)",err=90) choi
      
      IF (choi.eq.'n') then
         WRITE(*,*) 'Number of structure changes(max=99)?'
         READ(*,"(I2)",err=90) numb
         WRITE(*,*) ' '
         WRITE(*,*) '******** 0.0 value MUST be calculated for analysis'
         WRITE(*,*) ' '
         DO i=1,numb
            WRITE(*,*) 'Enter value ',i,' in %'
            READ(*,*) var(i)
         ENDDO

      ELSEIF (choi.eq.'y') then
         numb=5
         var(1)=-4.d0
         var(2)=-2.d0
         var(3)=0.d0
         var(4)=2.d0
         var(5)=4.d0
      ELSE 
         GOTO 90
      ENDIF

      CALL creatfrh(numb,var,a0,lat)

      ENDIF


      OPEN(unit=10,file='eos.job')
      WRITE(10,1050)'#!/bin/csh -f'
      WRITE(10,1050)'#Modify this script according to your needs'
      WRITE(10,1050)' '
      WRITE(10,1050)'set flist = `ls eos_*.struct | cut -c 1-9`'
      WRITE(10,1050)'cd ./eos'
      WRITE(10,1050)'foreach i ($flist)'
      WRITE(10,1050)'echo $i'
      WRITE(10,1050)'cp ../$i.struct ./eos.struct'
      WRITE(10,1050)'x_lapw dstart'
      WRITE(10,1050)'#x_lapw dstart -up'
      WRITE(10,1050)'#x_lapw dstart -dn'
      WRITE(10,1050)'#cp ../result/$i.clmsum ./eos.clmsum'
      WRITE(10,1050)'#cp ../result/$i.clmup ./eos.clmup'
      WRITE(10,1050)'#cp ../result/$i.clmdn ./eos.clmdn'
      WRITE(10,1050)'run_lapw -ec 0.0001'
      WRITE(10,1050)'	set stat = $status'
      WRITE(10,1050)'	if ($stat) then'
      WRITE(10,1050)'	 echo "ERROR status in" $i'
      WRITE(10,1050)'	 exit 1'
      WRITE(10,1050)'	endif'
      WRITE(10,1050)' '
      WRITE(10,1050)'#echo $i >> error'
      WRITE(10,1050)'#x lapw2 -p -qtl | & tee -a error'
      WRITE(10,1050)'#x tetra'
      WRITE(10,1050)'#mv eos.outputt $i.outputt'
      WRITE(10,1050)'#mv eos.qtl $i.qtl'
      WRITE(10,1050)'#mv eos.dos1 $i.dos1'
      WRITE(10,1050)'#mv eos.dos1ev $i.dos1ev'
      WRITE(10,1050)'save_lapw $i'
      WRITE(10,1050)'mv $i.* ../result'
      WRITE(10,1050)'end'

      CLOSE(10)

      OPEN(unit=20,file='tetra.job')
      WRITE(20,'(a)')'#!/bin/csh -f'
      WRITE(20,1050)'#Modify this script according to your needs'
      WRITE(20,1050)' '
      WRITE(20,1050)'set flist = `ls tetra_*.struct | cut -c 1-11`'
      WRITE(20,1050)'cd ./tetra'
      WRITE(20,1050)'foreach i ($flist)'
      WRITE(20,1050)'echo $i'
      WRITE(20,1050)'cp ../$i.struct ./tetra.struct'
      WRITE(20,1050)'x_lapw dstart'
      WRITE(20,1050)'#x_lapw dstart -up'
      WRITE(20,1050)'#x_lapw dstart -dn'
      WRITE(20,1050)'#cp ../result/$i.clmsum ./tetra.clmsum'
      WRITE(20,1050)'#cp ../result/$i.clmup ./tetra.clmup'
      WRITE(20,1050)'#cp ../result/$i.clmdn ./tetra.clmdn'
      WRITE(20,1050)'run_lapw -ec 0.0001'
      WRITE(20,1050)'	set stat = $status'
      WRITE(20,1050)'	if ($stat) then'
      WRITE(20,1050)'	    echo "ERROR status in" $i'
      WRITE(20,1050)'	    exit 1'
      WRITE(20,1050)'	endif'
      WRITE(20,1050)' '
      WRITE(20,1050)'#echo $i >> error'
      WRITE(20,1050)'#x lapw2 -p -qtl | & tee -a error'
      WRITE(20,1050)'#x tetra'
      WRITE(20,1050)'#mv tetra.outputt $i.outputt'
      WRITE(20,1050)'#mv tetra.qtl $i.qtl'
      WRITE(20,1050)'#mv tetra.dos1 $i.dos1'
      WRITE(20,1050)'#mv tetra.dos1ev $i.dos1ev'
      WRITE(20,1050)'save_lapw $i'
      WRITE(20,1050)'mv $i.* ../result'
      WRITE(20,1050)'end'

      CLOSE(20)
      
      OPEN(unit=30,file='rhomb.job')
      WRITE(30,'(a)')'#!/bin/csh -f'
      WRITE(30,1050)'#Modify this script according to your needs'
      WRITE(30,1050)'set flist = `ls rhomb_*.struct | cut -c 1-11`'
      WRITE(30,1050)'cd ./rhomb'
      WRITE(30,1050)'foreach i ($flist)'
      WRITE(30,1050)'echo $i'
      WRITE(30,1050)'cp ../$i.struct ./rhomb.struct'
      WRITE(30,1050)'x_lapw dstart'
      WRITE(30,1050)'#x_lapw dstart -up'
      WRITE(30,1050)'#x_lapw dstart -dn'
      WRITE(30,1050)'#cp ../result/$i.clmsum ./rhomb.clmsum'
      WRITE(30,1050)'#cp ../result/$i.clmup ./rhomb.clmup'
      WRITE(30,1050)'#cp ../result/$i.clmdn ./rhomb.clmdn'
      WRITE(30,1050)'run_lapw -ec 0.0001'
      WRITE(30,1050)'	set stat = $status'
      WRITE(30,1050)'	if ($stat) then'
      WRITE(30,1050)'	    echo "ERROR status in" $i'
      WRITE(30,1050)'	    exit 1'
      WRITE(30,1050)'	endif'
      WRITE(30,1050)' '
      WRITE(30,1050)'#echo $i >> error'
      WRITE(30,1050)'#x lapw2 -p -qtl | & tee -a error'
      WRITE(30,1050)'#x tetra'
      WRITE(30,1050)'#mv rhomb.outputt $i.outputt'
      WRITE(30,1050)'#mv rhomb.qtl $i.qtl'
      WRITE(30,1050)'#mv rhomb.dos1 $i.dos1'
      WRITE(30,1050)'#mv rhomb.dos1ev $i.dos1ev'
      WRITE(30,1050)'save_lapw $i'
      WRITE(30,1050)'mv $i.* ../result'
      WRITE(30,1050)'end'

      CLOSE(30)

 1000 FORMAT(A79)
 1010 FORMAT(A1)
 1020 FORMAT(F10.6)
 1050 FORMAT(a)

      GOTO 8001
 8000 write(*,*) 'No valid file init.struct'

      STOP

 8001 write(*,*) 'Structure files generated...'



      END


      SUBROUTINE creatfeos(numb,var,a0)
      INTEGER numb,i,j
      REAL*8 var(100),a0,a,eps
      REAL*8 at,alp
      CHARACTER*79 row
      CHARACTER*16 namout
      CHARACTER*5 indic
      
      DO i=1,numb
         eps=var(i)/100.d0
         WRITE(indic,"(F5.1)") var(i)
         namout='eos_'//indic//'.struct'
         DO j=1,16
            IF (namout(j:j).eq.' ') namout(j:j)= '_'
         ENDDO
         a=a0*((1.d0+eps)**(1.d0/3.d0))

         OPEN(unit=10,file='eos.templ',status='old',err=3000)
         OPEN(unit=20,file=namout)
         DO j=1,3
            READ(10,1000) row
            WRITE(20,1000) row
         ENDDO
         READ(10,1020) at,at,at,alp,alp,alp
         WRITE(20,1020) a,a,a,alp,alp,alp
         
         DO j=1,10000
            READ(10,1000,END=500) row
            WRITE(20,1000) row       
         ENDDO
 500     CONTINUE

         CLOSE(20)
         CLOSE(10)

      ENDDO

 1000 FORMAT(A79)
 1020 FORMAT(6F10.6)

      GOTO 3001
 3000 WRITE(*,*) 'No valid eos.templ file'
      STOP

 3001 WRITE(*,*) numb,' eos_xxx.x.struct files generated' 

      RETURN
      END



      SUBROUTINE creatftet(numb,var,a0)
      INTEGER numb,i,j
      REAL*8 var(100),a0,a,c,at,ct,eps,alp
      CHARACTER*18 namout
      CHARACTER*5 indic
      CHARACTER*79 row
      EXTERNAL c2te

      DO i=1,numb
         eps=var(i)/100.d0
         WRITE(indic,"(F5.1)") var(i)
         namout='tetra_'//indic//'.struct'
         DO j=1,18
            IF (namout(j:j).eq.' ') namout(j:j)= '_'
         ENDDO
         
         CALL c2te(a0,eps,at,ct)

         OPEN(unit=10,file='tetra.templ',status='old',err=3000)
         OPEN(unit=20,file=namout)
         DO j=1,3
            READ(10,1000) row
            WRITE(20,1000) row
         ENDDO
         READ(10,1020) a,a,c,alp,alp,alp
         WRITE(20,1020) at,at,ct,alp,alp,alp
         
         DO j=1,10000
            READ(10,1000,END=500) row
            WRITE(20,1000) row       
         ENDDO
 500     CONTINUE

         CLOSE(20)
         CLOSE(10)

      ENDDO

 1000 FORMAT(A79)
 1020 FORMAT(6F10.6)

      GOTO 3001
 3000 WRITE(*,*) 'No valid tetra.templ file'
      STOP

 3001 WRITE(*,*) numb,' tetra_xxx.x.struct files generated' 


      RETURN
      END


      SUBROUTINE creatfrh(numb,var,a0,lat)
      INTEGER numb,i,j
      REAL*8 var(100),a0,eps,ah,ch,a,c,alp
      CHARACTER*1 lat
      CHARACTER*18 namout
      CHARACTER*5 indic
      CHARACTER*79 row
      EXTERNAL c2rh

      DO i=1,numb
         eps=var(i)/100.d0
         WRITE(indic,"(F5.1)") var(i)
         namout='rhomb_'//indic//'.struct'
         DO j=1,18
            IF (namout(j:j).eq.' ') namout(j:j)= '_'
         ENDDO
         
         CALL c2rh(a0,eps,lat,ah,ch)

         OPEN(unit=10,file='rhomb.templ',status='old',err=3000)
         OPEN(unit=20,file=namout)
         DO j=1,3
            READ(10,1000) row
            WRITE(20,1000) row
         ENDDO
         READ(10,1020) a,a,c,alp,alp,alp
         WRITE(20,1020) ah,ah,ch,alp,alp,alp
         
         DO j=1,10000
            READ(10,1000,END=500) row
            WRITE(20,1000) row       
         ENDDO
 500     CONTINUE

         CLOSE(20)
         CLOSE(10)

      ENDDO

 1000 FORMAT(A79)
 1020 FORMAT(6F10.6)

      GOTO 3001
 3000 WRITE(*,*) 'No valid rhomb.templ file'
      STOP

 3001 WRITE(*,*) numb,' rhomb_xxx.x.struct files generated' 

      RETURN
      END



