      SUBROUTINE FOMAI3(forout,forcea,ftot,fvdrho,fsph,fnsp, &
                        fpsi,fequ,fsph2)                           
!                                                                       
!     output of forces
!
      USE param
      use defs
      use struk; USE com
      IMPLICIT REAL*8 (A-H,O-Z)

      logical     forout,forcea(0:3,*)
      real*8      ftot(0:3,ndif),fvdrho(0:3,ndif) &
        ,fsph(0:3,ndif),fnsp(0:3,ndif) &
        ,fpsi(0:3,ndif),fequ(0:3),fsph2(0:3,ndif)

        ia2=0 
        ia=0
        DO 5023 ja=1,nat
          
          mu=0
          ia1=ia2+1
          ia2=ia2+mult(ja)
          DO 5024 ia=ia1,ia2
            mu=mu+1
            DO 5025 ik=1,3
            if (forcea(ik,ja)) then
               fequ(ik)=( &
                  +fsph(ik,ia) &
                  +fsph2(ik,ia) &
                  +fnsp(ik,ia))/mult(ja)
               fpsi(ik,ja)=fpsi(ik,ja)+fequ(ik)
            else 
               fpsi(ik,ja)=0.0d0
            endif
 5025       CONTINUE
            call mag(fequ(0))   
            write(6,77) ja,mu,'+ EQU',(fequ(ik),ik=0,3)
!            if (forout) then
              call mag(fnsp(0,ia))   
              call mag(fsph(0,ia))   
              call mag(fsph2(0,ia))   
              write(6,77) ja,mu,'+ SPH',(fsph(ik,ia),ik=0,3)
              write(6,77) ja,mu,'+ SP2',(fsph2(ik,ia),ik=0,3)
              write(6,77) ja,mu,'+ NSP',(fnsp(ik,ia),ik=0,3)
!              write(6,77) ja,mu,'+ SUR',(fsur(ik,ia),ik=0,3)
              write(6,77) ja,mu,'> EQU',(fequ(ik),ik=0,3)
!            endif
            
 5024     CONTINUE
          
!     calculate partial forces 
          DO 5026 ik=1,3
            ftot(ik,ja)=fvdrho(ik,ja)+fpsi(ik,ja)
 5026     CONTINUE

          call mag(fvdrho(0,ja))
          call mag(fpsi(0,ja))
          call mag(ftot(0,ja))
        
 77   format(2i3,a7,4e15.7)
 78   FORMAT (7x,'VALENCE-FORCE IN mRy/a.u. = |F|',3x,'Fx',13x, &
              'Fy',13x,'Fz')
 79   FORMAT (':FVA',i3.3,':',1x,i3,'.ATOM',4f15.3)
          write(6,77) ja,1,'= PSI',(fpsi(ik,ja),ik=0,3)
          write(6,77) ja,1,'+ VDR',(fvdrho(ik,ja),ik=0,3)
          write(6,77) ja,1,'> TOT',(ftot(ik,ja),ik=0,3)
          write(21,78) 
          write(21,79) ja,ja,(ftot(ik,ja)*1000,ik=0,3)
          write(6,*)
          
 5023 CONTINUE
      RETURN
      END                                                               
