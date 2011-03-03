            subroutine workf1(kzz,potk,vkvxc,nk) 
        
            IMPLICIT REAL*8 (A,B,D-H,O-Z) 
!            INCLUDE 'param.inc' 
        
            COMPLEX*16       POTK,vkvxc(nk,2) 
            dimension        KZZ(3,Nk),POTK(Nk) 
!      c 
!      c      calculates electrostatic potential at z=0. 
!      c 
!      c      \frac{1}{\omega_{surf}} \int V( r_x, r_y, r_z=0 ) d^2r = 
!      c      \sum_{G_z} V( G_x=0, G_y=0, G_z) 
!      c 
!      c      attention: valid only for PC lattice + inversion 
!      c      geometry !!! 
            sup=0. 
            suv=0. 
            sup5=0. 
            suv5=0. 

            supy=0.
            suvy=0.
            sup5y=0.
            suv5y=0.
            supx=0
            suvx=0
            sup5x=0
            suv5x=0

            do ik=1,nk 
               ix=kzz(1,ik) 
               iy=kzz(2,ik) 
               iz=kzz(3,ik) 
               if (ix.eq.0.and.iy.eq.0) then 
!      c 
!      c...........star is not always multipied !! 
!                  if (iz.eq.0) then 
!                     suc1=dble(cvalue(ik)) 
!                     sup1=dble(potk(ik)) 
!                     suv1=dble(vkvxc(ik,1)) 
!                  endif 
!      c 
                  if (iz.eq.(iz/2)*2) then
                    phs=1.d0
                  else
                    phs=-1.0 
                  endif
!                  suc=suc+dble(cvalue(ik)) 
                  sup=sup+dble(potk(ik))
                  suv=suv+dble(vkvxc(ik,1))
                  sup5=sup5+dble(potk(ik))*phs 
                  suv5=suv5+dble(vkvxc(ik,1))*phs 
               endif 
!       Added section for y (CXZ case)
               if (ix.eq.0.and.iz.eq.0) then 
                  if (iy.eq.(iy/2)*2) then
                    phsy=1.d0
                  else
                    phsy=-1.0 
                  endif
!                  sucy=sucy+dble(cvalue(ik)) 
                  supy=supy+dble(potk(ik))
!                  print*, iy,supy,potk(ik)
                  suvy=suvy+dble(vkvxc(ik,1))
                  sup5y=sup5y+dble(potk(ik))*phsy 
                  suv5y=suv5y+dble(vkvxc(ik,1))*phsy 
               endif 
!       Added section for x
               if (iz.eq.0.and.iy.eq.0) then
                  if (ix.eq.(ix/2)*2) then
                    phsx=1.d0
                  else
                    phsx=-1.0
                  endif
!                  sucx=sucx+dble(cvalue(ik))
                  supx=supx+dble(potk(ik))
                  suvx=suvx+dble(vkvxc(ik,1))
                  sup5x=sup5x+dble(potk(ik))*phsx
                  suv5x=suv5x+dble(vkvxc(ik,1))*phsx
               endif

        
            enddo 
!      c 
!      c.....cvalue: star is not mutiplied; therefore the Gz=0 component 
!      c             has to be substracted once 
!!!               suc=2.*suc -suc1 
!      c 
!      c.....potk: is already muliplied by star 
               sup=sup 
!      c 
!      c        sup=2.*sup -sup1 
!      c        suv=2.*suv -suv1 
                suc=sup+suv
                suc5=sup5+suv5

                sucy=supy+suvy
                suc5y=sup5y+suv5y
                sucx=supx+suvx
                suc5x=sup5x+suv5x
               write(6,100) sup,sup5 
               write(21,100) sup,sup5
               write(6,102)supy,sup5y
               write(21,102)supy,sup5y
               write(6,101) suc,sup,suv ,suc5,sup5,suv5
               write(21,101) suc,sup,suv, suc5,sup5,suv5
               write(6,103) sucy,supy,suvy,suc5y,sup5y,suv5y
               write(21,103)sucy,supy,suvy,suc5y,sup5y,suv5y
               write(6,104) sucx,supx,suvx,suc5x,sup5x,suv5x
               write(21,104)sucx,supx,suvx,suc5x,sup5x,suv5x

  100          format(7x,'ELS_POTENTIAL_AT Z=0 and Z=0.5:',2f10.5) 
  102          format(7x,'ELS_POTENTIAL_AT Y=0 and Y=0.5:',2f10.5) 
  101          format(':VZERO:v0,v0c,v0x',3f10.5,' v5,v5c,v5x',3f10.5) 
  103          format(':VZERY:v0,v0c,v0x',3f10.5,' v5,v5c,v5x',3f10.5)
  104          format(':VZERX:v0,v0c,v0x',3f10.5,' v5,v5c,v5x',3f10.5)

               return 
               end
