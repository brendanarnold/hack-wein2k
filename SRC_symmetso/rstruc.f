      subroutine rstruc(title,lattic,nat,a,alpha,index2, &
       iatom,index1,pos,mult,name,iatnr,tm)
!-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      CHARACTER*4       LATTIC                                          
      CHARACTER*79      TITLE,NAME                                            
       dimension pos(3,*),index1(*),index2(*),tm(3)
       dimension name(*),mult(*),a(3),iatom(*)
       dimension alpha(3),br2(3,3)
!      common /iat/ iatnr(NATO)
      dimension iatnr(*)
!-----------------------------------------------------------------------
      pi=4.d0*atan(1.d0)
!.....START READING FILE STRUCT                                           
      READ(24,1510) TITLE 
      WRITE(6,1510) TITLE                                               
      READ(24,1511) LATTIC,NAT                                          
 17   FORMAT(6F10.7)                                                    
 1011 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
 1012 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,7X,I2,8X,I2)     
 1020 FORMAT(///,3X,'ERROR IN EBND  : MULT(JATOM)=0 ...',                &
             /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)           
 1050 FORMAT(A79)
 1509 FORMAT(///,'DETERMINATION OF POINTGROUP FOR ALL POSITIONS')
 1510 FORMAT(A79)                                                       
 1511 FORMAT(A4,23X,I3,/,13X,A4)                                        
!      IF(NAT.GT.NATO) STOP 'NATO TOO SMALL'                             
!     READ IN LATTICE CONSTANTS                                         
      READ(24,17) A(1),A(2),A(3),alpha(1),alpha(2),alpha(3)   
      if(alpha(1).eq.0.d0) alpha(1)=90.d0               
      if(alpha(2).eq.0.d0) alpha(2)=90.d0               
      if(alpha(3).eq.0.d0) alpha(3)=90.d0               
      if(lattic(1:1).eq.'H'.or.lattic(1:1).eq.'R') alpha(3)=90.d0
!     INDEX COUNTS ALL ATOMS IN THE UNIT CELL, JATOM COUNTS ONLY THE    
!     NONEQUIVALENT ATOMS                                               
      INDEX=0                                                           
      if(lattic(1:1).eq.'R') then
         br2(1,1)=1.d0/3.d0
         br2(2,1)=-1.d0/3.d0
         br2(3,1)=-1.d0/3.d0
         br2(1,2)=1.d0/3.d0
         br2(2,2)=2.d0/3.d0
         br2(3,2)=-1.d0/3.d0
         br2(1,3)=-2.d0/3.d0
         br2(2,3)=-1.d0/3.d0
         br2(3,3)=-1.d0/3.d0
                 TAU1=tm(1)
                 TAU2=tm(2)
                 TAU3=tm(3)
                 tm(1)=Br2(1,1)*TAU1+Br2(1,2)*TAU2+Br2(1,3)*TAU3
                 tm(2)=Br2(2,1)*TAU1+Br2(2,2)*TAU2+Br2(2,3)*TAU3
                 tm(3)=Br2(3,1)*TAU1+Br2(3,2)*TAU2+Br2(3,3)*TAU3
!         br2(1,1)=-1.d0/3.d0
!         br2(2,1)=1.d0/3.d0
!         br2(3,1)=1.d0/3.d0
!         br2(1,2)=2.d0/3.d0
!         br2(2,2)=1.d0/3.d0
!         br2(3,2)=1.d0/3.d0
!         br2(1,3)=-1.d0/3.d0
!         br2(2,3)=-2.d0/3.d0
!         br2(3,3)=1.d0/3.d0
         a(3)=1.d0
         a(1)=1.d0
         a(2)=1.d0
      end if
      DO 50 JATOM = 1,NAT                                               
         INDEX=INDEX+1
         iatom(index)=jatom               
!         IF(INDEX.GT.NDIF)then
!         write(6,*)' RSTRUC 50 index',index,' NDIF',NDIF
!         STOP 'NDIF TOO SMALL'                        
!         endif
         index1(jatom)=index
         READ(24,1012) IATNR(jatom),( POS(J,INDEX),J=1,3 ),MULT(JATOM)  
               if(lattic(1:1).eq.'P'.or.lattic(1:1).eq.'R'.or. &
                  lattic(1:1).eq.'H') then
                 index=index+1
                 iatom(index)=jatom                                           
                 pos(1,index)=pos(1,index-1)+1.0d0         
                 pos(2,index)=pos(2,index-1)         
                 pos(3,index)=pos(3,index-1)    
!                 if(pos(1,index).ge.1.0001) pos(1,index)=pos(1,index)-1.0
                 index=index+1
                 iatom(index)=jatom                       
                 pos(1,index)=pos(1,index-1)-1.0d0         
                 pos(2,index)=pos(2,index-1)+1.0d0         
                 pos(3,index)=pos(3,index-1)    
!                 if(pos(2,index).ge.1.0001) pos(2,index)=pos(2,index)-1.0
                 index=index+1
                 iatom(index)=jatom                                           
                 pos(1,index)=pos(1,index-1)
                 pos(2,index)=pos(2,index-1)-1.0d0         
                 pos(3,index)=pos(3,index-1)+1.0d0    
!                 if(pos(3,index).ge.1.0001) pos(3,index)=pos(3,index)-1.0
                   if(lattic(1:1).eq.'R') then
                   index=index+1
                   iatom(index)=jatom                                          
                   pos(1,index)=pos(1,index-1)
                   pos(2,index)=pos(2,index-1)-1.0d0         
                   pos(3,index)=pos(3,index-1)-1.0d0    
!                   if(pos(3,index).ge.1.0001) pos(3,index)=pos(3,index)-1.0
                   endif
               endif
               if(lattic(1:3).eq.'CXY') then
                 index=index+1
                 iatom(index)=jatom                                           
                 pos(1,index)=pos(1,index-1)+0.5d0         
                 pos(2,index)=pos(2,index-1)+0.5d0         
                 pos(3,index)=pos(3,index-1)    
                 if(pos(1,index).ge.0.9999) pos(1,index)=pos(1,index)-1.0
                 if(pos(2,index).ge.0.9999) pos(2,index)=pos(2,index)-1.0
               endif
               if(lattic(1:3).eq.'CXZ') then
                 index=index+1
                 iatom(index)=jatom                                            
                 pos(1,index)=pos(1,index-1)+0.5d0         
                 pos(2,index)=pos(2,index-1)             
                 pos(3,index)=pos(3,index-1)+0.5d0
                 if(pos(1,index).ge.0.9999) pos(1,index)=pos(1,index)-1.0
                 if(pos(3,index).ge.0.9999) pos(3,index)=pos(3,index)-1.0
               endif
               if(lattic(1:3).eq.'CYZ') then
                 index=index+1
                 iatom(index)=jatom                                            
                 pos(1,index)=pos(1,index-1)             
                 pos(2,index)=pos(2,index-1)+0.5d0         
                 pos(3,index)=pos(3,index-1)+0.5d0
                 if(pos(2,index).ge.0.9999) pos(2,index)=pos(2,index)-1.0
                 if(pos(3,index).ge.0.9999) pos(3,index)=pos(3,index)-1.0
               endif
               if(lattic(1:1).eq.'B') then
                 index=index+1
                 iatom(index)=jatom                                            
                 pos(1,index)=pos(1,index-1)+0.5d0         
                 pos(2,index)=pos(2,index-1)+0.5d0         
                 pos(3,index)=pos(3,index-1)+0.5d0
                 if(pos(1,index).ge.0.9999) pos(1,index)=pos(1,index)-1.0
                 if(pos(2,index).ge.0.9999) pos(2,index)=pos(2,index)-1.0
                 if(pos(3,index).ge.0.9999) pos(3,index)=pos(3,index)-1.0
               endif
               if(lattic(1:1).eq.'F') then
                 index=index+1
                 iatom(index)=jatom                                            
                 pos(1,index)=pos(1,index-1)+0.5d0         
                 pos(2,index)=pos(2,index-1)+0.5d0         
                 pos(3,index)=pos(3,index-1)        
                 if(pos(1,index).ge.0.9999) pos(1,index)=pos(1,index)-1.0
                 if(pos(2,index).ge.0.9999) pos(2,index)=pos(2,index)-1.0
                 index=index+1
                 iatom(index)=jatom                                           
                 pos(1,index)=pos(1,index-2)+0.5d0         
                 pos(2,index)=pos(2,index-2)        
                 pos(3,index)=pos(3,index-2)+0.5d0
                 if(pos(1,index).ge.0.9999) pos(1,index)=pos(1,index)-1.0
                 if(pos(3,index).ge.0.9999) pos(3,index)=pos(3,index)-1.0
                 index=index+1
                 iatom(index)=jatom                                           
                 pos(1,index)=pos(1,index-3)         
                 pos(2,index)=pos(2,index-3)+0.5d0         
                 pos(3,index)=pos(3,index-3)+0.5d0
                 if(pos(2,index).ge.0.9999) pos(2,index)=pos(2,index)-1.0
                 if(pos(3,index).ge.0.9999) pos(3,index)=pos(3,index)-1.0
               endif
         IF (MULT(JATOM).EQ.0) THEN                                     
            WRITE(6,1020) JATOM,INDEX,MULT(JATOM)                       
            STOP ' NNN: MULT EQ 0'                                      
         ENDIF                     
            DO 55 M=1,MULT(JATOM)-1                                     
               INDEX=INDEX+1                                            
               iatom(index)=jatom                      
!         IF(INDEX.GT.NDIF)then
!         STOP 'NDIF TOO SMALL'                        
!         endif
               READ(24,1011) IATNR1,( POS(J,INDEX),J=1,3)
               if(lattic(1:1).eq.'P'.or.lattic(1:1).eq.'R'.or. &
                  lattic(1:1).eq.'H') then
                 index=index+1
                 iatom(index)=jatom                                           
                 pos(1,index)=pos(1,index-1)+1.0d0         
                 pos(2,index)=pos(2,index-1)         
                 pos(3,index)=pos(3,index-1)    
!                 if(pos(1,index).ge.1.0001) pos(1,index)=pos(1,index)-1.0
                 index=index+1
                 iatom(index)=jatom                                           
                 pos(1,index)=pos(1,index-1)-1.0d0         
                 pos(2,index)=pos(2,index-1)+1.0d0         
                 pos(3,index)=pos(3,index-1)    
!                 if(pos(2,index).ge.1.0001) pos(2,index)=pos(2,index)-1.0
                 index=index+1
                 iatom(index)=jatom                                           
                 pos(1,index)=pos(1,index-1)
                 pos(2,index)=pos(2,index-1)-1.0d0         
                 pos(3,index)=pos(3,index-1)+1.0d0    
!                 if(pos(3,index).ge.1.0001) pos(3,index)=pos(3,index)-1.0
                   if(lattic(1:1).eq.'R') then
                   index=index+1
                   iatom(index)=jatom                                         
                   pos(1,index)=pos(1,index-1)
                   pos(2,index)=pos(2,index-1)-1.0d0         
                   pos(3,index)=pos(3,index-1)-1.0d0    
!                   if(pos(3,index).ge.1.0001) pos(3,index)=pos(3,index)-1.0
                   endif
               endif
               if(lattic(1:3).eq.'CXY') then
                 index=index+1
                 iatom(index)=jatom                                            
                 pos(1,index)=pos(1,index-1)+0.5         
                 pos(2,index)=pos(2,index-1)+0.5         
                 pos(3,index)=pos(3,index-1)    
                 if(pos(1,index).ge.0.9999) pos(1,index)=pos(1,index)-1.0
                 if(pos(2,index).ge.0.9999) pos(2,index)=pos(2,index)-1.0
               endif
               if(lattic(1:3).eq.'CXZ') then
                 index=index+1
                 iatom(index)=jatom                                            
                 pos(1,index)=pos(1,index-1)+0.5         
                 pos(2,index)=pos(2,index-1)             
                 pos(3,index)=pos(3,index-1)+0.5
                 if(pos(1,index).ge.0.9999) pos(1,index)=pos(1,index)-1.0
                 if(pos(3,index).ge.0.9999) pos(3,index)=pos(3,index)-1.0
               endif
               if(lattic(1:3).eq.'CYZ') then
                 index=index+1
                 iatom(index)=jatom                                            
                 pos(1,index)=pos(1,index-1)             
                 pos(2,index)=pos(2,index-1)+0.5         
                 pos(3,index)=pos(3,index-1)+0.5
                 if(pos(2,index).ge.0.9999) pos(2,index)=pos(2,index)-1.0
                 if(pos(3,index).ge.0.9999) pos(3,index)=pos(3,index)-1.0
               endif
               if(lattic(1:1).eq.'B') then
                 index=index+1
                 iatom(index)=jatom                                            
                 pos(1,index)=pos(1,index-1)+0.5         
                 pos(2,index)=pos(2,index-1)+0.5         
                 pos(3,index)=pos(3,index-1)+0.5
                 if(pos(1,index).ge.0.9999) pos(1,index)=pos(1,index)-1.0
                 if(pos(2,index).ge.0.9999) pos(2,index)=pos(2,index)-1.0
                 if(pos(3,index).ge.0.9999) pos(3,index)=pos(3,index)-1.0
               endif
               if(lattic(1:1).eq.'F') then
                 index=index+1
                 iatom(index)=jatom                                            
                 pos(1,index)=pos(1,index-1)+0.5         
                 pos(2,index)=pos(2,index-1)+0.5         
                 pos(3,index)=pos(3,index-1)         
                 if(pos(1,index).ge.0.9999) pos(1,index)=pos(1,index)-1.0
                 if(pos(2,index).ge.0.9999) pos(2,index)=pos(2,index)-1.0
                 index=index+1
                 iatom(index)=jatom                                            
                 pos(1,index)=pos(1,index-2)+0.5         
                 pos(2,index)=pos(2,index-2)         
                 pos(3,index)=pos(3,index-2)+0.5
                 if(pos(1,index).ge.0.9999) pos(1,index)=pos(1,index)-1.0
                 if(pos(3,index).ge.0.9999) pos(3,index)=pos(3,index)-1.0
                 index=index+1
                 iatom(index)=jatom                                            
                 pos(1,index)=pos(1,index-3)         
                 pos(2,index)=pos(2,index-3)+0.5         
                 pos(3,index)=pos(3,index-3)+0.5
                 if(pos(2,index).ge.0.9999) pos(2,index)=pos(2,index)-1.0
                 if(pos(3,index).ge.0.9999) pos(3,index)=pos(3,index)-1.0
               endif
  55        CONTINUE
         index2(jatom)=index
         mult(jatom)=index2(jatom)-index1(jatom)+1
         do 56 i=index1(jatom),index2(jatom)
            pos(1,i)=pos(1,i)*a(1)
            pos(2,i)=pos(2,i)*a(2)
            pos(3,i)=pos(3,i)*a(3)
               if((lattic(1:1).eq.'P'.or.lattic(1:3).eq.'CXZ').and. &
               alpha(3).ne.90.d0) then
!.....monoclinic lattices with gamma.ne.90
            pos(2,i)=pos(2,i)+pos(1,i)*cos(alpha(3)/180.d0*pi)
            pos(1,i)=pos(1,i)*sin(alpha(3)/180.d0*pi)
            if(pos(2,i).lt.-0.0001d0) pos(2,i)=pos(2,i)+a(2)
               endif
! cad......p monoclinic with alpha and beta
               if(lattic(1:1).eq.'P'.and.alpha(2).ne.90.d0) then
            pos(1,i)=pos(1,i)+pos(3,i)*cos(alpha(2)/180.d0*pi)
            pos(3,i)=pos(3,i)*sin(alpha(2)/180.d0*pi)
            if(pos(1,i).lt.-0.0001d0) pos(1,i)=pos(1,i)+a(1)
               endif
               if(lattic(1:1).eq.'P'.and.alpha(1).ne.90.d0) then
            pos(3,i)=pos(3,i)+pos(2,i)*cos(alpha(1)/180.d0*pi)
            pos(2,i)=pos(2,i)*sin(alpha(1)/180.d0*pi)
            if(pos(3,i).lt.-0.0001d0) pos(3,i)=pos(3,i)+a(3)
               endif
! cad
               if(lattic(1:1).eq.'R') then
                 TAU1=pos(1,i)
                 TAU2=pos(2,i)
                 TAU3=pos(3,i)
                 pos(1,i)=Br2(1,1)*TAU1+Br2(1,2)*TAU2+Br2(1,3)*TAU3
                 pos(2,i)=Br2(2,1)*TAU1+Br2(2,2)*TAU2+Br2(2,3)*TAU3
                 pos(3,i)=Br2(3,1)*TAU1+Br2(3,2)*TAU2+Br2(3,3)*TAU3
 1001        continue
             if(pos(1,i).lt.-0.0001d0) then
                pos(1,i)=pos(1,i)+1.0
                goto 1001
             endif
 1002        continue
             if(pos(2,i).lt.-0.0001d0) then
                pos(2,i)=pos(2,i)+1.0
                goto 1002
             endif
 1003        continue
             if(pos(3,i).lt.-0.0001d0) then
                pos(3,i)=pos(3,i)+1.0
                goto 1003
             endif
               endif
 56         continue
         READ(24,1050) NAME(JATOM) 
         read(24,*)           
         read(24,*)           
         read(24,*)           
         write(6,99) name(jatom),mult(jatom),index1(jatom),index2(jatom)      
  99  format(a10,': ',i3,' Atome, Index ',i3,' bis ',i3) 
 50   CONTINUE 
         return
         end
