      subroutine rewr(l1,l2)
!-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      CHARACTER*4       LATTIC,rela                                          
      CHARACTER*79      TITLE,  name                                        
      dimension pos(3)
      dimension alpha(3),a(3),rotold(3,3)
!-----------------------------------------------------------------------
!.....START READING FILE STRUCT                                           
 17   FORMAT(6F10.6)                                                    
 1011 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
 1012 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,7X,I2,8X,I2)     
 1510 FORMAT(A79)                                                       
 1511 FORMAT(A4,23X,I3,/,13X,A4)                                        
      READ(l1,1510) TITLE 
      READ(l1,1511) LATTIC,NAT,RELA                                          
      READ(l1,17) A(1),A(2),A(3),alpha(1),alpha(2),alpha(3)   
      write(l2,1510) TITLE 
      write(l2,1511) LATTIC,NAT,RELA                                          
      write(l2,17) A(1),A(2),A(3),alpha(1),alpha(2),alpha(3)   
      DO 50 JATOM = 1,NAT                                               
         READ(l1,1012) IATNR,( POS(J),J=1,3 ),MULT  
         write(l2,1012) IATNR,( POS(J),J=1,3 ),MULT  
            DO 55 M=1,MULT-1                                     
               read (l1,1011) IATNR1,( POS(J),J=1,3)
               write(l2,1011) IATNR1,( POS(J),J=1,3)
  55        CONTINUE
!        read/write radial info and rotmat  
         do i=1,4
         read (l1,1510) NAME 
         write(l2,1510) NAME 
         enddo
!         read(l1,1012) ((rotold(i,j),j=1,3),i=1,3)          
!         write(l2,1012) ((rotold(i,j),j=1,3),i=1,3)          
 50   CONTINUE 
         return
         end
