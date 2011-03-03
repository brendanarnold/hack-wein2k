      subroutine readcs(readcrit,readsurf)
      use crit
      implicit none
      include 'param.inc'

      real*8 cth,th,ph,wcth,wph,rs,lap,rho,x,y,z
      real*8 themin,themax,phimin,phimax,h,frmin,r0,dr0,weight
      real*8 br1,br2,br3,br4

      integer index,nsimax,nrs,nth,nph
      integer iat,ipos,i,j

      logical readcrit,readsurf
      
      character*4 ctemp

      COMMON /BRAV/ br1(3,3),br2(3,3),br3(3,3),br4(3,3)
!      COMMON /CRIT/ pc(3,NNPOS*NDIF),evpc(3,NNPOS*NDIF), &
!         zpc(3,3,NNPOS*NDIF),pcrb(3,NNPOS*NDIF),npc,icpc(NNPOS*NDIF), &
!         npcm3
      COMMON /INTEG/ cth(NGAUSS),th(NGAUSS),ph(NGAUSS),wcth(NGAUSS), &
         wph(NGAUSS),rs(NGAUSS,NGAUSS),nth,nph
      COMMON /SRF/ themin,themax,phimin,phimax,h,frmin,r0,dr0, &
         index,nsimax
      
      rewind(22)
      npc=0
      npcm3=0
      write(6,*) 'READCRIT'
      
 30   npc=npc+1
      read (22,*,END=20) ctemp,iat,ipos,pcrb(1,npc),pcrb(2,npc), &
         pcrb(3,npc),evpc(1,npc),evpc(2,npc),evpc(3,npc),icpc(npc), &
         lap,rho
      write(6,202) iat,ipos,pcrb(1,npc),pcrb(2,npc),pcrb(3,npc), &
         evpc(1,npc),evpc(2,npc),evpc(3,npc),icpc(npc),lap,rho
      do i=1,3
        pc(i,npc)=pcrb(1,npc)*br1(1,i)+pcrb(2,npc)*br1(2,i)+ &
           pcrb(3,npc)*br1(3,i)
      enddo
 202  format(':PC ',2I5,3F12.8,3E12.4,I4,2E12.4)
      if(icpc(npc).eq.-3) npcm3=npcm3+1
      goto 30
      
 20   npc=npc-1
      if(npc.gt.0) readcrit=.true.
      write(6,*) 'Number of PC = ',npc
      
      rewind(21)
      write(6,*) 'READSURF'
      read(21,*) index,x,y,z
      read(21,*) nth,themin,themax
      read(21,*) nph,phimin,phimax
      write(6,*) index,x,y,z
      write(6,*) nth,themin,themax
      write(6,*) nph,phimin,phimax
      do i=1,nth
        do j=1,nph
          rs(i,j)=0.d0
        enddo
      enddo
      nrs=0
      do i=1,nth
        do j=1,nph
          nrs=nrs+1
          read(21,*,END=21) th(i),ph(j),rs(i,j),weight
          write(6,200) i,j,th(i),ph(j),rs(i,j),weight
 200      format(2I4,2F12.8,2E16.8)
        enddo
      enddo
 21   continue
      nrs=nrs-1
      if(nrs.gt.0) readsurf=.true.
      write(6,*) 'Number of rs read = ',nrs

      return
      end

