SUBROUTINE efgsplit(jatom,lqmax,lefg,mefg,nefg,rhoefg,r)
  USE defs; USE struk
  USE param
  IMPLICIT NONE
  REAL(8)     :: R(NRAD),rhoefg(NRAD,5,23)
  INTEGER     :: jatom,lqmax,lefg(3,5,23),mefg(3,5,23),nefg(5)
  REAL(8)     :: vcol(-2:2,0:2),efg,efgfact,sqrt2
  REAL(8)     :: qmat(3,3),eigenvalue(3,3),eigenvector(3,3),work(9)
  REAL(8)     :: ev_tot(3,3)
  REAL(8)     :: efg_comb(0:3,0:3,-2:2,0:2)
  INTEGER     :: lq,ll,mm,jefg,l,lp,m,mp,jr,i,j,k,info
  SQRT2=SQRT(2.0D0) 
  efgfact=-2.d0*9.71736408
  vcol(-2:2,0:2)=ZERO
  efg_comb=ZERO
  DO LQ=1,LQMAX                                              
     LL=LEFG(1,LQ,1)                                                
     MM=MEFG(1,LQ,1)                                                
!     WRITE(6,3000) LL,MM                                            
     WRITE(6,3004)                                                  
     DO JEFG=1,NEFG(LQ)                                         
        LL=LEFG(1,LQ,JEFG)
        L =LEFG(2,LQ,JEFG)
        LP=LEFG(3,LQ,JEFG)
        MM=MEFG(1,LQ,JEFG)
        M =MEFG(2,LQ,JEFG)
        MP=MEFG(3,LQ,JEFG)
        !         QPART=?       
        DO JR=1,JRI(JATOM)
           RHOEFG(JR,LQ,JEFG)=RHOEFG(JR,LQ,JEFG)/R(JR)**3
        enddo
        CALL CHARGE(R,1.D0,1.D0,RHOEFG(1,LQ,JEFG),DX(JATOM), &
             JRI(JATOM),EFG) 
        if(mm.eq.0) then
        WRITE(6,3005) L,LL,LP,M,MM,MP,EFG*SQRT(15.D0/16.D0/PI)*4.d0*PI/5.d0*efgfact
        else
        WRITE(6,3005) L,LL,LP,M,MM,MP,EFG*SQRT(15.D0/16.D0/PI)*4.d0*PI/5.d0*efgfact*sqrt2
        endif
        vcol(ll,mm)=vcol(ll,mm)+efg
        IF(lp>l) THEN
           efg_comb(l,lp,ll,mm)=efg_comb(l,lp,ll,mm)+efg
        ELSE
           efg_comb(lp,l,ll,mm)=efg_comb(lp,l,ll,mm)+efg
        ENDIF
     ENDDO
  ENDDO

  WRITE(6,*) '-Total valence contribution to EFG tensor -'
  vcol=vcol*SQRT(15.D0/16.D0/PI)*4.d0*PI/5.d0*efgfact
  vcol(2,2)=vcol(2,2)*sqrt2
  vcol(2,1)=vcol(2,1)*sqrt2
  vcol(-2,1)=vcol(-2,1)*sqrt2
  vcol(-2,2)=vcol(-2,2)*sqrt2
  WRITE(6,3010) vcol(2,0),vcol(2,2),vcol(-2,2),vcol(2,1),vcol(-2,1)
  qmat(1,1)=Vcol(2,2)-Vcol(2,0)/SQRT(3.d0)
  qmat(2,2)=-Vcol(2,2)-Vcol(2,0)/SQRT(3.d0)
  qmat(3,3)=2.d0*Vcol(2,0)/SQRT(3.d0) 
  qmat(1,2)=vcol(-2,2)
  qmat(2,1)=vcol(-2,2)
  qmat(1,3)=vcol(2,1)
  qmat(3,1)=vcol(2,1)
  qmat(2,3)=vcol(-2,1)
  qmat(3,2)=vcol(-2,1)
  eigenvector=qmat
  CALL DSYEV('V','U',3,eigenvector,3,eigenvalue,work,9,info)
  IF(info>0) WRITE(6,*) 'Diagonalization failed'
  WRITE(6,1350) ((qmat(i,j),j=1,3),i=1,3)
  WRITE(6,1355) (eigenvalue(i,1),(eigenvector(i,j),j=1,3),i=1,3)

  ev_tot=eigenvector
  ! ___________ Components ____________________________________
  efg_comb=efg_comb*SQRT(15.D0/16.D0/PI)*4.d0*PI/5.d0*efgfact

  DO jefg=1,5
     IF(jefg==1) THEN 
        l=1; lp=1
        WRITE(6,*) '------ pp ------------------'
     ELSEIF(jefg==2) THEN
        l=0;lp=2
        WRITE(6,*) '------ sd ------------------'
     ELSEIF(jefg==3) THEN
        l=2;lp=2
        WRITE(6,*) '------ dd ------------------'
     ELSEIF(jefg==4) THEN
        l=1;lp=3
        WRITE(6,*) '------ pf ------------------'
     ELSE
        l=3;lp=3
        WRITE(6,*) '------ ff ------------------'
     ENDIF
     Efg_comb(l,lp,2,2)=Efg_comb(l,lp,2,2)*sqrt2
     Efg_comb(l,lp,2,1)=Efg_comb(l,lp,2,1)*sqrt2
     Efg_comb(l,lp,-2,1)=Efg_comb(l,lp,-2,1)*sqrt2
     Efg_comb(l,lp,-2,2)=Efg_comb(l,lp,-2,2)*sqrt2
     qmat(1,1)=Efg_comb(l,lp,2,2)-Efg_comb(l,lp,2,0)/SQRT(3.d0)
     qmat(2,2)=-Efg_comb(l,lp,2,2)-Efg_comb(l,lp,2,0)/SQRT(3.d0)
     qmat(3,3)=2.d0*Efg_comb(l,lp,2,0)/SQRT(3.d0) 
     qmat(1,2)=efg_comb(l,lp,-2,2)
     qmat(2,1)=efg_comb(l,lp,-2,2)
     qmat(1,3)=efg_comb(l,lp,2,1)
     qmat(3,1)=efg_comb(l,lp,2,1)
     qmat(2,3)=efg_comb(l,lp,-2,1)
     qmat(3,2)=efg_comb(l,lp,-2,1)
     eigenvector=qmat
     WRITE(6,1350) ((qmat(i,j),j=1,3),i=1,3)
     CALL DSYEV('V','U',3,eigenvector,3,eigenvalue,work,9,info)

     IF(info>0) WRITE(6,*) 'Diagonalization failed'
     WRITE(6,1355) (eigenvalue(i,1),(eigenvector(i,j),j=1,3),i=1,3)
     WRITE(6,*) 'partial QMAT projection on eigenvectors of total tensor'
     qmat=MATMUL(TRANSPOSE(ev_tot),MATMUL(qmat,ev_tot))
     WRITE(6,1350) ((qmat(i,j),j=1,3),i=1,3)
  ENDDO
106 FORMAT(3(10X,3F10.5,/)) 
1350 FORMAT(3(10X,3F11.5,/))
1355 FORMAT(/,9X,'Eigenvalues        Eigen vectors (rows) ',/,3(12X,F8.4,7X,3F8.4,/))
3000 FORMAT(//,3X,'EFG SPLIT:  LL=',I2,3X,'MM=',I2,                     &
       /,3X,'--------------------------',/)                       
3004 FORMAT(16X,'L',2X,'LL',2X,'LP',5X,'M ',2X,'MM',2X,'MP',4X,         &
       'CONTRIB TO EFG')                                          
3005 FORMAT(15X,I2,2X,I2,2X,I2,5X,I2,2X,I2,2X,I2,5X,F10.3)             
3010 FORMAT(/,15X,'TOTAL VCOUL 20 =',F12.3, &
       /,15X,'TOTAL VCOUL 22 =',F12.3, &
       /,15X,'TOTAL VCOUL 22M=',F12.3, &
       /,15X,'TOTAL VCOUL 21 =',F12.3, &
       /,15X,'TOTAL VCOUL 21M=',F12.3, &
       /,15x,'CARTH. TENSOR: ',&
       /,15x,'V_xx=V_22-V_20/sqrt3; V_yy=-V_22-V_20/sqrt3; V_zz=2*V_20/sqrt3',&
       /,15x,'V_xy=V_22M; V_xz=V_21; V_yz=V_21M')

END SUBROUTINE efgsplit
