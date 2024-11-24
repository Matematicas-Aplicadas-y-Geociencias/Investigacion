SUBROUTINE corrector_presion(corr_preso,xo,yo,zo,d_xuo,d_yvo,d_zwo,u_o,v_o,w_o,temp_o,b_o,au_o,av_o,aw_o,conv_po)
USE constantes
! USE mkl95_LAPACK
IMPLICIT NONE
INTEGER :: i,j,k,info
! ********************************
! Variables del flujo e inc'ognita
REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk+1), INTENT(inout) :: corr_preso
REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk+1), INTENT(out)   :: b_o
REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk+1), INTENT(in)    :: temp_o
REAL(kind=DBL), DIMENSION(mi,nj+1,lk+1),   INTENT(in)    :: u_o,au_o
REAL(kind=DBL), DIMENSION(mi+1,nj,lk+1),   INTENT(in)    :: v_o,av_o
REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk),   INTENT(in)    :: w_o,aw_o
REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk+1)                :: fcorr_pres
! ***************************
! Variables de la tridiagonal
REAL(kind=DBL), DIMENSION(mi+1) :: APi,Rxi
REAL(kind=DBL), DIMENSION(mi)   :: AW,AE
REAL(kind=DBL), DIMENSION(nj+1) :: APj,Ryj
REAL(kind=DBL), DIMENSION(nj)   :: AS,AN
REAL(kind=DBL), DIMENSION(lk+1) :: APk,Rzk
REAL(kind=DBL), DIMENSION(lk)   :: AB,AT
! ***************************
! Variables de interpolaci'on
REAL(kind=DBL) :: uw,ue,vs,vn,wb,wt,dw,de,ds,dn,db,dtop
REAL(kind=DBL) :: delta_x,delta_y,delta_z
REAL(kind=DBL) :: a_w,a_e,a_s,a_n,a_b,a_top,conv_po,rel_po
!*******************************************************
!Variables de la malla, volumen de control e incrementos
REAL(kind=DBL), DIMENSION(mi+1), INTENT(in) :: xo
REAL(kind=DBL), DIMENSION(nj+1), INTENT(in) :: yo
REAL(kind=DBL), DIMENSION(lk+1), INTENT(in) :: zo
REAL(kind=DBL), DIMENSION(mi-1), INTENT(in) :: d_xuo
REAL(kind=DBL), DIMENSION(nj-1), INTENT(in) :: d_yvo
REAL(kind=DBL), DIMENSION(lk-1), INTENT(in) :: d_zwo
! ********************
! Variables auxiliares
REAL(kind=DBL) alpha,beta,kappa,sigma
! ***********************************
! relajaci'on para la presi'on
rel_po = 0.85_DBL
! ***********************************
! auxiliar para calcular dcorr_preso
! 2 continue
fcorr_pres = corr_preso
! *************************
! Comienza regi'on paralela
! !$OMP PARALLEL
! !$OMP SECTIONS
! ****************************
! Secci'on TDMA en direccion x
! !$OMP SECTION
!$OMP  PARALLEL DO DEFAULT(NONE)&
!$OMP& PRIVATE(AW,APi,AE,Rxi,uw,ue,vs,vn,wb,wt,a_w,a_e,a_s,a_n,a_b,a_top,delta_x,delta_y,delta_z,alpha,beta,kappa,sigma)&
!$OMP& SHARED(u_o,v_o,w_o,b_o,au_o,av_o,aw_o,corr_preso,d_xuo,d_yvo,d_zwo,rel_po)
DO k = 2, lk
  DO j = 2, nj
    !***********************
    !Condiciones de frontera
    APi(1) = 1._DBL
    AE(1)  = cero
    Rxi(1) = cero
    !************
    APi(mi+1) = 1._DBL
    AW(mi)    = cero
    Rxi(mi+1) = cero
    !*************************
    !Llenado de la matriz en x
    DO i = 2, mi
      !**************************
      !Interpolaciones necesarias
      ue  = u_o(i,j,k)
      uw  = u_o(i-1,j,k)
      vn  = v_o(i,j,k)
      vs  = v_o(i,j-1,k)
      wt  = w_o(i,j,k)
      wb  = w_o(i,j,k-1)
      a_e = au_o(i,j,k)
      a_w = au_o(i-1,j,k)
      a_n = av_o(i,j,k)
      a_s = av_o(i,j-1,k)
      a_top = aw_o(i,j,k)
      a_b = aw_o(i,j,k-1)
      delta_x = d_xuo(i-1)
      delta_y = d_yvo(j-1)
      delta_z = d_zwo(k-1)
      !************************************************
      AE(i)     =-(delta_y*delta_y*delta_z*delta_z/a_e)
      IF (i == mi) AE(i) = cero
      AW(i-1)   =-(delta_y*delta_y*delta_z*delta_z/a_w)
      IF (i == 2) AW(i-1) = cero
      alpha     = delta_x*delta_x*delta_z*delta_z/a_n
      IF (j ==nj) alpha = cero
      beta      = delta_x*delta_x*delta_z*delta_z/a_s
      IF (j == 2) beta = cero
      kappa   = delta_x*delta_x*delta_y*delta_y/a_top
      IF (k == lk) kappa = cero
      sigma   = delta_x*delta_x*delta_y*delta_y/a_b
      IF (k == 2) sigma = cero
      APi(i)    =(-AW(i-1)-AE(i)+alpha+beta+kappa+sigma) / rel_po
      b_o(i,j,k)= (uw-ue)*delta_y*delta_z+(vs-vn)*delta_x*delta_z+(wb-wt)*delta_x*delta_y
      Rxi(i)    = alpha*corr_preso(i,j+1,k)+beta*corr_preso(i,j-1,k)+kappa*corr_preso(i,j,k+1)+sigma*corr_preso(i,j,k-1)+b_o(i,j,k) + &
                  APi(i)*(1._DBL-rel_po)*corr_preso(i,j,k)
      !************************************************
    END DO
    !**************************
    CALL tri(AW,APi,AE,Rxi,mi+1)
    DO i = 1, mi+1
      corr_preso(i,j,k) = Rxi(i)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO
! *******************************
! Secci'on de TDMA en direccion y
! !$OMP SECTION
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(AS,APj,AN,Ryj,uw,ue,vs,vn,wb,wt,a_w,a_e,a_s,a_n,a_b,a_top,delta_x,delta_y,delta_z,alpha,beta,kappa,sigma) SHARED(u_o,v_o,w_o,b_o,au_o,av_o,aw_o,corr_preso,d_xuo,d_yvo,d_zwo,rel_po)
DO k = 2, lk
  DO i = 2, mi
    !***********************
    !Condiciones de frontera
    APj(1) = 1._DBL
    AN(1)  = cero
    Ryj(1) = cero
    !*****************
    APj(nj+1) = 1._DBL
    AS(nj)    = cero
    Ryj(nj+1) = cero
    !*************************
    !Llenado de la matriz en x
    DO j = 2, nj
      !**************************
      !Interpolaciones necesarias
      ue  = u_o(i,j,k)
      uw  = u_o(i-1,j,k)
      vn  = v_o(i,j,k)
      vs  = v_o(i,j-1,k)
      wt  = w_o(i,j,k)
      wb  = w_o(i,j,k-1)
      a_e = au_o(i,j,k)
      a_w = au_o(i-1,j,k)
      a_n = av_o(i,j,k)
      a_s = av_o(i,j-1,k)
      a_top = aw_o(i,j,k)
      a_b = aw_o(i,j,k-1)
      delta_x = d_xuo(i-1)
      delta_y = d_yvo(j-1)
      delta_z = d_zwo(k-1)
      !**********************************************
      alpha = delta_y*delta_y*delta_z*delta_z/a_e
      IF (i == mi) alpha = cero
      beta = delta_y*delta_y*delta_z*delta_z/a_w
      IF (i == 2) beta = cero
      AN(j) =-(delta_x*delta_x*delta_z*delta_z/a_n)
      IF (j ==nj) AN(j) = cero
      AS(j-1) =-(delta_x*delta_x*delta_z*delta_z/a_s)
      IF (j == 2) AS(j-1) = cero
      kappa   = delta_x*delta_x*delta_y*delta_y/a_top
      IF (k == lk) kappa = cero
      sigma   = delta_x*delta_x*delta_y*delta_y/a_b
      IF (k == 2) sigma = cero
      APj(j)    = (alpha+beta-AN(j)-AS(j-1)+kappa+sigma) / rel_po
      b_o(i,j,k)= (uw-ue)*delta_y*delta_z+(vs-vn)*delta_x*delta_z+(wb-wt)*delta_x*delta_y
      Ryj(j)    = alpha*corr_preso(i+1,j,k)+beta*corr_preso(i-1,j,k)+kappa*corr_preso(i,j,k+1)+sigma*corr_preso(i,j,k-1)+b_o(i,j,k)+&
                  APj(j)*(1._DBL-rel_po)*corr_preso(i,j,k)
      !**********************************************
    END DO
    !**************************
    CALL tri(AS,APj,AN,Ryj,nj+1)
    DO j = 1, nj+1
      corr_preso(i,j,k) = Ryj(j)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO
! *******************************
! Secci'on de TDMA en direccion z
! !$OMP SECTION
!$OMP  PARALLEL DO DEFAULT(NONE)&
!$OMP& PRIVATE(AB,APk,AT,Rzk,uw,ue,vs,vn,wb,wt,a_w,a_e,a_s,a_n,a_b,a_top,delta_x,delta_y,delta_z,alpha,beta,kappa,sigma)&
!$OMP& SHARED(u_o,v_o,w_o,b_o,au_o,av_o,aw_o,corr_preso,d_xuo,d_yvo,d_zwo,rel_po)
DO j = 2, nj
  DO i = 2, mi
    !***********************
    !Condiciones de frontera
    APk(1) = 1._DBL
    AT(1)  = cero
    Rzk(1) = cero
    !*****************
    APk(lk+1) = 1._DBL
    AB(lk)    = cero
    Rzk(lk+1) = cero
    !*************************
    !Llenado de la matriz en x
    DO k = 2, lk
      !**************************
      !Interpolaciones necesarias
      ue  = u_o(i,j,k)
      uw  = u_o(i-1,j,k)
      vn  = v_o(i,j,k)
      vs  = v_o(i,j-1,k)
      wt  = w_o(i,j,k)
      wb  = w_o(i,j,k-1)
      a_e = au_o(i,j,k)
      a_w = au_o(i-1,j,k)
      a_n = av_o(i,j,k)
      a_s = av_o(i,j-1,k)
      a_top = aw_o(i,j,k)
      a_b = aw_o(i,j,k-1)
      delta_x = d_xuo(i-1)
      delta_y = d_yvo(j-1)
      delta_z = d_zwo(k-1)
      !**********************************************
      alpha     = delta_y*delta_y*delta_z*delta_z/a_e
      IF (i == mi) alpha = cero
      beta      = delta_y*delta_y*delta_z*delta_z/a_w
      IF (i == 2) beta = cero
      kappa     = delta_x*delta_x*delta_z*delta_z/a_n
      IF (j ==nj) kappa = cero
      sigma     = delta_x*delta_x*delta_z*delta_z/a_s
      IF (j == 2) sigma = cero
      AT(k)     =-(delta_x*delta_x*delta_y*delta_y/a_top)
      IF (k == lk) AT(k) = cero
      AB(k-1)   =-(delta_x*delta_x*delta_y*delta_y/a_b)
      IF (k == 2) AB(k-1) = cero
      APk(k)    = (alpha+beta+kappa+sigma-AT(k)-AB(k-1)) / rel_po
      b_o(i,j,k)= (uw-ue)*delta_y*delta_z+(vs-vn)*delta_x*delta_z+(wb-wt)*delta_x*delta_y
      Rzk(k)    = alpha*corr_preso(i+1,j,k)+beta*corr_preso(i-1,j,k)+kappa*corr_preso(i,j+1,k)+sigma*corr_preso(i,j-1,k)+b_o(i,j,k)+&
                  APk(k)*(1._DBL-rel_po)*corr_preso(i,j,k)
      ! **********************************************
    END DO
    ! **************************
    CALL tri(AB,APk,AT,Rzk,lk+1)
    DO k = 1, lk+1
      corr_preso(i,j,k) = Rzk(k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO
! ! *******************
! ! Finalizan secciones
! ! !$OMP END SECTIONS
! ! ****************************
! ! Finaliza la regi'on paralela
! ! !$OMP END PARALLEL
! IF(MAXVAL(DABS(fcorr_pres-corr_preso)) < conv_po)THEN
! !   WRITE(*,*) MAXVAL(DABS(corr_preso-corr_presx)), MAXVAL(DABS(corr_preso-corr_presy)), conv_po
!   GOTO 3
! ELSE
!   WRITE(*,*) MAXVAL(DABS(corr_preso-fcorr_pres))
! !   itera = itera + 1
! !   WRITE(*,*) maxval(dabs((teta_xo-teta_yo)))
!   GOTO 2
! ENDIF
! 3 CONTINUE
END SUBROUTINE
