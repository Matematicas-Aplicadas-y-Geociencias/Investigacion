SUBROUTINE corrector_presion(corr_preso,d_xuo,d_yvo,u_o,v_o,b_o,au_o,av_o,dcorr_preso)
USE malla
! USE mkl95_LAPACK
IMPLICIT NONE
INTEGER :: i,j,info
!*********************************
! Variables del flujo e inc'ognita
REAL(kind=DBL), DIMENSION(mi+1,nj+1), INTENT(inout) :: corr_preso
REAL(kind=DBL), DIMENSION(mi+1,nj+1), INTENT(out)   :: b_o,dcorr_preso
REAL(kind=DBL), DIMENSION(mi,nj+1),   INTENT(in)    :: u_o,au_o
REAL(kind=DBL), DIMENSION(mi+1,nj),   INTENT(in)    :: v_o,av_o
REAL(kind=DBL), DIMENSION(mi+1,nj+1) :: fcorr_preso
!*********************************
!Variables de la tridiagonal
REAL(kind=DBL), DIMENSION(mi+1) :: ACi,Rxi
REAL(kind=DBL), DIMENSION(mi)   :: AI,AD
REAL(kind=DBL), DIMENSION(nj+1) :: ACj,Ryj
REAL(kind=DBL), DIMENSION(nj)   :: AS,AN
!*********************
!Variables de interpolaci'on
REAL(kind=DBL) :: ui,ud,vn,vs,delta_x,delta_y,a_i,a_d,a_s,a_n
!*********************
!Variables de la malla, volumen de control, incremento de tiempo y nums Peclet, Reynolds y Richardson
REAL(kind=DBL), DIMENSION(mi-1), INTENT(in) :: d_xuo
REAL(kind=DBL), DIMENSION(nj-1), INTENT(in) :: d_yvo
!********************
!Variables auxiliares
REAL(kind=DBL) :: alpha,beta,rel
!**********************************
!auxiliar para calcular dcorr_preso
! DO
fcorr_preso = corr_preso
! corr_preso = 0._DBL
rel        = 0.85_DBL
!*******************
!TDMA en direccion x
!$OMP  PARALLEL DO DEFAULT(NONE)&
!$OMP& PRIVATE(AI,ACi,AD,Rxi,ui,ud,vs,vn,a_i,a_d,a_s,a_n,delta_x,delta_y,alpha,beta,info)&
!$OMP& SHARED(u_o,v_o,b_o,au_o,av_o,corr_preso,d_xuo,d_yvo,rel,fcorr_preso)
DO j = 2, nj
  !***********************
  !Condiciones de frontera
  ACi(1) = 1._DBL
  AD(1)  = cero
  Rxi(1) = cero
  !************
  ACi(mi+1) = 1._DBL
  AI(mi)    = cero
  Rxi(mi+1) = cero
  !*************************
  !Llenado de la matriz en x
  DO i = 2, mi
    !**************************
    !Interpolaciones necesarias
    ud  = u_o(i,j)
    ui  = u_o(i-1,j)
    vs  = v_o(i,j-1)
    vn  = v_o(i,j)
    a_d = au_o(i,j)
    a_i = au_o(i-1,j)
    a_s = av_o(i,j-1)
    a_n = av_o(i,j)
    delta_x = d_xuo(i-1)
    delta_y = d_yvo(j-1)
    !******************************
    AI(i-1) =-(delta_y*delta_y/a_i)
    AD(i)   =-(delta_y*delta_y/a_d)
    alpha   = delta_x*delta_x/a_s
    beta    = delta_x*delta_x/a_n
    ACi(i)  = (-AI(i-1)-AD(i)+alpha+beta) / rel
    b_o(i,j)= (ui-ud)*delta_y+(vs-vn)*delta_x
    Rxi(i)  = alpha*corr_preso(i,j-1)+beta*corr_preso(i,j+1)+b_o(i,j)+(1._DBL-rel)*ACi(i)*corr_preso(i,j)
    !**************************
  END DO
  !**************************
!   CALL gtsv(AI,ACi,AD,Rxi,info)
  CALL tri(AI,ACi,AD,Rxi,mi+1)
  DO i = 1, mi+1
    corr_preso(i,j) = Rxi(i)
  END DO
END DO
!$OMP END PARALLEL DO
!*******************
!TDMA en direccion y
!$OMP  PARALLEL DO DEFAULT(NONE)&
!$OMP& PRIVATE(AS,ACj,AN,Ryj,ui,ud,vs,vn,a_i,a_d,a_s,a_n,delta_x,delta_y,alpha,beta,info)&
!$OMP& SHARED(u_o,v_o,b_o,au_o,av_o,corr_preso,d_xuo,d_yvo,rel,fcorr_preso)
DO i = 2, mi
  !***********************
  !Condiciones de frontera
  ACj(1) = 1._DBL
  AN(1)  = cero
  Ryj(1) = cero
  !**********
  ACj(nj+1) = 1._DBL
  AS(nj)    = cero
  Ryj(nj+1) = cero
  !*************************
  !Llenado de la matriz en x
  DO j = 2, nj
    !**************************
    !Interpolaciones necesarias
    ud  = u_o(i,j)
    ui  = u_o(i-1,j)
    vs  = v_o(i,j-1)
    vn  = v_o(i,j)
    a_d = au_o(i,j)
    a_i = au_o(i-1,j)
    a_s = av_o(i,j-1)
    a_n = av_o(i,j)
    delta_x = d_xuo(i-1)
    delta_y = d_yvo(j-1)
    !************************
    alpha   = delta_y*delta_y/a_i
    beta    = delta_y*delta_y/a_d
    AS(j-1) =-(delta_x*delta_x/a_s)
    AN(j)   =-(delta_x*delta_x/a_n)
    ACj(j)  = (alpha+beta-AS(j-1)-AN(j)) / rel
    b_o(i,j)= (ui-ud)*delta_y+(vs-vn)*delta_x
    Ryj(j)  = alpha*corr_preso(i-1,j)+beta*corr_preso(i+1,j)+b_o(i,j)+(1._DBL-rel)*ACj(j)*corr_preso(i,j)
  END DO
  !***************************
!   CALL gtsv(AS,ACj,AN,Ryj,info)
  CALL tri(AS,ACj,AN,Ryj,nj+1)
  DO j = 1, nj+1
    corr_preso(i,j) = Ryj(j)
  END DO
END DO
!$OMP END PARALLEL DO
WHERE(corr_preso /= cero)
  dcorr_preso = (corr_preso-fcorr_preso)!/corr_preso
ELSEWHERE
  dcorr_preso = corr_preso-fcorr_preso
END WHERE
! IF(MAXVAL(DABS(dcorr_preso))<1.e-3_DBL)EXIT
! WRITE(*,*) 'corrector presion ', MAXVAL(DABS(dcorr_preso)), MAXVAL(DABS(b_o))
! END DO
END SUBROUTINE
