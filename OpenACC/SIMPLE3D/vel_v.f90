SUBROUTINE vel_v(xo,yvo,zo,fex,fez,d_yvo,d2_yvo,d_xuo,d_zwo,u_o,v_o,v_anto,w_o,pres_o,gamma_vo,dt_o,dv_o,av_o,rel_vo)
USE constantes
! USE mkl95_LAPACK
IMPLICIT NONE
INTEGER :: i,j,k,info
!**********************************
! Variables del flujo e inc'ognitas
REAL(kind=DBL), DIMENSION(mi+1,nj,lk+1), INTENT(inout) :: v_o
REAL(kind=DBL), DIMENSION(mi+1,nj,lk+1), INTENT(out)   :: av_o,dv_o
REAL(kind=DBL), DIMENSION(mi+1,nj,lk+1), INTENT(in)    :: v_anto,gamma_vo
REAL(kind=DBL), DIMENSION(mi,nj+1,lk+1), INTENT(in)    :: u_o
REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk), INTENT(in)    :: w_o
REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk+1), INTENT(in)  :: pres_o
REAL(kind=DBL), DIMENSION(mi+1,nj,lk+1)                :: fv,v_1
!***************************
!Variables de la tridiagonal
REAL(kind=DBL), DIMENSION(mi+1) :: APi,Rxi
REAL(kind=DBL), DIMENSION(mi)   :: AW,AE
REAL(kind=DBL), DIMENSION(nj)   :: APj,Ryj
REAL(kind=DBL), DIMENSION(nj-1) :: AS,AN
REAL(kind=DBL), DIMENSION(lk+1) :: APk,Rzk
REAL(kind=DBL), DIMENSION(lk)   :: AB,AT
!***************************
!Variables de interpolaci'on
REAL(kind=DBL) :: uw,ue,vs,vn,wb,wt,dw,de,ds,dn,db,dtop,gamma_w,gamma_e,gamma_s,gamma_n,gamma_b,gamma_top
REAL(kind=DBL) :: delta_x,delta_y,delta_z
!*************************************************************
!Variables de la malla,volumen de control,incremento de tiempo
REAL(kind=DBL), DIMENSION(mi+1), INTENT(in) :: xo,fex
REAL(kind=DBL), DIMENSION(nj),   INTENT(in) :: yvo
REAL(kind=DBL), DIMENSION(lk+1), INTENT(in) :: zo,fez
REAL(kind=DBL), DIMENSION(mi-1), INTENT(in) :: d_xuo
REAL(kind=DBL), DIMENSION(nj-1), INTENT(in) :: d_yvo,d2_yvo
REAL(kind=DBL), DIMENSION(lk-1), INTENT(in) :: d_zwo
REAL(kind=DBL), INTENT(in)                  :: dt_o,rel_vo
!********************
!Variables auxiliares
REAL(kind=DBL) :: alpha,beta,kappa,sigma
!*************************
!auxiliar para calcular dv
fv = v_o
!*******************
!TDMA en direccion x
!$OMP  PARALLEL DO DEFAULT(NONE)&
!$OMP& PRIVATE(AW,APi,AE,Rxi,uw,ue,vs,vn,wb,wt,dw,de,ds,dn,db,dtop,gamma_w,gamma_e,gamma_s,gamma_n,gamma_b,gamma_top,delta_x,delta_y,delta_z,alpha,beta,kappa,sigma)&
!$OMP& SHARED(v_o,fv,v_anto,pres_o,gamma_vo,av_o,u_o,w_o,d_yvo,d2_yvo,xo,zo,fex,fez,d_xuo,d_zwo,rel_vo,dt_o)
DO k = 2, lk
  DO j = 2, nj-1
    !*************************
    !Llenado de la matriz en x
    DO i = 2, mi
      !**************************
      !Interpolaciones necesarias
      ue = u_o(i,j,k)  +d_yvo(j-1)/d2_yvo(j)*(u_o(i,j+1,k)-u_o(i,j,k))
      uw = u_o(i-1,j,k)+d_yvo(j-1)/d2_yvo(j)*(u_o(i-1,j+1,k)-u_o(i-1,j,k))
      vn = (v_o(i,j+1,k)+v_o(i,j,k))/2._DBL
      vs = (v_o(i,j,k)+v_o(i,j-1,k))/2._DBL
      wt = w_o(i,j,k)  +d_yvo(j-1)/d2_yvo(j)*(w_o(i,j+1,k)-w_o(i,j,k))
      wb = w_o(i,j,k-1)+d_yvo(j-1)/d2_yvo(j)*(w_o(i,j+1,k-1)-w_o(i,j,k-1))
      de = xo(i+1)-xo(i)
      dw = xo(i)-xo(i-1)
      dn = d_yvo(j)
      ds = d_yvo(j-1)
      dtop = zo(k+1)-zo(k)
      db = zo(k)-zo(k-1)
      gamma_e = 1._DBL/((1._DBL-fex(i+1)) / gamma_vo(i,j,k) + fex(i+1) / gamma_vo(i+1,j,k))
      gamma_w = 1._DBL/((1._DBL-fex(i)) / gamma_vo(i-1,j,k) + fex(i) / gamma_vo(i,j,k))
      gamma_n = 2._DBL * gamma_vo(i,j+1,k) * gamma_vo(i,j,k) / (gamma_vo(i,j+1,k) + gamma_vo(i,j,k))
      gamma_s = 2._DBL * gamma_vo(i,j-1,k) * gamma_vo(i,j,k) / (gamma_vo(i,j-1,k) + gamma_vo(i,j,k))
      gamma_top = 1._DBL/((1._DBL-fez(k+1)) / gamma_vo(i,j,k) + fez(k+1) / gamma_vo(i,j,k+1))
      gamma_b = 1._DBL/((1._DBL-fez(k)) / gamma_vo(i,j,k-1) + fez(k) / gamma_vo(i,j,k))
      delta_x = d_xuo(i-1)
      delta_y = d2_yvo(j)/2._DBL
      delta_z = d_zwo(k-1)
      !*******************
      AE(i)   =-(gamma_e*delta_y*delta_z/de*DMAX1(cero,(1._DBL-0.1_DBL*DABS(ue*de/gamma_e))**5._DBL)+&
      &DMAX1(cero,-ue*delta_y*delta_z))
      AW(i-1) =-(gamma_w*delta_y*delta_z/dw*DMAX1(cero,(1._DBL-0.1_DBL*DABS(uw*dw/gamma_w))**5._DBL)+&
      &DMAX1(cero,uw*delta_y*delta_z))
      alpha   = gamma_n*delta_x*delta_z/dn*DMAX1(cero,(1._DBL-0.1_DBL*DABS(vn*dn/gamma_n))**5._DBL)+&
      &DMAX1(cero,-vn*delta_x*delta_z)
      beta    = gamma_s*delta_x*delta_z/ds*DMAX1(cero,(1._DBL-0.1_DBL*DABS(vs*ds/gamma_s))**5._DBL)+&
      &DMAX1(cero,vs*delta_x*delta_z)
      kappa   = gamma_top*delta_x*delta_y/dtop*DMAX1(cero,(1._DBL-0.1_DBL*DABS(wt*dtop/gamma_top))**5._DBL)+&
      &DMAX1(cero,-wt*delta_x*delta_y)
      sigma   = gamma_b*delta_x*delta_y/db*DMAX1(cero,(1._DBL-0.1_DBL*DABS(wb*db/gamma_b))**5._DBL)+&
      &DMAX1(cero, wb*delta_x*delta_y)
      APi(i)  =(-AE(i)-AW(i-1)+alpha+beta+kappa+sigma+delta_x*delta_y*delta_z/dt_o) / rel_vo
      Rxi(i)  = alpha*v_o(i,j+1,k)+beta*v_o(i,j-1,k)+kappa*v_o(i,j,k+1)+sigma*v_o(i,j,k-1)+delta_x*delta_y*delta_z*v_anto(i,j,k)/dt_o+&
      &(pres_o(i,j,k)-pres_o(i,j+1,k))*delta_x*delta_z+APi(i)*(1._DBL-rel_vo)*v_o(i,j,k)
      av_o(i,j,k) = APi(i) * rel_vo
      !****************************
    END DO
    !***********************
    !Condiciones de frontera
    APi(1)      = 1._DBL !* AW(1)
    AE(1)       = cero
    Rxi(1)      = cero !* AW(1)
    av_o(1,j,k) = 1.e40_DBL !APi(1)
    !*******************************
    APi(mi+1)      = 1._DBL !* AE(mi)
    AW(mi)         = cero
    Rxi(mi+1)      = cero !* AE(mi)
    av_o(mi+1,j,k) = 1.e40_DBL !APi(mi+1)
    !***************************
    CALL tri(AW,APi,AE,Rxi,mi+1)
    DO i = 1, mi+1
      v_o(i,j,k) = Rxi(i)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO
!********************
!TDMA en direccion  y
!$OMP  PARALLEL DO DEFAULT(NONE)&
!$OMP& PRIVATE(AS,APj,AN,Ryj,uw,ue,vs,vn,wb,wt,dw,de,ds,dn,db,dtop,gamma_w,gamma_e,gamma_s,gamma_n,gamma_b,gamma_top,delta_x,delta_y,delta_z,alpha,beta,kappa,sigma)&
!$OMP& SHARED(v_o,fv,v_anto,pres_o,gamma_vo,av_o,u_o,w_o,d_yvo,d2_yvo,xo,zo,fex,fez,d_xuo,d_zwo,rel_vo,dt_o)
DO k = 2, lk
  DO i = 2, mi
    !*************************
    !Llenado de la matriz en x
    DO j = 2, nj-1
      !**************************
      !Interpolaciones necesarias
      ue = u_o(i,j,k)+d_yvo(j-1)/d2_yvo(j)*(u_o(i,j+1,k)-u_o(i,j,k))
      uw = u_o(i-1,j,k)+d_yvo(j-1)/d2_yvo(j)*(u_o(i-1,j+1,k)-u_o(i-1,j,k))
      vn = (v_o(i,j+1,k)+v_o(i,j,k))/2._DBL
      vs = (v_o(i,j,k)+v_o(i,j-1,k))/2._DBL
      wt = w_o(i,j,k)  +d_yvo(j-1)/d2_yvo(j)*(w_o(i,j+1,k)-w_o(i,j,k))
      wb = w_o(i,j,k-1)+d_yvo(j-1)/d2_yvo(j)*(w_o(i,j+1,k-1)-w_o(i,j,k-1))
      de = xo(i+1)-xo(i)
      dw = xo(i)-xo(i-1)
      dn = d_yvo(j)
      ds = d_yvo(j-1)
      dtop = zo(k+1)-zo(k)
      db = zo(k)-zo(k-1)
      gamma_e = 1._DBL/((1._DBL-fex(i+1)) / gamma_vo(i,j,k) + fex(i+1) / gamma_vo(i+1,j,k))
      gamma_w = 1._DBL/((1._DBL-fex(i)) / gamma_vo(i-1,j,k) + fex(i) / gamma_vo(i,j,k))
      gamma_n = 2._DBL * gamma_vo(i,j+1,k) * gamma_vo(i,j,k) / (gamma_vo(i,j+1,k) + gamma_vo(i,j,k))
      gamma_s = 2._DBL * gamma_vo(i,j-1,k) * gamma_vo(i,j,k) / (gamma_vo(i,j-1,k) + gamma_vo(i,j,k))
      gamma_top = 1._DBL/((1._DBL-fez(k+1)) / gamma_vo(i,j,k) + fez(k+1) / gamma_vo(i,j,k+1))
      gamma_b = 1._DBL/((1._DBL-fez(k)) / gamma_vo(i,j,k-1) + fez(k) / gamma_vo(i,j,k))
      delta_x = d_xuo(i-1)
      delta_y = d2_yvo(j)/2._DBL
      delta_z = d_zwo(k-1)
      !*******************
      alpha   = gamma_e*delta_y*delta_z/de*DMAX1(cero,(1._DBL-0.1_DBL*DABS(ue*de/gamma_e))**5._DBL)+&
      &DMAX1(cero,-ue*delta_y*delta_z)
      beta    = gamma_w*delta_y*delta_z/dw*DMAX1(cero,(1._DBL-0.1_DBL*DABS(uw*dw/gamma_w))**5._DBL)+&
      &DMAX1(cero,uw*delta_y*delta_z)
      AN(j)   =-(gamma_n*delta_x*delta_z/dn*DMAX1(cero,(1._DBL-0.1_DBL*DABS(vn*dn/gamma_n))**5._DBL)+&
      &DMAX1(cero,-vn*delta_x*delta_z))
      AS(j-1) =-(gamma_s*delta_x*delta_z/ds*DMAX1(cero,(1._DBL-0.1_DBL*DABS(vs*ds/gamma_s))**5._DBL)+&
      &DMAX1(cero,vs*delta_x*delta_z))
      kappa   = gamma_top*delta_x*delta_y/dtop*DMAX1(cero,(1._DBL-0.1_DBL*DABS(wt*dtop/gamma_top))**5._DBL)+&
      &DMAX1(cero,-wt*delta_x*delta_y)
      sigma   = gamma_b*delta_x*delta_y/db*DMAX1(cero,(1._DBL-0.1_DBL*DABS(wb*db/gamma_b))**5._DBL)+&
      &DMAX1(cero, wb*delta_x*delta_y)
      APj(j)  =(alpha+beta-AN(j)-AS(j-1)+kappa+sigma+delta_x*delta_y*delta_z/dt_o) / rel_vo
      Ryj(j)  = alpha*v_o(i+1,j,k)+beta*v_o(i-1,j,k)+kappa*v_o(i,j,k+1)+sigma*v_o(i,j,k-1)+delta_x*delta_y*delta_z*v_anto(i,j,k)/dt_o+&
      &(pres_o(i,j,k)-pres_o(i,j+1,k))*delta_x*delta_z+APj(j)*(1._DBL-rel_vo)*v_o(i,j,k)
      av_o(i,j,k) = APj(j) * rel_vo
      !****************************
    END DO
    !***********************
    !Condiciones de frontera
    APj(1)      = 1._DBL !* AS(1)
    AN(1)       = cero
    Ryj(1)      = cero !* AS(1)
    av_o(i,1,k) = 1.e40_DBL !APj(1)
    !*******************************
    APj(nj)      = 1._DBL !* AN(nj-1)
    AS(nj-1)     = cero
    Ryj(nj)      = cero !* AN(nj-1)
    av_o(i,nj,k) = 1.e40_DBL !APj(nj)
    !*************************
    CALL tri(AS,APj,AN,Ryj,nj)
    DO j = 1, nj
      v_o(i,j,k) = Ryj(j)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO
!********************
!TDMA en direccion  z
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(AB,APk,AT,Rzk,uw,ue,vs,vn,wb,wt,dw,de,ds,dn,db,dtop,gamma_w,gamma_e,gamma_s,gamma_n,gamma_b,gamma_top,delta_x,delta_y,delta_z,alpha,beta,kappa,sigma) SHARED(v_o,fv,v_anto,pres_o,gamma_vo,av_o,u_o,w_o,d_yvo,d2_yvo,xo,zo,fex,fez,d_xuo,d_zwo,rel_vo,dt_o)
DO j = 2, nj-1
  DO i = 2, mi
    !*************************
    !Llenado de la matriz en x
    DO k = 2, lk
      !**************************
      !Interpolaciones necesarias
      ue = u_o(i,j,k)+d_yvo(j-1)/d2_yvo(j)*(u_o(i,j+1,k)-u_o(i,j,k))
      uw = u_o(i-1,j,k)+d_yvo(j-1)/d2_yvo(j)*(u_o(i-1,j+1,k)-u_o(i-1,j,k))
      vn = (v_o(i,j+1,k)+v_o(i,j,k))/2._DBL
      vs = (v_o(i,j,k)+v_o(i,j-1,k))/2._DBL
      wt = w_o(i,j,k)  +d_yvo(j-1)/d2_yvo(j)*(w_o(i,j+1,k)-w_o(i,j,k))
      wb = w_o(i,j,k-1)+d_yvo(j-1)/d2_yvo(j)*(w_o(i,j+1,k-1)-w_o(i,j,k-1))
      de = xo(i+1)-xo(i)
      dw = xo(i)-xo(i-1)
      dn = d_yvo(j)
      ds = d_yvo(j-1)
      dtop = zo(k+1)-zo(k)
      db = zo(k)-zo(k-1)
      gamma_e = 1._DBL/((1._DBL-fex(i+1)) / gamma_vo(i,j,k) + fex(i+1) / gamma_vo(i+1,j,k))
      gamma_w = 1._DBL/((1._DBL-fex(i)) / gamma_vo(i-1,j,k) + fex(i) / gamma_vo(i,j,k))
      gamma_n = 2._DBL * gamma_vo(i,j+1,k) * gamma_vo(i,j,k) / (gamma_vo(i,j+1,k) + gamma_vo(i,j,k))
      gamma_s = 2._DBL * gamma_vo(i,j-1,k) * gamma_vo(i,j,k) / (gamma_vo(i,j-1,k) + gamma_vo(i,j,k))
      gamma_top = 1._DBL/((1._DBL-fez(k+1)) / gamma_vo(i,j,k) + fez(k+1) / gamma_vo(i,j,k+1))
      gamma_b = 1._DBL/((1._DBL-fez(k)) / gamma_vo(i,j,k-1) + fez(k) / gamma_vo(i,j,k))
      delta_x = d_xuo(i-1)
      delta_y = d2_yvo(j)/2._DBL
      delta_z = d_zwo(k-1)
      !*******************
      alpha   = gamma_e*delta_y*delta_z/de*DMAX1(cero,(1._DBL-0.1_DBL*DABS(ue*de/gamma_e))**5._DBL)+&
      &DMAX1(cero,-ue*delta_y*delta_z)
      beta    = gamma_w*delta_y*delta_z/dw*DMAX1(cero,(1._DBL-0.1_DBL*DABS(uw*dw/gamma_w))**5._DBL)+&
      &DMAX1(cero,uw*delta_y*delta_z)
      kappa   = gamma_n*delta_x*delta_z/dn*DMAX1(cero,(1._DBL-0.1_DBL*DABS(vn*dn/gamma_n))**5._DBL)+&
      &DMAX1(cero,-vn*delta_x*delta_z)
      sigma   = gamma_s*delta_x*delta_z/ds*DMAX1(cero,(1._DBL-0.1_DBL*DABS(vs*ds/gamma_s))**5._DBL)+&
      &DMAX1(cero,vs*delta_x*delta_z)
      AT(k)   =-(gamma_top*delta_x*delta_y/dtop*DMAX1(cero,(1._DBL-0.1_DBL*DABS(wt*dtop/gamma_top))**5._DBL)+&
      &DMAX1(cero,-wt*delta_x*delta_y))
      AB(k-1) =-(gamma_b*delta_x*delta_y/db*DMAX1(cero,(1._DBL-0.1_DBL*DABS(wb*db/gamma_b))**5._DBL)+&
      &DMAX1(cero, wb*delta_x*delta_y))
      APk(k)  =(alpha+beta+kappa+sigma-AT(k)-AB(k-1)+delta_x*delta_y*delta_z/dt_o) / rel_vo
      Rzk(k)  = alpha*v_o(i+1,j,k)+beta*v_o(i-1,j,k)+kappa*v_o(i,j+1,k)+sigma*v_o(i,j-1,k)+delta_x*delta_y*delta_z*v_anto(i,j,k)/dt_o+&
      &(pres_o(i,j,k)-pres_o(i,j+1,k))*delta_x*delta_z+APk(k)*(1._DBL-rel_vo)*v_o(i,j,k)
      av_o(i,j,k) = APk(k) * rel_vo
      !****************************
    END DO
    !***********************
    !Condiciones de frontera
    APk(1)      = 1._DBL !* AB(1)
    AT(1)       = cero
    Rzk(1)      = cero !* AB(1)
    av_o(i,j,1) = 1.e40_DBL !APk(1)
    !*******************************
    APk(lk+1)      = 1._DBL !* AT(lk)
    AB(lk)         = cero
    Rzk(lk+1)      = cero !* AT(lk)
    av_o(i,j,lk+1) = 1.e40_DBL !APk(lk+1)
    !***************************
    CALL tri(AB,APk,AT,Rzk,lk+1)
    DO k = 1, lk+1
      v_o(i,j,k) = Rzk(k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO
!************
dv_o = v_o-fv
END SUBROUTINE
