SUBROUTINE vel_w(xo,yo,zwo,fex,fey,d_zwo,d2_zwo,d_xuo,d_yvo,w_o,w_anto,u_o,v_o,pres_o,temp_o,gamma_wo,Ri_o,dt_o,dw_o,aw_o,rel_vo)
USE constantes
! USE mkl95_LAPACK
IMPLICIT NONE
INTEGER :: i,j,k,info
!**********************************
! Variables del flujo e inc'ognitas
REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk),   INTENT(inout) :: w_o
REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk),   INTENT(out)   :: aw_o,dw_o
REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk),   INTENT(in)    :: w_anto,gamma_wo,Ri_o
REAL(kind=DBL), DIMENSION(mi,nj+1,lk+1),   INTENT(in)    :: u_o
REAL(kind=DBL), DIMENSION(mi+1,nj,lk+1),   INTENT(in)    :: v_o
REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk+1), INTENT(in)    :: pres_o,temp_o
REAL(kind=DBL), DIMENSION(mi+1,nj+1,lk)                  :: fw
!***************************
!Variables de la tridiagonal
REAL(kind=DBL), DIMENSION(mi+1) :: APi,Rxi
REAL(kind=DBL), DIMENSION(mi)   :: AW,AE
REAL(kind=DBL), DIMENSION(nj+1) :: APj,Ryj
REAL(kind=DBL), DIMENSION(nj)   :: AS,AN
REAL(kind=DBL), DIMENSION(lk)   :: APk,Rzk
REAL(kind=DBL), DIMENSION(lk-1) :: AB,AT
!***************************
!Variables de interpolaci'on
REAL(kind=DBL) :: uw,ue,vs,vn,wb,wt,dw,de,ds,dn,db,dtop,gamma_w,gamma_e,gamma_s,gamma_n,gamma_b,gamma_top
REAL(kind=DBL) :: delta_x,delta_y,delta_z,temp_int
!*************************************************************
!Variables de la malla,volumen de control,incremento de tiempo
REAL(kind=DBL), DIMENSION(lk),   INTENT(in) :: zwo
REAL(kind=DBL), DIMENSION(mi+1), INTENT(in) :: xo,fex
REAL(kind=DBL), DIMENSION(nj+1), INTENT(in) :: yo,fey
REAL(kind=DBL), DIMENSION(lk-1), INTENT(in) :: d_zwo,d2_zwo
REAL(kind=DBL), DIMENSION(mi-1), INTENT(in) :: d_xuo
REAL(kind=DBL), DIMENSION(nj-1), INTENT(in) :: d_yvo
REAL(kind=DBL), INTENT(in)                  :: dt_o,rel_vo
!********************
!Variables auxiliares
REAL(kind=DBL) :: alpha,beta,kappa,sigma
!*************************
!auxiliar para calcular du
fw = w_o
!********************
!TDMA en direccion x
!$OMP  PARALLEL DO DEFAULT(NONE)&
!$OMP& PRIVATE(AW,APi,AE,Rxi,uw,ue,vs,vn,wb,wt,dw,de,ds,dn,db,dtop,gamma_w,gamma_e,gamma_s,gamma_n,gamma_b,gamma_top,delta_x,delta_y,delta_z,alpha,beta,kappa,sigma,temp_int)&
!$OMP& SHARED(w_o,w_anto,pres_o,temp_o,gamma_wo,aw_o,u_o,v_o,d_zwo,d2_zwo,xo,yo,fex,fey,d_xuo,d_yvo,Ri_o,rel_vo,dt_o)
DO k = 2, lk-1
  DO j = 2, nj
    !*************************
    !Llenado de la matriz en x
    DO i = 2, mi
      !**************************
      !Interpolaciones necesarias
      ue = u_o(i,j,k)+d_zwo(k-1)/d2_zwo(k)*(u_o(i,j,k+1)-u_o(i,j,k))
      uw = u_o(i-1,j,k)+d_zwo(k-1)/d2_zwo(k)*(u_o(i-1,j,k+1)-u_o(i-1,j,k))
      vn = v_o(i,j,k)+d_zwo(k-1)/d2_zwo(k)*(v_o(i,j,k+1)-v_o(i,j,k))
      vs = v_o(i,j-1,k)+d_zwo(k-1)/d2_zwo(k)*(v_o(i,j-1,k+1)-v_o(i,j-1,k))
      wt = (w_o(i,j,k)+w_o(i,j,k+1))/2._DBL
      wb = (w_o(i,j,k-1)+w_o(i,j,k))/2._DBL
      de = xo(i+1)-xo(i)
      dw = xo(i)-xo(i-1)
      dn = yo(j+1)-yo(j)
      ds = yo(j)-yo(j-1)
      dtop = d_zwo(k)
      db = d_zwo(k-1)
      gamma_e = 1._DBL/((1._DBL-fex(i+1)) / gamma_wo(i,j,k) + fex(i+1) / gamma_wo(i+1,j,k))
      gamma_w = 1._DBL/((1._DBL-fex(i)) / gamma_wo(i-1,j,k) + fex(i) / gamma_wo(i,j,k))
      gamma_n = 1._DBL/((1._DBL-fey(j+1)) / gamma_wo(i,j,k) + fey(j+1) / gamma_wo(i,j+1,k))
      gamma_s = 1._DBL/((1._DBL-fey(j)) / gamma_wo(i,j-1,k) + fey(j) / gamma_wo(i,j,k))
      gamma_top = 2._DBL * gamma_wo(i,j,k+1) * gamma_wo(i,j,k) / (gamma_wo(i,j,k+1) + gamma_wo(i,j,k))
      gamma_b = 2._DBL * gamma_wo(i,j,k-1) * gamma_wo(i,j,k) / (gamma_wo(i,j,k-1) + gamma_wo(i,j,k))
      delta_x = d_xuo(i-1)
      delta_y = d_yvo(j-1)
      delta_z = d2_zwo(k)/2._DBL
      temp_int = (temp_o(i,j,k)+temp_o(i,j,k+1))/2._DBL
      !****************************
      AE(i)   =-(gamma_e*delta_y*delta_z/de*DMAX1(cero,(1._DBL-0.1_DBL*DABS(ue*de/gamma_e))**5._DBL)+&
      &DMAX1(cero,-ue*delta_y*delta_z))
      AW(i-1) =-(gamma_w*delta_y*delta_z/dw*DMAX1(cero,(1._DBL-0.1_DBL*DABS(uw*dw/gamma_w))**5._DBL)+&
      &DMAX1(cero,uw*delta_y*delta_z))
      alpha   = gamma_n*delta_x*delta_z/dn*DMAX1(cero,(1._DBL-0.1_DBL*DABS(vn*dn/gamma_n))**5._DBL)+&
      &DMAX1(cero,-vn*delta_x*delta_z)
      beta    = gamma_s*delta_x*delta_z/ds*DMAX1(cero,(1._DBL-0.1_DBL*DABS(vs*ds/gamma_s))**5._DBL)+&
      &DMAX1(cero, vs*delta_x*delta_z)
      kappa   = gamma_top*delta_x*delta_y/dtop*DMAX1(cero,(1._DBL-0.1_DBL*DABS(wt*dtop/gamma_top))**5._DBL)+&
      &DMAX1(cero,-wt*delta_x*delta_y)
      sigma   = gamma_b*delta_x*delta_y/db*DMAX1(cero,(1._DBL-0.1_DBL*DABS(wb*db/gamma_b))**5._DBL)+&
      &DMAX1(cero, wb*delta_x*delta_y)
      APi(i)  =(-AE(i)-AW(i-1)+alpha+beta+kappa+sigma+delta_x*delta_y*delta_z/dt_o) / rel_vo
      Rxi(i)  = alpha*w_o(i,j+1,k)+beta*w_o(i,j-1,k)+kappa*w_o(i,j,k+1)+sigma*w_o(i,j,k-1)+delta_x*delta_y*delta_z*w_anto(i,j,k)/dt_o+&
      &(pres_o(i,j,k)-pres_o(i,j,k+1))*delta_x*delta_y+APi(i)*(1._DBL-rel_vo)*w_o(i,j,k)+Ri_o(i,j,k)*temp_int*delta_x*delta_y*delta_z
      aw_o(i,j,k) = APi(i) * rel_vo
      !****************************
    END DO
    !***********************
    !Condiciones de frontera
    APi(1)      = 1._DBL !* AW(1)
    AE(1)       = cero
    Rxi(1)      = cero !* AW(1)
    aw_o(1,j,k) = 1.e40_DBL !APi(1)
    !***********************
    APi(mi+1)      = 1._DBL !* AE(mi)
    AW(mi)         = cero
    Rxi(mi+1)      = cero !* AE(mi)
    aw_o(mi+1,j,k) = 1.e40_DBL !APi(mi+1)
    !***********************
    CALL tri(AW,APi,AE,Rxi,mi+1)
    DO i = 1, mi+1
      w_o(i,j,k) = Rxi(i)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO
!********************
!TDMA en direccion  y
!$OMP  PARALLEL DO DEFAULT(NONE)&
!$OMP& PRIVATE(AS,APj,AN,Ryj,uw,ue,vs,vn,wb,wt,dw,de,ds,dn,db,dtop,gamma_w,gamma_e,gamma_s,gamma_n,gamma_b,gamma_top,delta_x,delta_y,delta_z,alpha,beta,kappa,sigma,temp_int)&
!$OMP& SHARED(w_o,w_anto,pres_o,temp_o,gamma_wo,aw_o,u_o,v_o,d_zwo,d2_zwo,xo,yo,fex,fey,d_xuo,d_yvo,Ri_o,rel_vo,dt_o)
DO k = 2, lk-1
  DO i = 2, mi
    !*************************
    !Llenado de la matriz en y
    DO j = 2, nj
      !**************************
      !Interpolaciones necesarias
      ue = u_o(i,j,k)+d_zwo(k-1)/d2_zwo(k)*(u_o(i,j,k+1)-u_o(i,j,k))
      uw = u_o(i-1,j,k)+d_zwo(k-1)/d2_zwo(k)*(u_o(i-1,j,k+1)-u_o(i-1,j,k))
      vn = v_o(i,j,k)+d_zwo(k-1)/d2_zwo(k)*(v_o(i,j,k+1)-v_o(i,j,k))
      vs = v_o(i,j-1,k)+d_zwo(k-1)/d2_zwo(k)*(v_o(i,j-1,k+1)-v_o(i,j-1,k))
      wt = (w_o(i,j,k)+w_o(i,j,k+1))/2._DBL
      wb = (w_o(i,j,k-1)+w_o(i,j,k))/2._DBL
      de = xo(i+1)-xo(i)
      dw = xo(i)-xo(i-1)
      dn = yo(j+1)-yo(j)
      ds = yo(j)-yo(j-1)
      dtop = d_zwo(k)
      db = d_zwo(k-1)
      gamma_e = 1._DBL/((1._DBL-fex(i+1)) / gamma_wo(i,j,k) + fex(i+1) / gamma_wo(i+1,j,k))
      gamma_w = 1._DBL/((1._DBL-fex(i)) / gamma_wo(i-1,j,k) + fex(i) / gamma_wo(i,j,k))
      gamma_n = 1._DBL/((1._DBL-fey(j+1)) / gamma_wo(i,j,k) + fey(j+1) / gamma_wo(i,j+1,k))
      gamma_s = 1._DBL/((1._DBL-fey(j)) / gamma_wo(i,j-1,k) + fey(j) / gamma_wo(i,j,k))
      gamma_top = 2._DBL * gamma_wo(i,j,k+1) * gamma_wo(i,j,k) / (gamma_wo(i,j,k+1) + gamma_wo(i,j,k))
      gamma_b = 2._DBL * gamma_wo(i,j,k-1) * gamma_wo(i,j,k) / (gamma_wo(i,j,k-1) + gamma_wo(i,j,k))
      delta_x = d_xuo(i-1)
      delta_y = d_yvo(j-1)
      delta_z = d2_zwo(k)/2._DBL
      temp_int = (temp_o(i,j,k)+temp_o(i,j,k+1))/2._DBL
      !****************************
      alpha   = gamma_e*delta_y*delta_z/de*DMAX1(cero,(1._DBL-0.1_DBL*DABS(ue*de/gamma_e))**5._DBL)+&
      &DMAX1(cero,-ue*delta_y*delta_z)
      beta    = gamma_w*delta_y*delta_z/dw*DMAX1(cero,(1._DBL-0.1_DBL*DABS(uw*dw/gamma_w))**5._DBL)+&
      &DMAX1(cero,uw*delta_y*delta_z)
      AN(j)   =-(gamma_n*delta_x*delta_z/dn*DMAX1(cero,(1._DBL-0.1_DBL*DABS(vn*dn/gamma_n))**5._DBL)+&
      &DMAX1(cero,-vn*delta_x*delta_z))
      AS(j-1) =-(gamma_s*delta_x*delta_z/ds*DMAX1(cero,(1._DBL-0.1_DBL*DABS(vs*ds/gamma_s))**5._DBL)+&
      &DMAX1(cero, vs*delta_x*delta_z))
      kappa   = gamma_top*delta_x*delta_y/dtop*DMAX1(cero,(1._DBL-0.1_DBL*DABS(wt*dtop/gamma_top))**5._DBL)+&
      &DMAX1(cero,-wt*delta_x*delta_y)
      sigma   = gamma_b*delta_x*delta_y/db*DMAX1(cero,(1._DBL-0.1_DBL*DABS(wb*db/gamma_b))**5._DBL)+&
      &DMAX1(cero, wb*delta_x*delta_y)
      APj(j)  =(alpha+beta-AN(j)-AS(j-1)+kappa+sigma+delta_x*delta_y*delta_z/dt_o) / rel_vo
      Ryj(j)  = alpha*w_o(i+1,j,k)+beta*w_o(i-1,j,k)+kappa*w_o(i,j,k+1)+sigma*w_o(i,j,k-1)+delta_x*delta_y*delta_z*w_anto(i,j,k)/dt_o+&
      &(pres_o(i,j,k)-pres_o(i,j,k+1))*delta_x*delta_y+APj(j)*(1._DBL-rel_vo)*w_o(i,j,k)+Ri_o(i,j,k)*temp_int*delta_x*delta_y*delta_z
      aw_o(i,j,k) = APj(j) * rel_vo
      !****************************
    END DO
    !***********************
    !Condiciones de frontera
    APj(1)      = 1._DBL !* AS(1)
    AN(1)       = cero
    Ryj(1)      = cero !* AS(1)
    aw_o(i,1,k) = APj(1)
    !***********************
    APj(nj+1)      = 1._DBL !* AN(nj)
    AS(nj)         = cero
    Ryj(nj+1)      = cero !* AN(nj)
    aw_o(i,nj+1,k) = APj(nj+1)
    !***********************
    CALL tri(AS,APj,AN,Ryj,nj+1)
    DO j = 1, nj+1
      w_o(i,j,k) = Ryj(j)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO
!********************
!TDMA en direccion  z
!$OMP  PARALLEL DO DEFAULT(NONE)&
!$OMP& PRIVATE(AB,APk,AT,Rzk,uw,ue,vs,vn,wb,wt,dw,de,ds,dn,db,dtop,gamma_w,gamma_e,gamma_s,gamma_n,gamma_b,gamma_top,delta_x,delta_y,delta_z,alpha,beta,kappa,sigma,temp_int)&
!$OMP& SHARED(w_o,w_anto,pres_o,temp_o,gamma_wo,aw_o,u_o,v_o,d_zwo,d2_zwo,xo,yo,fex,fey,d_xuo,d_yvo,Ri_o,rel_vo,dt_o,zwo,fw)
DO j = 2, nj
  DO i = 2, mi
    !*************************
    !Llenado de la matriz en z
    DO k = 2, lk-1
      !**************************
      !Interpolaciones necesarias
      ue = u_o(i,j,k)+d_zwo(k-1)/d2_zwo(k)*(u_o(i,j,k+1)-u_o(i,j,k))
      uw = u_o(i-1,j,k)+d_zwo(k-1)/d2_zwo(k)*(u_o(i-1,j,k+1)-u_o(i-1,j,k))
      vn = v_o(i,j,k)+d_zwo(k-1)/d2_zwo(k)*(v_o(i,j,k+1)-v_o(i,j,k))
      vs = v_o(i,j-1,k)+d_zwo(k-1)/d2_zwo(k)*(v_o(i,j-1,k+1)-v_o(i,j-1,k))
      wt = (w_o(i,j,k)+w_o(i,j,k+1))/2._DBL
      wb = (w_o(i,j,k-1)+w_o(i,j,k))/2._DBL
      de = xo(i+1)-xo(i)
      dw = xo(i)-xo(i-1)
      dn = yo(j+1)-yo(j)
      ds = yo(j)-yo(j-1)
      dtop = d_zwo(k)
      db = d_zwo(k-1)
      gamma_e = 1._DBL/((1._DBL-fex(i+1)) / gamma_wo(i,j,k) + fex(i+1) / gamma_wo(i+1,j,k))
      gamma_w = 1._DBL/((1._DBL-fex(i)) / gamma_wo(i-1,j,k) + fex(i) / gamma_wo(i,j,k))
      gamma_n = 1._DBL/((1._DBL-fey(j+1)) / gamma_wo(i,j,k) + fey(j+1) / gamma_wo(i,j+1,k))
      gamma_s = 1._DBL/((1._DBL-fey(j)) / gamma_wo(i,j-1,k) + fey(j) / gamma_wo(i,j,k))
      gamma_top = 2._DBL * gamma_wo(i,j,k+1) * gamma_wo(i,j,k) / (gamma_wo(i,j,k+1) + gamma_wo(i,j,k))
      gamma_b = 2._DBL * gamma_wo(i,j,k-1) * gamma_wo(i,j,k) / (gamma_wo(i,j,k-1) + gamma_wo(i,j,k))
      delta_x = d_xuo(i-1)
      delta_y = d_yvo(j-1)
      delta_z = d2_zwo(k)/2._DBL
      temp_int = (temp_o(i,j,k)+temp_o(i,j,k+1))/2._DBL
      !****************************
      alpha   = gamma_e*delta_y*delta_z/de*DMAX1(cero,(1._DBL-0.1_DBL*DABS(ue*de/gamma_e))**5._DBL)+&
      &DMAX1(cero,-ue*delta_y*delta_z)
      beta    = gamma_w*delta_y*delta_z/dw*DMAX1(cero,(1._DBL-0.1_DBL*DABS(uw*dw/gamma_w))**5._DBL)+&
      &DMAX1(cero,uw*delta_y*delta_z)
      kappa   = gamma_n*delta_x*delta_z/dn*DMAX1(cero,(1._DBL-0.1_DBL*DABS(vn*dn/gamma_n))**5._DBL)+&
      &DMAX1(cero,-vn*delta_x*delta_z)
      sigma   = gamma_s*delta_x*delta_z/ds*DMAX1(cero,(1._DBL-0.1_DBL*DABS(vs*ds/gamma_s))**5._DBL)+&
      &DMAX1(cero, vs*delta_x*delta_z)
      AT(k)   =-(gamma_top*delta_x*delta_y/dtop*DMAX1(cero,(1._DBL-0.1_DBL*DABS(wt*dtop/gamma_top))**5._DBL)+&
      &DMAX1(cero,-wt*delta_x*delta_y))
      AB(k-1) =-(gamma_b*delta_x*delta_y/db*DMAX1(cero,(1._DBL-0.1_DBL*DABS(wb*db/gamma_b))**5._DBL)+&
      &DMAX1(cero, wb*delta_x*delta_y))
      APk(k)  =(alpha+beta+kappa+sigma-AT(k)-AB(k-1)+delta_x*delta_y*delta_z/dt_o) / rel_vo
      Rzk(k)  = alpha*w_o(i+1,j,k)+beta*w_o(i-1,j,k)+kappa*w_o(i,j+1,k)+sigma*w_o(i,j-1,k)+delta_x*delta_y*delta_z*w_anto(i,j,k)/dt_o+&
      &(pres_o(i,j,k)-pres_o(i,j,k+1))*delta_x*delta_y+APk(k)*(1._DBL-rel_vo)*w_o(i,j,k)+Ri_o(i,j,k)*temp_int*delta_x*delta_y*delta_z
      aw_o(i,j,k) = APk(k) * rel_vo
      !****************************
    END DO
    !***********************
    !Condiciones de frontera
    APk(1)      = 1._DBL !* AB(1)
    AT(1)       = 0._DBL !cero
    Rzk(1)      = 0._DBL
    ! (fw(i,j,3)*(zwo(2)-zwo(1))*(zwo(2)-zwo(1))/((zwo(3)-zwo(2))*(2*zwo(1)-zwo(2)-zwo(3)))+&
    !              &fw(i,j,2)*(1._DBL-(zwo(2)-zwo(1))*(zwo(2)-zwo(1))/((zwo(3)-zwo(2))*(2._DBL*zwo(1)-zwo(2)-zwo(3))))) !*AB(1) !cero * AB(1)
    aw_o(i,j,1) = 1.e40_DBL !APk(1)
    !*******************************
    APk(lk)      = 1._DBL !* AT(lk-1)
    AB(lk-1)     = 0._DBL !cero
    Rzk(lk)      = 0._DBL !-1._DBL
!                    (fw(i,j,lk-1)*((zwo(lk)-zwo(lk-2))*(zwo(lk)-zwo(lk-1))/((zwo(lk-1)-zwo(lk-2))*(2*zwo(lk)-zwo(lk-2)-zwo(lk-1)))+&
!                    (zwo(lk)-zwo(lk-2))/(2*zwo(lk)-zwo(lk-2)-zwo(lk-1)))-fw(i,j,lk-2)*((zwo(lk)-zwo(lk-1))*(zwo(lk)-zwo(lk-1))/&
!                    ((zwo(lk-1)-zwo(lk-2))*(2._DBL*zwo(lk)-zwo(lk-2)-zwo(lk-1))))) !* AT(lk-1)
    aw_o(i,j,lk) = 1.e40_DBL !APk(lk)
    !***************************
    CALL tri(AB,APk,AT,Rzk,lk)
    DO k = 1, lk
      w_o(i,j,k) = Rzk(k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO
!************
dw_o = w_o-fw
END SUBROUTINE
