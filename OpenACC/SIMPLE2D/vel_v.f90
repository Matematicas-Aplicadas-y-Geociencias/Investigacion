SUBROUTINE vel_v(xo,fex,d_yvo,d2_yvo,d_xuo,u_o,v_o,v_anto,pres_o,gamma_vo,dt_o,dv_o,av_o,rel_vo)
USE malla
! USE mkl95_LAPACK
IMPLICIT NONE
INTEGER :: i,j,k,info
!**********************************
! Variables del flujo e inc'ognitas
REAL(kind=DBL), DIMENSION(mi+1,nj), INTENT(inout) :: v_o
REAL(kind=DBL), DIMENSION(mi+1,nj), INTENT(out)   :: av_o,dv_o
REAL(kind=DBL), DIMENSION(mi+1,nj), INTENT(in)    :: v_anto,gamma_vo
REAL(kind=DBL), DIMENSION(mi,nj+1), INTENT(in)    :: u_o
REAL(kind=DBL), DIMENSION(mi+1,nj+1), INTENT(in)  :: pres_o
REAL(kind=DBL), DIMENSION(mi+1,nj) :: fv
!**********************************
!Variables de la tridiagonal
REAL(kind=DBL), DIMENSION(mi+1) :: ACi,Rxi
REAL(kind=DBL), DIMENSION(mi)   :: AI,AD
REAL(kind=DBL), DIMENSION(nj)   :: ACj,Ryj
REAL(kind=DBL), DIMENSION(nj-1) :: AS,AN
!*********************
!Variables de interpolaci'on
REAL(kind=DBL) :: ui,ud,vn,vs,di,dd,ds,dn,gamma_i,gamma_d,gamma_s,gamma_n,delta_x,delta_y
!*********************
!Variables de la malla,volumen de control,incremento de tiempo y num Richardson
REAL(kind=DBL), DIMENSION(mi+1), INTENT(in) :: xo,fex
REAL(kind=DBL), DIMENSION(mi-1), INTENT(in) :: d_xuo
REAL(kind=DBL), DIMENSION(nj-1), INTENT(in) :: d_yvo,d2_yvo
REAL(kind=DBL), INTENT(in) :: dt_o,rel_vo
!********************
!Variables auxiliares
REAL(kind=DBL) :: alpha,beta
!*************************
!auxiliar para calcular dv
fv = v_o
!********************
!TDMA en direccion x
!$OMP  PARALLEL DO DEFAULT(NONE)&
!$OMP& PRIVATE(AI,ACi,AD,Rxi,ui,ud,vs,vn,di,dd,ds,dn,gamma_i,gamma_d,gamma_s,gamma_n,delta_x,delta_y,alpha,beta,info)&
!$OMP& SHARED(v_o,fv,v_anto,pres_o,gamma_vo,av_o,u_o,d_xuo,d2_yvo,xo,fex,d_yvo,rel_vo,dt_o)
DO j = 2, nj-1
  !***********************
  !Condiciones de frontera
  ACi(1) = 1._DBL
  AD(1)  = cero
  Rxi(1) = cero
  av_o(1,j) = 1.e40_DBL !ACi(1)
  !*************************
  !Llenado de la matriz en x
  DO i = 2, mi
    !**************************
    !Interpolaciones necesarias
    ud = u_o(i,j)+d_yvo(j-1)/d2_yvo(j)*(u_o(i,j+1)-u_o(i,j))
    ui = u_o(i-1,j)+d_yvo(j-1)/d2_yvo(j)*(u_o(i-1,j+1)-u_o(i-1,j))
    vn = (v_o(i,j+1)+v_o(i,j))/2._DBL
    vs = (v_o(i,j)+v_o(i,j-1))/2._DBL
    di = xo(i)-xo(i-1)
    dd = xo(i+1)-xo(i)
    ds = d_yvo(j-1)
    dn = d_yvo(j)
    gamma_i = 1._DBL/((1._DBL-fex(i)) / gamma_vo(i-1,j) + fex(i) / gamma_vo(i,j))
    gamma_d = 1._DBL/((1._DBL-fex(i+1)) / gamma_vo(i,j) + fex(i+1) / gamma_vo(i+1,j))
    gamma_s = 2._DBL * gamma_vo(i,j-1) * gamma_vo(i,j) / (gamma_vo(i,j-1) + gamma_vo(i,j))
    gamma_n = 2._DBL * gamma_vo(i,j+1) * gamma_vo(i,j) / (gamma_vo(i,j+1) + gamma_vo(i,j))
    delta_x = d_xuo(i-1)
    delta_y = d2_yvo(j)/2._DBL
    !************************
    AI(i-1) =-(gamma_i*delta_y/di*DMAX1(cero,(1._DBL-0.1_DBL*dabs(ui*di/gamma_i))**5._DBL)+&
    &DMAX1(cero,ui*delta_y))
    AD(i)   =-(gamma_d*delta_y/dd*DMAX1(cero,(1._DBL-0.1_DBL*dabs(ud*dd/gamma_d))**5._DBL)+&
    &DMAX1(cero,-ud*delta_y))
    alpha   = gamma_s*delta_x/ds*DMAX1(cero,(1._DBL-0.1_DBL*dabs(vs*ds/gamma_s))**5._DBL)+&
    &DMAX1(cero,vs*delta_x)
    beta    = gamma_n*delta_x/dn*DMAX1(cero,(1._DBL-0.1_DBL*dabs(vn*dn/gamma_n))**5._DBL)+&
    &DMAX1(cero,-vn*delta_x)
    ACi(i)  =(-AI(i-1)-AD(i)+alpha+beta+delta_x*delta_y/dt_o) / rel_vo
    Rxi(i)  = alpha*v_o(i,j-1)+beta*v_o(i,j+1)+delta_x*delta_y*v_anto(i,j)/dt_o+&
    &(pres_o(i,j)-pres_o(i,j+1))*delta_x+ACi(i)*(1._DBL-rel_vo)*v_o(i,j)
    av_o(i,j) = ACi(i) * rel_vo
  END DO
  !***********************
  !Condiciones de frontera
  ACi(mi+1) = 1._DBL
  AI(mi)    = cero !-1._DBL
  Rxi(mi+1) = cero
  av_o(mi+1,j) = 1.e40_DBL !ACi(mi+1)
  !***************************
  CALL tri(AI,ACi,AD,Rxi,mi+1)
  DO i = 1, mi+1
    v_o(i,j) = Rxi(i)
  END DO
END DO
!$OMP END PARALLEL DO
!$OMP  PARALLEL DO DEFAULT(NONE)&
!$OMP& PRIVATE(AS,ACj,AN,Ryj,ui,ud,vs,vn,di,dd,ds,dn,gamma_i,gamma_d,gamma_s,gamma_n,delta_x,delta_y,alpha,beta,info)&
!$OMP& SHARED(v_o,v_anto,pres_o,gamma_vo,av_o,u_o,d_xuo,d2_yvo,xo,fex,d_yvo,rel_vo,dt_o)
DO i = 2, mi
  !***********************
  !Condiciones de frontera
  ACj(1) = 1._DBL
  AN(1)  = cero
  Ryj(1) = cero
  av_o(i,1) = 1.e40_DBL !ACj(1)
  !*************************
  !Llenado de la matriz en x
  DO j = 2, nj-1
    !**************************
    !Interpolaciones necesarias
    ud = u_o(i,j)+d_yvo(j-1)/d2_yvo(j)*(u_o(i,j+1)-u_o(i,j))
    ui = u_o(i-1,j)+d_yvo(j-1)/d2_yvo(j)*(u_o(i-1,j+1)-u_o(i-1,j))
    vn = (v_o(i,j+1)+v_o(i,j))/2._DBL
    vs = (v_o(i,j)+v_o(i,j-1))/2._DBL
    di = xo(i)-xo(i-1)
    dd = xo(i+1)-xo(i)
    ds = d_yvo(j-1)
    dn = d_yvo(j)
    gamma_i = 1._DBL/((1._DBL-fex(i)) / gamma_vo(i-1,j) + fex(i) / gamma_vo(i,j))
    gamma_d = 1._DBL/((1._DBL-fex(i+1)) / gamma_vo(i,j) + fex(i+1) / gamma_vo(i+1,j))
    gamma_s = 2._DBL * gamma_vo(i,j-1) * gamma_vo(i,j) / (gamma_vo(i,j-1) + gamma_vo(i,j))
    gamma_n = 2._DBL * gamma_vo(i,j+1) * gamma_vo(i,j) / (gamma_vo(i,j+1) + gamma_vo(i,j))
    delta_x = d_xuo(i-1)
    delta_y = d2_yvo(j)/2._DBL
    !************************
    alpha   = gamma_i*delta_y/di*DMAX1(cero,(1._DBL-0.1_DBL*dabs(ui*di/gamma_i))**5._DBL)+&
    &DMAX1(cero,ui*delta_y)
    beta    = gamma_d*delta_y/dd*DMAX1(cero,(1._DBL-0.1_DBL*dabs(ud*dd/gamma_d))**5._DBL)+&
    &DMAX1(cero,-ud*delta_y)
    AS(j-1) =-(gamma_s*delta_x/ds*DMAX1(cero,(1._DBL-0.1_DBL*dabs(vs*ds/gamma_s))**5._DBL)+&
    &DMAX1(cero,vs*delta_x))
    AN(j)   =-(gamma_n*delta_x/dn*DMAX1(cero,(1._DBL-0.1_DBL*dabs(vn*dn/gamma_n))**5._DBL)+&
    &DMAX1(cero,-vn*delta_x))
    ACj(j)  =(-AS(j-1)-AN(j)+alpha+beta+delta_x*delta_y/dt_o) / rel_vo
    Ryj(j)  = alpha*v_o(i-1,j)+beta*v_o(i+1,j)+delta_x*delta_y*v_anto(i,j)/dt_o+&
    &(pres_o(i,j)-pres_o(i,j+1))*delta_x+ACj(j)*(1._DBL-rel_vo)*v_o(i,j)
    av_o(i,j) = ACj(j) * rel_vo
  END DO
  !***********************
  !Condiciones de frontera
  ACj(nj)  = 1._DBL
  AS(nj-1) = cero
  Ryj(nj)  = cero
  av_o(i,nj) = 1.e40_DBL !ACj(nj) 
  !*************************
  CALL tri(AS,ACj,AN,Ryj,nj)
  DO j = 1, nj
    v_o(i,j) = Ryj(j)
  END DO
END DO
!$OMP END PARALLEL DO
WHERE(v_o /= cero)
  dv_o = (v_o-fv)/v_o
ELSEWHERE
  dv_o = v_o-fv
END WHERE
END SUBROUTINE
