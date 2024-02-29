SUBROUTINE vel_u(yo,fey,d_xuo,d2_xuo,d_yvo,u_o,u_anto,v_o,temp_o,pres_o,gamma_uo,Ri_o,dt_o,du_o,au_o,rel_vo)
USE constantes
! USE mkl95_LAPACK
IMPLICIT NONE
INTEGER :: i,j,k,info
!**********************************
! Variables del flujo e inc'ognitas
REAL(kind=DBL), DIMENSION(mi,nj+1), INTENT(inout) :: u_o
REAL(kind=DBL), DIMENSION(mi,nj+1), INTENT(out)   :: au_o,du_o
REAL(kind=DBL), DIMENSION(mi,nj+1), INTENT(in)    :: u_anto,gamma_uo,Ri_o
REAL(kind=DBL), DIMENSION(mi+1,nj), INTENT(in)    :: v_o
REAL(kind=DBL), DIMENSION(mi+1,nj+1), INTENT(in)  :: temp_o,pres_o
REAL(kind=DBL), DIMENSION(mi,nj+1) :: fu
!**********************************
!Variables de la tridiagonal
REAL(kind=DBL), DIMENSION(mi)   :: ACi,Rxi
REAL(kind=DBL), DIMENSION(mi-1) :: AI,AD
REAL(kind=DBL), DIMENSION(nj+1) :: ACj,Ryj
REAL(kind=DBL), DIMENSION(nj)   :: AS,AN
!*********************
!Variables de interpolaci'on
REAL(kind=DBL) :: ui,ud,vn,vs,di,dd,ds,dn,gamma_i,gamma_d,gamma_s,gamma_n
REAL(kind=DBL) :: delta_x,delta_y,temp_int
!*********************
!Variables de la malla,volumen de control,incremento de tiempo y num Richardson
REAL(kind=DBL), DIMENSION(nj+1), INTENT(in) :: yo,fey
REAL(kind=DBL), DIMENSION(mi-1), INTENT(in) :: d_xuo,d2_xuo
REAL(kind=DBL), DIMENSION(nj-1), INTENT(in) :: d_yvo
REAL(kind=DBL), INTENT(in) :: dt_o,rel_vo
!********************
!Variables auxiliares
REAL(kind=DBL) :: alpha,beta
!*************************
!auxiliar para calcular du
fu = u_o
!********************
!TDMA en direccion x
!$OMP  PARALLEL DO DEFAULT(NONE)&
!$OMP& PRIVATE(AI,ACi,AD,Rxi,ui,ud,vs,vn,di,dd,ds,dn,gamma_i,gamma_d,gamma_s,gamma_n,delta_x,delta_y,temp_int,alpha,beta,info)&
!$OMP& SHARED(u_o,fu,u_anto,pres_o,gamma_uo,au_o,v_o,d_xuo,d2_xuo,yo,fey,d_yvo,temp_o,rel_vo,Ri_o,dt_o)
DO j = 2, nj
  !***********************
  !Condiciones de frontera
  ACi(1) = 1._DBL
  AD(1)  = cero
  Rxi(1) = 1._DBL !-6.d0*dfloat(j-1)/dfloat(nj)*(dfloat(j-1)/dfloat(nj)-1.d0)! 1._DBL !
  au_o(1,j) = 1.e40_DBL !ACi(1)
  !*************************
  !Llenado de la matriz en x
  DO i = 2, mi-1
    !**************************
    !Interpolaciones necesarias
    ui = (u_o(i-1,j)+u_o(i,j))/2._DBL
    ud = (u_o(i,j)+u_o(i+1,j))/2._DBL
    vs = v_o(i,j-1)+d_xuo(i-1)/d2_xuo(i)*(v_o(i+1,j-1)-v_o(i,j-1))
    vn = v_o(i,j)  +d_xuo(i-1)/d2_xuo(i)*(v_o(i+1,j)-v_o(i,j))
    di = d_xuo(i-1)
    dd = d_xuo(i)
    ds = yo(j)-yo(j-1)
    dn = yo(j+1)-yo(j)
    gamma_i = 2._DBL * gamma_uo(i-1,j) * gamma_uo(i,j) / (gamma_uo(i-1,j) + gamma_uo(i,j))
    gamma_d = 2._DBL * gamma_uo(i+1,j) * gamma_uo(i,j) / (gamma_uo(i+1,j) + gamma_uo(i,j))
    gamma_s = 1._DBL/((1._DBL-fey(j)) / gamma_uo(i,j-1) + fey(j) / gamma_uo(i,j))
    gamma_n = 1._DBL/((1._DBL-fey(j+1)) / gamma_uo(i,j) + fey(j+1) / gamma_uo(i,j+1))
    delta_x = d2_xuo(i)/2._DBL
    delta_y = d_yvo(j-1)
    temp_int = (1._DBL - di/(2._DBL*delta_x)) * temp_o(i,j) + di/(2._DBL*delta_x) * temp_o(i+1,j)
!     (temp_o(i,j)+temp_o(i+1,j))/2._DBL
    !************************
    AI(i-1) =-(gamma_i*delta_y/di*DMAX1(cero,(1._DBL-0.1_DBL*dabs(ui*di/gamma_i))**5._DBL)+&
    &DMAX1(cero,ui*delta_y))
    AD(i)   =-(gamma_d*delta_y/dd*DMAX1(cero,(1._DBL-0.1_DBL*dabs(ud*dd/gamma_d))**5._DBL)+&
    &DMAX1(cero,-ud*delta_y))
    alpha   = gamma_s*delta_x/ds*DMAX1(cero,(1._DBL-0.1_DBL*dabs(vs*ds/gamma_s))**5._DBL)+&
    &DMAX1(cero, vs*delta_x)
    beta    = gamma_n*delta_x/dn*DMAX1(cero,(1._DBL-0.1_DBL*dabs(vn*dn/gamma_n))**5._DBL)+&
    &DMAX1(cero,-vn*delta_x)
    ACi(i)  = (-AI(i-1)-AD(i)+alpha+beta+delta_x*delta_y/dt_o) / rel_vo
    Rxi(i)  = alpha*u_o(i,j-1)+beta*u_o(i,j+1)-delta_x*delta_y*Ri_o(i,j)*temp_int+&
    &delta_x*delta_y*u_anto(i,j)/dt_o+(pres_o(i,j)-pres_o(i+1,j))*delta_y+&
    &ACi(i)*(1._DBL-rel_vo)*u_o(i,j)
    au_o(i,j) = ACi(i) * rel_vo
  END DO
  !***********************
  !Condiciones de frontera
  ACi(mi)  = 1._DBL
  AI(mi-1) =-1._DBL !
  Rxi(mi)  = cero
  au_o(mi,j) = 1.e40_DBL !ACi(mi)
  !*******************
  CALL tri(AI,ACi,AD,Rxi,mi)
  DO i = 1, mi
    u_o(i,j) = Rxi(i)
  END DO
END DO
!$OMP END PARALLEL DO
!********************
!TDMA en direccion y
!$OMP  PARALLEL DO DEFAULT(NONE)&
!$OMP& PRIVATE(AS,ACj,AN,Ryj,ui,ud,vs,vn,di,dd,ds,dn,gamma_i,gamma_d,gamma_s,gamma_n,delta_x,delta_y,temp_int,alpha,beta,info)&
!$OMP& SHARED(gamma_uo,au_o,u_o,u_anto,pres_o,v_o,d_xuo,d2_xuo,yo,fey,d_yvo,temp_o,rel_vo,Ri_o,dt_o)
DO i = 2, mi-1
  !***********************
  !Condiciones de frontera
  ACj(1) = 1._DBL
  AN(1)  = cero
  Ryj(1) = cero
  au_o(i,1) =1.e40_DBL !ACj(1)
  !*************************
  !Llenado de la matriz en y
  DO j = 2, nj
    !**************************
    !Interpolaciones necesarias
    ui = (u_o(i-1,j)+u_o(i,j))/2._DBL
    ud = (u_o(i,j)+u_o(i+1,j))/2._DBL
    vs = v_o(i,j-1)+d_xuo(i-1)/d2_xuo(i)*(v_o(i+1,j-1)-v_o(i,j-1))
    vn = v_o(i,j)+d_xuo(i-1)/d2_xuo(i)*(v_o(i+1,j)-v_o(i,j))
    di = d_xuo(i-1)
    dd = d_xuo(i)
    ds = yo(j)-yo(j-1)
    dn = yo(j+1)-yo(j)
    gamma_i = 2._DBL * gamma_uo(i-1,j) * gamma_uo(i,j) / (gamma_uo(i-1,j) + gamma_uo(i,j))
    gamma_d = 2._DBL * gamma_uo(i+1,j) * gamma_uo(i,j) / (gamma_uo(i+1,j) + gamma_uo(i,j))
    gamma_s = 1._DBL/((1._DBL-fey(j)) / gamma_uo(i,j-1) + fey(j) / gamma_uo(i,j))
    gamma_n = 1._DBL/((1._DBL-fey(j+1)) / gamma_uo(i,j) + fey(j+1) / gamma_uo(i,j+1))
    delta_x = d2_xuo(i)/2._DBL
    delta_y = d_yvo(j-1)
    temp_int= (1._DBL - di/(2._DBL*delta_x)) * temp_o(i,j) + di/(2._DBL*delta_x) * temp_o(i+1,j)
!     (temp_o(i,j)+temp_o(i+1,j))/2._DBL
    !************************
    alpha  = gamma_i*delta_y/di*DMAX1(cero,(1._DBL-0.1_DBL*dabs(ui*di/gamma_i))**5._DBL)+&
    &DMAX1(cero,ui*delta_y)
    beta   = gamma_d*delta_y/dd*DMAX1(cero,(1._DBL-0.1_DBL*dabs(ud*dd/gamma_d))**5._DBL)+&
    &DMAX1(cero,-ud*delta_y)
    AS(j-1)=-(gamma_s*delta_x/ds*DMAX1(cero,(1._DBL-0.1_DBL*dabs(vs*ds/gamma_s))**5._DBL)+&
    &DMAX1(cero,vs*delta_x))
    AN(j)  =-(gamma_n*delta_x/dn*DMAX1(cero,(1._DBL-0.1_DBL*dabs(vn*dn/gamma_n))**5._DBL)+&
    &DMAX1(cero,-vn*delta_x))
    ACj(j) = (alpha+beta-AS(j-1)-AN(j)+delta_x*delta_y/dt_o) / rel_vo
    Ryj(j) = alpha*u_o(i-1,j)+beta*u_o(i+1,j)-delta_x*delta_y*Ri_o(i,j)*temp_int+&
    &delta_x*delta_y*u_anto(i,j)/dt_o+(pres_o(i,j)-pres_o(i+1,j))*delta_y+&
    &ACj(j)*(1._DBL-rel_vo)*u_o(i,j)
    au_o(i,j) = ACj(j) * rel_vo
  END DO
  !***********************
  !Condiciones de frontera
  ACj(nj+1) = 1._DBL
  AS(nj)    = cero
  Ryj(nj+1) = cero
  au_o(i,nj+1) = 1.e40_DBL !ACj(nj+1)
  !*************************
  CALL tri(AS,ACj,AN,Ryj,nj+1)
  DO j = 1, nj+1
    u_o(i,j) = Ryj(j)
  END DO
END DO
!$OMP END PARALLEL DO
WHERE(u_o /= cero)
  du_o = (u_o-fu)/u_o
ELSEWHERE
  du_o = u_o-fu
END WHERE
END SUBROUTINE
