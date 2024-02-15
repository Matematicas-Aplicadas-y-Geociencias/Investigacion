SUBROUTINE temperatura(xo,yo,fexu,feyv,d_xuo,d_yvo,u_o,v_o,temp_o,temp_anto,gamma_to,dt_o,dtemp_o,i_oo,i_1o,rel_to)
USE constantes
! USE mkl95_LAPACK
IMPLICIT NONE
INTEGER :: i,j,k,info
INTEGER, INTENT(in) :: i_oo,i_1o
!*********************************
! Variables del flujo e inc'ognita
REAL(kind=DBL), DIMENSION(mi+1,nj+1), INTENT(inout) :: temp_o
REAL(kind=DBL), DIMENSION(mi+1,nj+1), INTENT(out) :: dtemp_o
REAL(kind=DBL), DIMENSION(mi+1,nj+1), INTENT(in)  :: temp_anto,gamma_to
REAL(kind=DBL), DIMENSION(mi,nj+1),   INTENT(in)  :: u_o
REAL(kind=DBL), DIMENSION(mi+1,nj),   INTENT(in)  :: v_o
REAL(kind=DBL), DIMENSION(mi+1,nj+1) :: ftemp
!***************************
!Variables de la tridiagonal
REAL(kind=DBL), DIMENSION(mi+1) :: ACi,Rxi
REAL(kind=DBL), DIMENSION(mi)   :: AI,AD
REAL(kind=DBL), DIMENSION(nj+1) :: ACj,Ryj
REAL(kind=DBL), DIMENSION(nj)   :: AS,AN
!***************************
!Variables de interpolaci'on
REAL(kind=DBL) :: ui,ud,vn,vs,di,dd,ds,dn,gamma_i,gamma_d,gamma_s,gamma_n
REAL(kind=DBL) :: delta_x,delta_y
!***************************
!Variables de la malla, volumen de control, incremento de tiempo
REAL(kind=DBL), DIMENSION(mi+1), INTENT(in) :: xo,fexu
REAL(kind=DBL), DIMENSION(nj+1), INTENT(in) :: yo,feyv
REAL(kind=DBL), DIMENSION(mi-1), INTENT(in) :: d_xuo
REAL(kind=DBL), DIMENSION(nj-1), INTENT(in) :: d_yvo
REAL(kind=DBL), INTENT(in) :: dt_o,rel_to
!********************
!Variables auxiliares
REAL(kind=DBL) alpha, beta
!****************************
!auxiliar para calcular dtemp
ftemp = temp_o
!*******************
!TDMA en direccion x
!$OMP  PARALLEL DO DEFAULT(NONE)&
!$OMP& PRIVATE(AI,ACi,AD,Rxi,ui,ud,vs,vn,di,dd,ds,dn,gamma_i,gamma_d,gamma_s,gamma_n,delta_x,delta_y,alpha,beta,info)&
!$OMP& SHARED(u_o,v_o,temp_o,ftemp,temp_anto,d_xuo,d_yvo,xo,fexu,feyv,yo,dt_o,gamma_to,rel_to)
DO j = 2, nj
  !***********************
  !Condiciones de frontera
  ACi(1) =-1._DBL
  AD(1)  = 1._DBL
  Rxi(1) = cero
!   ACi(1) = (2._DBL*xo(1)-xo(2)-xo(3))/((xo(2)-xo(1))*(xo(3)-xo(1)))
!   AD(1)  = (xo(3)-xo(1))/((xo(3)-xo(2))*(xo(2)-xo(1)))
!   Rxi(1) = (xo(2)-xo(1))/((xo(3)-xo(2))*(xo(3)-xo(1)))*ftemp(3,j)
  !************
  ACi(mi+1) = 1._DBL
  AI(mi)    =-1._DBL
  Rxi(mi+1) = cero
!   AI(mi)    = (xo(mi-1)-xo(mi+1))/((xo(mi)-xo(mi-1))*(xo(mi+1)-xo(mi)))
!   ACi(mi+1) = (2._DBL*xo(mi+1)-xo(mi)-xo(mi-1))/((xo(mi+1)-xo(mi-1))*(xo(mi+1)-xo(mi)))
!   Rxi(mi+1) = ((xo(mi)-xo(mi+1))/((xo(mi)-xo(mi-1))*(xo(mi+1)-xo(mi-1))))*ftemp(mi-1,j)
  !*************************
  !Llenado de la matriz en x
  DO i = 2, mi
    !**************************
    !Interpolaciones necesarias
    ud = u_o(i,j)
    ui = u_o(i-1,j)
    vn = v_o(i,j)
    vs = v_o(i,j-1)
    di = xo(i)-xo(i-1)
    dd = xo(i+1)-xo(i)
    ds = yo(j)-yo(j-1)
    dn = yo(j+1)-yo(j)
    gamma_i = 1._DBL/((1._DBL-fexu(i)) / gamma_to(i-1,j) + fexu(i) / gamma_to(i,j))
    gamma_d = 1._DBL/((1._DBL-fexu(i+1)) / gamma_to(i,j) + fexu(i+1) / gamma_to(i+1,j))
    gamma_s = 1._DBL/((1._DBL-feyv(j)) / gamma_to(i,j-1) + feyv(j) / gamma_to(i,j))
    gamma_n = 1._DBL/((1._DBL-feyv(j+1)) / gamma_to(i,j) + feyv(j+1) / gamma_to(i,j+1))
!     WRITE(*,*) i,j,gamma_i,gamma_d,gamma_s,gamma_n
    delta_x = d_xuo(i-1)
    delta_y = d_yvo(j-1)
    !************************
    AI(i-1) =-(gamma_i*delta_y/di*DMAX1(cero,(1._DBL-0.1_DBL*dabs(ui*di/gamma_i))**5._DBL)+&
    &DMAX1(cero,ui*delta_y))
    AD(i)   =-(gamma_d*delta_y/dd*DMAX1(cero,(1._DBL-0.1_DBL*dabs(ud*dd/gamma_d))**5._DBL)+&
    &DMAX1(cero,-ud*delta_y))
    alpha   = gamma_s*delta_x/ds*DMAX1(cero,(1._DBL-0.1_DBL*dabs(vs*ds/gamma_s))**5._DBL)+&
    &DMAX1(cero,vs*delta_x)
    beta    = gamma_n*delta_x/dn*DMAX1(cero,(1._DBL-0.1_DBL*dabs(vn*dn/gamma_n))**5._DBL)+&
    &DMAX1(cero,-vn*delta_x)
    ACi(i)  =(-AI(i-1)-AD(i)+alpha+beta+delta_x*delta_y/dt_o) / rel_to
    Rxi(i)  = alpha*temp_o(i,j-1)+beta*temp_o(i,j+1)+delta_x*delta_y*temp_anto(i,j)/dt_o+&
    &ACi(i)*(1._DBL-rel_to)*temp_o(i,j)
    !**************************
  END DO
  CALL tri(AI,ACi,AD,Rxi,mi+1)
  DO i = 1, mi+1
    temp_o(i,j) = Rxi(i)
  END DO
END DO
!$OMP END PARALLEL DO
!*******************
!TDMA en direccion y
!$OMP  PARALLEL DO DEFAULT(NONE)&
!$OMP& PRIVATE(AS,ACj,AN,Ryj,ui,ud,vs,vn,di,dd,ds,dn,gamma_i,gamma_d,gamma_s,gamma_n,delta_x,delta_y,alpha,beta,info)&
!$OMP& SHARED(u_o,v_o,temp_o,ftemp,temp_anto,gamma_to,d_xuo,d_yvo,xo,fexu,feyv,yo,dt_o,i_oo,i_1o,rel_to)
DO i = 2, mi
!***********************
!Condiciones de frontera
  IF (i>=i_oo .and. i<=i_1o) THEN
    ACj(1) = 1._DBL
    AN(1)  = cero
    Ryj(1) = 1._DBL
  ELSE
    ACj(1) =-1._DBL
    AN(1)  = 1._DBL
    Ryj(1) = cero
!     ACj(1) = (2._DBL*yo(1)-yo(2)-yo(3))/((yo(3)-yo(1))*(yo(2)-yo(1)))
!     AN(1)  = (yo(3)-yo(1))/((yo(2)-yo(1))*(yo(3)-yo(2)))
!     Ryj(1) = ((yo(2)-yo(1))/((yo(3)-yo(2))*(yo(3)-yo(1))))*ftemp(i,3)
  END IF    
  !**********
  IF (i>=i_oo .and. i<=i_1o) THEN
    ACj(nj+1) = 1._DBL
    AS(nj)    = cero
    Ryj(nj+1) = 1._DBL
  ELSE
    AS(nj)    =-1._DBL
    ACj(nj+1) = 1._DBL
    Ryj(nj+1) = cero
!     AS(nj)    = (yo(nj-1)-yo(nj+1))/((yo(nj)-yo(nj-1))*(yo(nj+1)-yo(nj)))
!     ACj(nj+1) = (2._DBL*yo(nj+1)-yo(nj)-yo(nj-1))/((yo(nj+1)-yo(nj-1))*(yo(nj+1)-yo(nj)))
!     Ryj(nj+1) = ((yo(nj)-yo(nj+1))/((yo(nj)-yo(nj-1))*(yo(nj+1)-yo(nj-1))))*ftemp(i,nj-1)
  END IF 
  !*************************
  !Llenado de la matriz en x
  DO j=2,nj
    !**************************
    !Interpolaciones necesarias
    ud = u_o(i,j)
    ui = u_o(i-1,j)
    vn = v_o(i,j)
    vs = v_o(i,j-1)
    di = xo(i)-xo(i-1)
    dd = xo(i+1)-xo(i)
    ds = yo(j)-yo(j-1)
    dn = yo(j+1)-yo(j)
    gamma_i = 1._DBL/((1._DBL-fexu(i)) / gamma_to(i-1,j) + fexu(i) / gamma_to(i,j))
    gamma_d = 1._DBL/((1._DBL-fexu(i+1)) / gamma_to(i,j) + fexu(i+1) / gamma_to(i+1,j))
    gamma_s = 1._DBL/((1._DBL-feyv(j)) / gamma_to(i,j-1) + feyv(j) / gamma_to(i,j))
    gamma_n = 1._DBL/((1._DBL-feyv(j+1)) / gamma_to(i,j) + feyv(j+1) / gamma_to(i,j+1))
    delta_x = d_xuo(i-1)
    delta_y = d_yvo(j-1)
    !************************
    alpha   = gamma_i*delta_y/di*DMAX1(cero,(1._DBL-0.1_DBL*dabs(ui*di/gamma_i))**5._DBL)+&
    &DMAX1(cero,ui*delta_y)
    beta    = gamma_d*delta_y/dd*DMAX1(cero,(1._DBL-0.1_DBL*dabs(ud*dd/gamma_d))**5._DBL)+&
    &DMAX1(cero,-ud*delta_y)
    AS(j-1) =-(gamma_s*delta_x/ds*DMAX1(cero,(1._DBL-0.1_DBL*dabs(vs*ds/gamma_s))**5._DBL)+&
    &DMAX1(cero,vs*delta_x))
    AN(j)   =-(gamma_n*delta_x/dn*DMAX1(cero,(1._DBL-0.1_DBL*dabs(vn*dn/gamma_n))**5._DBL)+&
    &DMAX1(cero,-vn*delta_x))
    ACj(j)  = (alpha+beta-AS(j-1)-AN(j)+delta_x*delta_y/dt_o) / rel_to
    Ryj(j)  = alpha*temp_o(i-1,j)+beta*temp_o(i+1,j)+delta_x*delta_y*temp_anto(i,j)/dt_o+&
    &ACj(j)*(1._DBL-rel_to)*temp_o(i,j)
    !**************************
  END DO
  CALL tri(AS,ACj,AN,Ryj,nj+1)
  DO j = 1, nj+1
    temp_o(i,j) = Ryj(j)
  END DO
END DO
!$OMP END PARALLEL DO
WHERE(temp_o /= cero)
  dtemp_o = (ftemp-temp_o)/temp_o
ELSEWHERE
  dtemp_o = ftemp-temp_o
ENDWHERE
END SUBROUTINE
