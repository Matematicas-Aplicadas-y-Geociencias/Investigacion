SUBROUTINE residuotemp(xo,yo,d_xuo,d_yvo,u_o,v_o,temp_o,temp_anto,Restemp_o,gamma_to,dt_o)
USE constantes
IMPLICIT NONE
INTEGER :: i,j,k
!*********************************
! Variables del flujo e inc'ognita
REAL(kind=DBL), DIMENSION(mi+1,nj+1), INTENT(inout) :: Restemp_o
REAL(kind=DBL), DIMENSION(mi+1,nj+1), INTENT(in)  :: temp_o,temp_anto
REAL(kind=DBL), DIMENSION(mi,nj+1),   INTENT(in)  :: u_o
REAL(kind=DBL), DIMENSION(mi+1,nj),   INTENT(in)  :: v_o
REAL(kind=DBL), DIMENSION(mi+1,nj+1) :: ftemp
!***************************
!Variables de la tridiagonal
REAL(kind=DBL) :: AC,AI,AD,AS,AN
!***************************
!Variables de interpolaci'on
REAL(kind=DBL) :: ui,ud,vn,vs,di,dd,ds,dn,delta_x,delta_y
!***************************
!Variables de la malla, volumen de control, incremento de tiempo
REAL(kind=DBL), DIMENSION(mi+1), INTENT(in) :: xo
REAL(kind=DBL), DIMENSION(nj+1), INTENT(in) :: yo
REAL(kind=DBL), DIMENSION(mi-1), INTENT(in) :: d_xuo
REAL(kind=DBL), DIMENSION(nj-1), INTENT(in) :: d_yvo
REAL(kind=DBL), INTENT(in) :: dt_o,gamma_to
!********************
!C'alculo de residuos
Restemp_o=0._DBL
DO j=2,nj
  DO i=2,mi
    !**************************
    !Interpolaciones necesarias
    ud=u_o(i,j)
    ui=u_o(i-1,j)
    vn=v_o(i,j)
    vs=v_o(i,j-1)
    di=xo(i)-xo(i-1)
    dd=xo(i+1)-xo(i)
    ds=yo(j)-yo(j-1)
    dn=yo(j+1)-yo(j)
    delta_x=d_xuo(i-1)
    delta_y=d_yvo(j-1)
    !************************
    AI=gamma_to*delta_y/di*DMAX1(0._DBL,(1._DBL-0.1_DBL*dabs(ui*di/gamma_to))**5._DBL)+&
    &DMAX1(0._DBL,ui*delta_y)
    AD=gamma_to*delta_y/dd*DMAX1(0._DBL,(1._DBL-0.1_DBL*dabs(ud*dd/gamma_to))**5._DBL)+&
    &DMAX1(0._DBL,-ud*delta_y)
    AS=gamma_to*delta_x/ds*DMAX1(0._DBL,(1._DBL-0.1_DBL*dabs(vs*ds/gamma_to))**5._DBL)+&
    &DMAX1(0._DBL,vs*delta_x)
    AN=gamma_to*delta_x/dn*DMAX1(0._DBL,(1._DBL-0.1_DBL*dabs(vn*dn/gamma_to))**5._DBL)+&
    &DMAX1(0._DBL,-vn*delta_x)
    AC=AI+AD+AN+AS+delta_x*delta_y/dt_o
    Restemp_o(i,j)=AI*temp_o(i-1,j)+AD*temp_o(i+1,j)+AS*temp_o(i,j-1)+AN*temp_o(i,j+1)+&
    &delta_x*delta_y*temp_anto(i,j)/dt_o-AC*temp_o(i,j)
    !**************************
  END DO
END DO
END SUBROUTINE
